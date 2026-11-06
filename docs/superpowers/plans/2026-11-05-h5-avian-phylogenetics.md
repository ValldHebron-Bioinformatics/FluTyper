# H5 Avian Phylogenetics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend FluTyper's phylogenetics/Auspice module to the AVIAN protocol and the `H5` subtype, then use it to test, segment by segment, whether 25CR6244 (H5 LPAIV) has reassorted with HPAIV Gs/GD clade 2.3.4.4b viruses.

**Architecture:** The module already partitions samples by subtype and builds one tree per segment per subtype. AVIAN needs a subtype key that matches on the **H tag alone** (`H5`), because the LPAI context pools H5N1/H5N2/H5N3/H5N4 — grouping by `h_tag+n_tag` would split it. Everything else is reached by extending the existing per-subtype registry (`phyloKnownSubtypes`/`phyloSubtypeTag`/`phyloSubtypeDefaults`/`resolvePhyloConfig`), adding one `if ('H5' in phylo_subtypes)` block in `main.nf` mirroring the two that exist, and making the tip annotation fields declarative so avian fields (pathotype, clade, genotype, host) reach both the figure panels and the Auspice colourings. The three datasets are three pipeline runs with different `--inputFasta`/`--outDir`.

**Tech Stack:** Nextflow DSL2, nf-test, Python 3 (pandas, BioPython), MAFFT 7.525, IQ-TREE 2.3.6, TreeTime 0.11.4, TreeCluster 1.0.4, augur/Auspice (envs/nextstrain.yaml), R/ggtree.

**Spec:** `docs/superpowers/specs/2026-11-05-h5-avian-reassortment-25CR6244.md`

## Global Constraints

- Protocol resources for AVIAN live in `protocols/AVIAN/v1`; the whole-genome reference file there is named **`Avian_references.fasta`**, not `AVIAN_references.fasta`.
- The only AVIAN reference strain with all 8 segments is `A/American_wigeon/South_Dakota/23-033792-006-original/2023` (H5N1, HPAI, clade 2.3.4.4b). It is the alignment coordinate reference, the sole root strain and the outgroup for `H5`.
- AVIAN reference headers are `{SUBTYPE}_{SEGMENT}_{STRAIN}_{EPI_ISL}_{PATHOTYPE}`; AVIAN CDS reference headers are `{H5|N1|H7|…}_{GENE}_REFERENCE_CDS_{ACC}_{PATHOTYPE}` — the CDS prefix is the **H tag or N tag**, never `H5N1`.
- Segments analysed for 25CR6244: `HA,PB2,PB1,NP,MP,NS` (it has no PA and no NA). `--phyloSegmentOrder` keeps all 8; `--phyloSegmentTrees` selects the 6.
- `--phyloWholeGenome` stays `false`: segment trees are what answer a reassortment question.
- No molecular clock rate, date or coordinate is ever invented. A subtype/segment with no row in `RESOURCES/auspice_clock_rates.tsv` must keep the module's existing behaviour: divergence-only export plus a WARNING.
- `phyloKnownSubtypes()` in `nf_pipeline/subworkflows/Phylogenetics.nf` and `PHYLO_KNOWN_SUBTYPES` in `nf_pipeline/modules/MergeReports.nf` must stay in sync, and no entry may be a prefix of another entry's `_SEGMENT` suffix.
- Python helper scripts live in `scripts/` and are tested with stdlib `unittest` (`python3 -m unittest discover -s tests/scripts`); Nextflow code is tested with nf-test.

## Review Focus

Input classes the spec implies that no task's happy path exercises. Each has its test added to the owning task.

1. **A GISAID export whose `Pathogenicity` cell is empty** (9 of 1369 rows). The screening must classify them from sequence and the dataset builder must never silently drop them — Task 1, Task 2.
2. **A sample with fewer than 8 segments** (25CR6244 has 6). It must stay in the 6 trees it has data for and be excluded only from PA and NA, not from everything — Task 9.
3. **An `H5` sample whose N tag could not be called** (`Nx`; 2 isolates are bare `A / H5`). H-tag-only matching must still include it — Task 4.
4. **A metadata row whose `LOCATION` does not geocode.** The location must be kept as written with a WARNING, never replaced by a guess or dropped — Task 1.
5. **An `H5` tree with no clock rate row.** Export must fall back to divergence-only with a WARNING, not fail and not guess a rate — Task 8.

---

### Task 1: GISAID export → FluTyper metadata, with avian fields

`scripts/gisaid_to_metadata.py` already converts a GISAID EpiFlu export into `--metadata` + `--coordinates`. It cannot read this export: it strips `EPI_ISL_` from `Isolate_Id` while this download's FASTA headers keep it, so `g[g["Isolate_Id"].isin(labs)]` selects zero rows. It also drops the avian columns the rest of the plan needs.

**Files:**
- Modify: `scripts/gisaid_to_metadata.py:93-120` (`main`), and add the parsing helpers below
- Create: `tests/scripts/test_gisaid_to_metadata.py`
- Create: `tests/data/gisaid_avian_sample.csv` (12 synthetic rows with the real column names, covering: HPAI 2.3.4.4b, LPAI H5N3, blank `Pathogenicity`, a `Location` with no city level, and a `Location` that cannot geocode)

**Interfaces:**
- Produces: `isolate_key(raw_id: str) -> str` returning the `EPI_ISL_`-prefixed form of either spelling; `read_export(path: str) -> pandas.DataFrame` accepting `.xls`, `.xlsx` and `.csv`; `main()` writing a metadata CSV with columns `ID,DATE,LOCATION,COUNTRY,HOST,SUBTYPE,PATHOGENICITY,CLADE,GENOTYPE,AGE GROUP,SEX,ORIGINATING_LAB`.
- Consumes: nothing.

- [ ] **Step 1: Write the failing tests**

```python
class TestIsolateKey(unittest.TestCase):
    def test_accepts_both_spellings(self):
        self.assertEqual(isolate_key("EPI_ISL_20298047"), "EPI_ISL_20298047")
        self.assertEqual(isolate_key("20298047"), "EPI_ISL_20298047")

class TestMain(unittest.TestCase):
    def test_keeps_rows_whose_fasta_ids_carry_the_epi_isl_prefix(self):
        meta = pandas.read_csv(self.out_csv)
        self.assertEqual(len(meta), 12)
        self.assertEqual(list(meta.columns), ["ID", "DATE", "LOCATION", "COUNTRY", "HOST",
                                              "SUBTYPE", "PATHOGENICITY", "CLADE", "GENOTYPE",
                                              "AGE GROUP", "SEX", "ORIGINATING_LAB"])

    def test_blank_pathogenicity_becomes_empty_not_dropped(self):
        meta = pandas.read_csv(self.out_csv, keep_default_na=False)
        self.assertEqual(len(meta[meta["PATHOGENICITY"] == ""]), 1)

    def test_ungeocodable_location_is_kept_as_written_and_warned(self):
        meta = pandas.read_csv(self.out_csv)
        self.assertIn("Nowhereville", list(meta["LOCATION"]))
        self.assertIn("could not geocode", self.stderr)
```

Geocoding is network I/O: the test monkeypatches `gisaid_to_metadata._nominatim` with a dict-backed stub that returns `None` for `Nowhereville`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: FAIL — `NameError: name 'isolate_key' is not defined`

- [ ] **Step 3: Implement `isolate_key`, `read_export` and the new columns in `scripts/gisaid_to_metadata.py`**

`read_export` dispatches on the suffix: `.csv` → `pandas.read_csv`, otherwise `pandas.read_excel`. Key both `labs` and `Isolate_Id` through `isolate_key` before the `isin` filter. `parse_location` already returns `(place, region, country)`: emit `country` as the new `COUNTRY` column so Auspice's geo join has it without a second lookup. `LOCATION` must keep the deepest level even when it did not geocode — the existing `for place in reversed(...)` loop falls through to `""` and has to return the deepest non-empty level instead.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: PASS

- [ ] **Step 5: Add `xlrd >=2.0.1` to `FluTyper_env.yaml`**

The export is a real legacy `.xls` (Composite Document File V2); pandas needs `xlrd` for it.

- [ ] **Step 6: Commit**

```bash
git add scripts/gisaid_to_metadata.py tests/scripts/test_gisaid_to_metadata.py tests/data/gisaid_avian_sample.csv FluTyper_env.yaml
git commit -m "feat(scripts): read EPI_ISL-prefixed GISAID exports and carry avian fields"
```

---

### Task 2: Dataset builder

**Files:**
- Create: `scripts/build_phylo_datasets.py`
- Create: `tests/scripts/test_build_phylo_datasets.py`

**Interfaces:**
- Consumes: Task 1's metadata CSV; Task 3's pathotype table (`ID,inferred_subtype,pathotype`); the two input FASTAs.
- Produces: `classify(meta_row: dict, called: dict[str, str]) -> str` returning `"LPAI"`, `"HPAI"` or `"UNKNOWN"`, where `called` is FluTyper's own pathotype per ID and **wins over** the GISAID field; `build(dataset: str, ...) -> dict` writing `<outdir>/<dataset>.fasta` and `<outdir>/<dataset>_metadata.csv` and returning the per-dataset counts.
- CLI: `build_phylo_datasets.py <gisaid_meta.csv> <pathotype.csv> <gisaid.fasta> <own.fasta> <own_meta.csv> <focus_id> <outdir>`

- [ ] **Step 1: Write the failing tests**

```python
class TestClassify(unittest.TestCase):
    def test_own_pathotype_call_wins_over_the_gisaid_field(self):
        row = {"ID": "EPI_ISL_1", "PATHOGENICITY": "HPAI", "CLADE": "2.3.4.4b"}
        self.assertEqual(classify(row, {"EPI_ISL_1": "LPAI"}), "LPAI")

    def test_blank_gisaid_field_is_resolved_by_the_call(self):
        row = {"ID": "EPI_ISL_2", "PATHOGENICITY": "", "CLADE": ""}
        self.assertEqual(classify(row, {"EPI_ISL_2": "LPAI"}), "LPAI")

    def test_no_call_and_no_field_is_unknown_never_assumed(self):
        row = {"ID": "EPI_ISL_3", "PATHOGENICITY": "", "CLADE": ""}
        self.assertEqual(classify(row, {}), "UNKNOWN")

class TestBuild(unittest.TestCase):
    def test_focus_sample_is_in_every_dataset(self):
        for ds in ("lpai", "hpai", "combined"):
            ids = {h.split("_")[0] for h in fasta_headers(self.outdir / f"{ds}.fasta")}
            self.assertIn("25CR6244", ids)

    def test_hpai_dataset_is_every_2344b_isolate(self):
        self.assertEqual(self.counts["hpai"]["context"], self.n_clade_2344b)

    def test_combined_is_the_union_of_both_contexts(self):
        self.assertEqual(self.counts["combined"]["context"],
                         self.counts["lpai"]["context"] + self.counts["hpai"]["context"])

    def test_unknown_isolates_are_listed_not_silently_dropped(self):
        self.assertIn("UNKNOWN", self.report)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: FAIL — module `build_phylo_datasets` not found

- [ ] **Step 3: Implement `scripts/build_phylo_datasets.py`**

Three datasets: `lpai` (every `LPAI`), `hpai` (every isolate whose `CLADE` is `2.3.4.4b`), `combined` (the union). Each one also gets `focus_id`'s records from the own-sample FASTA and its row from the own metadata, with `PATHOGENICITY` set to its own called pathotype. Records are selected from the GISAID FASTA by `isolate_key` of the first `|` field, so each header's segment and lab fields are copied through untouched. `UNKNOWN` isolates go into no dataset and are listed in `<outdir>/dataset_report.txt` alongside the per-dataset counts.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/build_phylo_datasets.py tests/scripts/test_build_phylo_datasets.py
git commit -m "feat(scripts): build LPAI/HPAI/combined H5 context datasets"
```

---

### Task 3: Screening run — pathotype from sequence, compared with GISAID

The spec's first deliverable, and Task 2's input. No pipeline code changes.

**Files:**
- Create: `scripts/compare_pathotype.py`
- Create: `tests/scripts/test_compare_pathotype.py`
- Create: `docs/superpowers/results/2026-11-05-screening.md`
- Produces (not committed): `RESULTS/screening/`

**Interfaces:**
- Produces: `compare(called: pandas.DataFrame, gisaid: pandas.DataFrame) -> pandas.DataFrame` with columns `ID,inferred_subtype,pathotype_called,pathotype_gisaid,clade,genotype,agreement`, where `agreement` is one of `agree`, `disagree`, `only_called`, `only_gisaid`, `neither`.

- [ ] **Step 1: Write the failing test**

```python
def test_agreement_classes(self):
    out = compare(called_df, gisaid_df).set_index("ID")["agreement"].to_dict()
    self.assertEqual(out, {"EPI_ISL_1": "agree", "EPI_ISL_2": "disagree",
                           "EPI_ISL_3": "only_called", "EPI_ISL_4": "only_gisaid",
                           "EPI_ISL_5": "neither"})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: FAIL — module `compare_pathotype` not found

- [ ] **Step 3: Implement `compare` in `scripts/compare_pathotype.py`**

- [ ] **Step 4: Run the test to verify it passes**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: PASS

- [ ] **Step 5: Build the screening input and metadata**

```bash
cat docs/fastas/AVIAN_JUNY_NOV_2025.fasta docs/fastas/cresaid.fasta > /tmp/screening_input.fasta
python3 scripts/gisaid_to_metadata.py docs/AVIAN_JUNY_NOV_2025.xls \
    docs/fastas/AVIAN_JUNY_NOV_2025.fasta /tmp/gisaid_meta.csv /tmp/gisaid_coords.tsv
```

Then append `docs/metadata.csv`'s rows to `/tmp/gisaid_meta.csv` under the same header, leaving the avian columns empty for the own samples.

- [ ] **Step 6: Run the screening**

```bash
nextflow run nf_pipeline/main.nf --protocol AVIAN \
    --inputFasta /tmp/screening_input.fasta --outDir RESULTS/screening \
    --metadata /tmp/gisaid_meta.csv --coordinates /tmp/gisaid_coords.tsv
```

Expected: `RESULTS/screening/final_genotyping_results.csv` has one row per isolate — 1369 GISAID plus 21 own.

- [ ] **Step 7: Compare and write up `docs/superpowers/results/2026-11-05-screening.md`**

Run `scripts/compare_pathotype.py` over the genotyping results and `/tmp/gisaid_meta.csv`.
Expected: the 9 blank-`Pathogenicity` isolates appear as `only_called` or `neither`, and every `disagree` row is listed with its HA cleavage site. Record 25CR6244's own called pathotype and subtype, and state that N1 is unconfirmed because the NA segment is absent.

- [ ] **Step 8: Commit**

```bash
git add scripts/compare_pathotype.py tests/scripts/test_compare_pathotype.py docs/superpowers/results/2026-11-05-screening.md
git commit -m "feat(scripts): compare called and GISAID pathotypes; add H5 screening results"
```

---

### Task 4: `H5` subtype registry and AVIAN protocol gate

**Files:**
- Modify: `nf_pipeline/subworkflows/Phylogenetics.nf` — the subtype filters at `:63` and `:80`, `validatePhyloParams` (`:622`), `phyloKnownSubtypes` (`:807`), `phyloSubtypeDefaults` (`:813`), `phyloSubtypeTag` (`:826`)
- Modify: `nf_pipeline/modules/MergeReports.nf:26` (`PHYLO_KNOWN_SUBTYPES`)
- Modify: `nf_pipeline/nextflow.config` — add `phyloSubtypeConfig.H5` and `aaPositionsH5`
- Test: `tests/modules/phylogeny/PhyloFunctions.nf.test`, `tests/modules/MergeReports.nf.test`

**Interfaces:**
- Produces: `phyloKnownSubtypes()` returns `['H3N2', 'H1N1pdm09', 'H5']`; `phyloSubtypeTag('H5')` returns `'H5'`; `phyloSubtypeMatches(String tag, String h_tag, String n_tag) -> boolean`, true when `tag` is an H-only tag equal to `h_tag`, or when `tag == h_tag + n_tag`; `resolvePhyloConfig('H5', p)` returns `A/American_wigeon/South_Dakota/23-033792-006-original/2023` as `reference`, `referenceFallback` and `outgroup`, and as the only entry of `rootStrains`.
- Consumes: nothing.

- [ ] **Step 1: Write the failing tests**

```groovy
test("phyloSubtypeTag: H5 matches on the H tag alone") {
    function "phyloSubtypeTag"
    when { function { """input[0] = "H5" """ } }
    then { assert function.success; assert function.result == "H5" }
}

test("phyloSubtypeMatches: H5 includes every N tag") {
    function "phyloSubtypeMatches"
    when { function { """input[0] = "H5"; input[1] = "H5"; input[2] = "N3" """ } }
    then { assert function.success; assert function.result == true }
}

test("phyloSubtypeMatches: H5 includes a sample whose N tag is Nx") {
    function "phyloSubtypeMatches"
    when { function { """input[0] = "H5"; input[1] = "H5"; input[2] = "Nx" """ } }
    then { assert function.success; assert function.result == true }
}

test("phyloSubtypeMatches: H3N2 still requires both tags") {
    function "phyloSubtypeMatches"
    when { function { """input[0] = "H3N2"; input[1] = "H3"; input[2] = "N1" """ } }
    then { assert function.success; assert function.result == false }
}

test("validatePhyloParams: AVIAN with H5 is accepted") {
    function "validatePhyloParams"
    then { assert function.success; assert function.result == [] }
}

test("validatePhyloParams: AVIAN with H3N2 is rejected") {
    then { assert function.result.any { it.contains("H3N2") && it.contains("AVIAN") } }
}

test("validatePhyloParams: HUMAN with H5 is rejected") {
    then { assert function.result.any { it.contains("H5") && it.contains("HUMAN") } }
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `nf-test test tests/modules/phylogeny/PhyloFunctions.nf.test --verbose`
Expected: FAIL — `Unknown phylogenetics subtype 'H5'`, and `phyloSubtypeMatches` is not defined

- [ ] **Step 3: Implement the registry changes**

Replace the single protocol check in `validatePhyloParams` with a per-subtype protocol map (`H3N2`, `H1N1pdm09` → `HUMAN`; `H5` → `AVIAN`), so each error names both the subtype and the protocol it needs and a mixed `--phyloSubtype` is rejected. Add `phyloSubtypeMatches` and call it at both filter sites in place of `"${h_tag}${n_tag}" == subtype_tag`. `phyloSubtypeDefaults()['H5']` and `nextflow.config`'s `phyloSubtypeConfig.H5` both carry the American wigeon strain. `aaPositionsH5` defaults to `HA1:160,PB2:627,PB2:701` — the T160A glycosylation site and the two canonical PB2 mammalian-adaptation positions.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `nf-test test tests/modules/phylogeny/PhyloFunctions.nf.test tests/modules/MergeReports.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add nf_pipeline/subworkflows/Phylogenetics.nf nf_pipeline/modules/MergeReports.nf nf_pipeline/nextflow.config tests/modules/phylogeny/PhyloFunctions.nf.test tests/modules/MergeReports.nf.test
git commit -m "feat(phylo): register the H5 subtype and allow the AVIAN protocol"
```

---

### Task 5: AVIAN references and CDS prefixes in preflight

**Files:**
- Modify: `nf_pipeline/subworkflows/Phylogenetics.nf:38` (reference genome filename)
- Modify: `nf_pipeline/modules/phylogeny/PhyloPreflight.nf:123` (`ref_subtype`) and `:181-185` (CDS prefix filter)
- Test: `tests/modules/phylogeny/PhyloPreflight.nf.test`
- Create: `tests/data/phylo/avian_references.fasta`, `tests/data/phylo/avian_cds_references.fasta` (the 8 segments and the genes of the American wigeon strain, plus one decoy `H7_HA_…` record)

**Interfaces:**
- Consumes: `phyloSubtypeDefaults()['H5']` from Task 4.
- Produces: `PhyloPreflight` emitting `reference_segments.tsv` with all 8 segments resolved to the American wigeon strain, and a codon map covering `HA1`, `HA2`, `NA`, `PB2`, `PB1`, `PB1-F2`, `PA`, `PA-X`, `NP`, `M1`, `M2`, `NS1`, `NS2` for `H5`.

- [ ] **Step 1: Write the failing test**

```groovy
test("PhyloPreflight: AVIAN H5 resolves all 8 segments and maps the H5/N1 CDS prefixes") {
    then {
        assert process.success
        def segs = path(process.out.reference_table[0]).readLines()
        assert segs.size() == 9  // header + 8 segments
        assert segs.drop(1).every { it.contains("American_wigeon") }
        def codons = path(process.out.coordinate_table[0]).text
        assert codons.contains("HA1") && codons.contains("NA") && codons.contains("PA-X")
        assert !path(process.out.log[0]).text.contains("WARNING")
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `nf-test test tests/modules/phylogeny/PhyloPreflight.nf.test --verbose`
Expected: FAIL — the codon map is empty, because no CDS record's first `_` field equals `H5N1` or `H5`

- [ ] **Step 3: Implement the reference and CDS prefix resolution**

In `Phylogenetics.nf`, resolve the reference genome by trying `${params.protocol}_references.fasta` and then the capitalised spelling (`Avian_references.fasta`), failing with both names when neither exists — renaming the shipped file would break the other modules that read it. In `PhyloPreflight`, derive the accepted CDS prefixes from the subtype: `H1N1pdm09` → `{"H1N1pdm09", "H1N1"}` as today, and an H-only subtype such as `H5` → `{"H5"}` plus the N tags that appear in that strain's reference headers, so `N1_NA_…` is accepted for the NA gene. Keep `parts[0] in accepted` as the membership test.

- [ ] **Step 4: Run the test to verify it passes**

Run: `nf-test test tests/modules/phylogeny/PhyloPreflight.nf.test`
Expected: PASS

- [ ] **Step 5: Run the HUMAN preflight and alignment tests to confirm no regression**

Run: `nf-test test tests/modules/phylogeny/PhyloPreflight.nf.test tests/modules/phylogeny/PhyloAlignment.nf.test`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add nf_pipeline/subworkflows/Phylogenetics.nf nf_pipeline/modules/phylogeny/PhyloPreflight.nf tests/modules/phylogeny/PhyloPreflight.nf.test tests/data/phylo/avian_references.fasta tests/data/phylo/avian_cds_references.fasta
git commit -m "feat(phylo): resolve AVIAN references and H5/N1 CDS prefixes in preflight"
```

---

### Task 6: Wire `H5` into the workflow

**Files:**
- Modify: `nf_pipeline/main.nf:70-105` (aa-position resolution), `:425-470` (subtype blocks), `:548-549` (publish declarations), `:636-646` (output targets)
- Test: `tests/modules/phylogeny/Phylogenetics.nf.test`

**Interfaces:**
- Consumes: `phyloSubtypeTag`, `resolvePhyloConfig`, `phyloKnownSubtypes` from Task 4.
- Produces: `PHYLOGENETICS_H5`, instantiated exactly like `PHYLOGENETICS_H3N2`, with its `outputs` published to `<outDir>/phylogeny/H5` through a `phylogeny_h5` target.

- [ ] **Step 1: Write the failing test**

```groovy
test("PHYLOGENETICS: an H5 subworkflow run builds one tree per requested segment") {
    then {
        assert workflow.success
        assert workflow.out.outputs.any { it.contains("phylogeny/H5") }
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `nf-test test tests/modules/phylogeny/Phylogenetics.nf.test --verbose`
Expected: FAIL — no `H5` branch exists

- [ ] **Step 3: Add the `H5` block, publish declaration and output target**

Mirror the `if ('H3N2' in phylo_subtypes)` block: a `ch_phylogeny_h5` channel, `phylo_positions_h5 = parseAaPositions(params.aaPositionsH5, params.aaGene)` with the same `--aa-positions-file` override, mixing into `ch_phylo_errors`, `ch_phylo_report`, `ch_phylo_tree_problems` and `ch_auspice`, and an output target `phylogeny_h5 { path { "${projectDir}/../${params.outDir}/phylogeny/H5" }; mode "copy" }`.

- [ ] **Step 4: Run the test to verify it passes**

Run: `nf-test test tests/modules/phylogeny/Phylogenetics.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add nf_pipeline/main.nf tests/modules/phylogeny/Phylogenetics.nf.test
git commit -m "feat(phylo): build and publish H5 segment trees"
```

---

### Task 7: Avian tip annotations in the figure and in Auspice

`PhyloReport` builds `annotations.tsv` from the declarative `RESOURCES/phylo_panels.tsv`, but `PhyloAuspice` hardcodes `ANN_FIELDS = ["subclade", "cluster", "season", "age_group", "sex"]`, so a new column reaches the figure and never the Auspice colourings.

**Files:**
- Modify: `nf_pipeline/modules/phylogeny/PhyloAuspice.nf:113` (`ANN_FIELDS`)
- Create: `RESOURCES/phylo_panels_avian.tsv`, `RESOURCES/auspice_config_avian.json`
- Modify: `nf_pipeline/nextflow.config` (`phyloPanels`, `auspiceConfig` defaults per protocol)
- Test: `tests/modules/phylogeny/PhyloAuspice.nf.test`, `tests/modules/phylogeny/AuspiceConfig.nf.test`

**Interfaces:**
- Produces: `ANN_FIELDS` read from the staged panel spec's `column` values (minus the `@aa_positions` placeholder) instead of a literal list, so `pathotype`, `clade`, `genotype` and `host` become colourings, filters and valid `--auspiceTraits` columns.

- [ ] **Step 1: Write the failing test**

```groovy
test("PhyloAuspice: avian annotation columns become colourings and filters") {
    then {
        assert process.success
        def ds = new JsonSlurper().parse(file(process.out.dataset[0]))
        def keys = ds.meta.colorings.collect { it.key }
        assert keys.containsAll(["pathotype", "clade", "genotype", "host"])
        assert ds.meta.filters.containsAll(["pathotype", "clade", "genotype", "host"])
        assert ds.meta.geo_resolutions.collect { it.key }.contains("country")
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `nf-test test tests/modules/phylogeny/PhyloAuspice.nf.test --verbose`
Expected: FAIL — the colourings contain `age_group` and `sex` and none of the avian keys

- [ ] **Step 3: Derive `ANN_FIELDS` from the panel spec and add the avian resources**

`phylo_panels_avian.tsv` panels: `subclade` (required), `pathotype`, `clade`, `genotype`, `host`, `period`, `aa` (required). `auspice_config_avian.json` mirrors the existing config with `pathotype`, `clade`, `genotype` and `host` in place of `age_group`/`sex`, `color_by` defaulting to `pathotype`, and `geo_resolutions` `country`/`division`/`location`. `nextflow.config` selects the avian spec and config when `params.protocol == 'AVIAN'`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `nf-test test tests/modules/phylogeny/PhyloAuspice.nf.test tests/modules/phylogeny/AuspiceConfig.nf.test tests/modules/phylogeny/PhyloReport.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add nf_pipeline/modules/phylogeny/PhyloAuspice.nf nf_pipeline/nextflow.config RESOURCES/phylo_panels_avian.tsv RESOURCES/auspice_config_avian.json tests/modules/phylogeny/PhyloAuspice.nf.test tests/modules/phylogeny/AuspiceConfig.nf.test
git commit -m "feat(auspice): drive tip annotations from the panel spec and add avian resources"
```

---

### Task 8: H5 clock rates, or an honest divergence-only fallback

**Files:**
- Modify: `RESOURCES/auspice_clock_rates.tsv`
- Test: `tests/modules/phylogeny/PhyloAuspice.nf.test`

**Interfaces:**
- Consumes: nothing.
- Produces: either 8 `H5` rows each with a citable `source`, or no `H5` rows at all.

- [ ] **Step 1: Write the failing test**

```groovy
test("PhyloAuspice: an H5 tree with no clock rate exports divergence-only and warns") {
    then {
        assert process.success
        def log = path(process.out.step_log[0]).text
        assert log.contains("WARNING") && log.contains("divergence")
        def ds = new JsonSlurper().parse(file(process.out.dataset[0]))
        assert !ds.meta.panels.contains("frequencies")
    }
}
```

- [ ] **Step 2: Run the test**

Run: `nf-test test tests/modules/phylogeny/PhyloAuspice.nf.test --verbose`
Expected: PASS is acceptable here — the module documents this fallback. A FAIL means the documented behaviour is not real and must be fixed before any rate is added.

- [ ] **Step 3: Look up per-segment H5 clade 2.3.4.4b clock rates in `nextstrain/avian-flu`**

Read the `clock_rate`/`clock_std_dev` values in that repository's workflow config at a pinned commit. Add one row per segment with `generic = no` and the repository, file and commit in `source`, exactly as the existing H3N2 rows cite `nextstrain/seasonal-flu core.smk@a1e704c7`. **If the values cannot be found at a pinned, citable location, add no rows** and record in the commit message that H5 exports are divergence-only by design.

- [ ] **Step 4: Run the test to confirm the chosen behaviour**

Run: `nf-test test tests/modules/phylogeny/PhyloAuspice.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add RESOURCES/auspice_clock_rates.tsv tests/modules/phylogeny/PhyloAuspice.nf.test
git commit -m "feat(auspice): H5 clock rates (or documented divergence-only export)"
```

---

### Task 9: End-to-end AVIAN phylogenetics test

**Files:**
- Create: `tests/phylogeny_avian.nf.test`
- Create: `tests/data/AVIAN_phylo.fasta` (6 synthetic H5 samples: 4 H5N1 HPAI, 1 H5N3 LPAI, and `SIXSEG`, a six-segment sample standing in for 25CR6244), `tests/data/avianmetadata.csv`

**Interfaces:**
- Consumes: everything from Tasks 4–8.

- [ ] **Step 1: Write the failing test**

```groovy
test("AVIAN H5: segment trees are built, and the six-segment sample is in six of them") {
    when { params {
        protocol = "AVIAN"; phylogenetics = true; phyloSubtype = "H5"
        phyloSegmentTrees = "HA,PB2,PB1,NP,MP,NS"; phyloWholeGenome = false
    } }
    then {
        assert workflow.success
        ["HA","PB2","PB1","NP","MP","NS"].each { seg ->
            assert path("$outputDir/phylogeny/H5/$seg/tree.nwk").text.contains("SIXSEG")
        }
        assert !outputDir.resolve("phylogeny/H5").toFile().list().toList().contains("PA")
        def qc = path("$outputDir/phylo_coverage_qc.tsv").readLines()
        assert qc.findAll { it.startsWith("SIXSEG") && it.endsWith("PASS") }.size() == 6
        PhyloTestHelpers.checkQuietTerminal(workflow.stdout)
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `nf-test test tests/phylogeny_avian.nf.test --verbose`
Expected: FAIL until the fixtures exist and every earlier task is in

- [ ] **Step 3: Build the fixtures from the real references**

Derive the synthetic samples from `protocols/AVIAN/v1/Avian_references.fasta`'s American wigeon segments with a handful of substitutions each, so the alignment and the codon map are real while the files stay small.

- [ ] **Step 4: Run the whole suite**

Run: `nf-test test tests/modules/*.nf.test tests/modules/phylogeny/*.nf.test tests/phylogeny.nf.test tests/phylogeny_avian.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/phylogeny_avian.nf.test tests/data/AVIAN_phylo.fasta tests/data/avianmetadata.csv
git commit -m "test(phylo): end-to-end AVIAN H5 segment tree run"
```

---

### Task 10: Build the three datasets and run the LPAI context

**Files:**
- Produces (not committed): `/tmp/h5_datasets/`, `RESULTS/h5_lpai/`

- [ ] **Step 1: Build the datasets**

```bash
python3 scripts/build_phylo_datasets.py /tmp/gisaid_meta.csv \
    RESULTS/screening/final_genotyping_results.csv \
    docs/fastas/AVIAN_JUNY_NOV_2025.fasta docs/fastas/cresaid.fasta docs/metadata.csv \
    25CR6244 /tmp/h5_datasets
```

Expected: `dataset_report.txt` showing roughly 34 LPAI, 1328 HPAI and 1362 combined, each including 25CR6244, with every `UNKNOWN` isolate listed.

- [ ] **Step 2: Run the LPAI context**

```bash
nextflow run nf_pipeline/main.nf --protocol AVIAN \
    --inputFasta /tmp/h5_datasets/lpai.fasta --outDir RESULTS/h5_lpai \
    --metadata /tmp/h5_datasets/lpai_metadata.csv --coordinates /tmp/gisaid_coords.tsv \
    --phylogenetics true --auspice true --phyloSubtype H5 \
    --phyloSegmentTrees HA,PB2,PB1,NP,MP,NS --bootstrap 1000
```

Expected: `RESULTS/h5_lpai/phylogeny/H5/<SEG>/tree.nwk` for the 6 segments, 6 Auspice datasets in `RESULTS/h5_lpai/auspice/`, and no `WARNING` in `phylo.log` other than the documented clock-rate one.

- [ ] **Step 3: Check 25CR6244 is a tip of all six trees**

Expected: each `tree.nwk` contains `25CR6244`, and `phylo_coverage_qc.tsv` shows it PASS for exactly those 6 segments.

---

### Task 11: Run the HPAI and combined contexts

**Files:**
- Produces (not committed): `RESULTS/h5_hpai/`, `RESULTS/h5_combined/`

- [ ] **Step 1: Run the HPAI context**

The Task 10 Step 2 command with `--inputFasta /tmp/h5_datasets/hpai.fasta --outDir RESULTS/h5_hpai`. About 1328 tips per segment; `executor.cpus` is 18 and `big_mem` is capped at 12 CPUs / 10 GB in `nextflow.config`, so expect hours per segment and run it in the background.

- [ ] **Step 2: Run the combined context**

The same command with `--inputFasta /tmp/h5_datasets/combined.fasta --outDir RESULTS/h5_combined`.

- [ ] **Step 3: Verify both runs**

Expected: 6 trees and 6 Auspice datasets each; no ignored-task entries for the phylogenetics steps in `pipeline_errors.log`; the terminal `PHYLOGENETICS(H5):` line reporting 6 segment trees with the expected tip counts.

---

### Task 12: Reassortment read-out

**Files:**
- Create: `scripts/segment_placement.py`
- Create: `tests/scripts/test_segment_placement.py`
- Create: `docs/superpowers/results/2026-11-05-reassortment-25CR6244.md`

**Interfaces:**
- Produces: `placement(tree_path: str, focus: str, pathotype_by_tip: dict[str, str]) -> dict` returning, for one tree, the focus tip's sister clade, that clade's pathotype composition, and the bootstrap support of the smallest ancestral node containing the focus tip and at least one `HPAI` tip (`None` when there is none).

- [ ] **Step 1: Write the failing test**

```python
def test_focus_inside_an_hpai_clade_is_reported_with_its_support(self):
    # ((focus,HPAI_a)95,(LPAI_b,LPAI_c)80);
    out = placement(self.tree, "focus", {"HPAI_a": "HPAI", "LPAI_b": "LPAI", "LPAI_c": "LPAI"})
    self.assertEqual(out["sister_pathotypes"], {"HPAI": 1})
    self.assertEqual(out["hpai_node_support"], 95.0)

def test_focus_among_lpai_only_reports_no_hpai_node(self):
    out = placement(self.lpai_tree, "focus", {"LPAI_b": "LPAI", "LPAI_c": "LPAI"})
    self.assertIsNone(out["hpai_node_support"])
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: FAIL — module `segment_placement` not found

- [ ] **Step 3: Implement `placement` with `Bio.Phylo`**

- [ ] **Step 4: Run the test to verify it passes**

Run: `python3 -m unittest discover -s tests/scripts -v`
Expected: PASS

- [ ] **Step 5: Run it over the six combined trees and write the read-out**

One row per segment: sister clade, its pathotype and genotype composition, and the bootstrap support of the smallest clade joining 25CR6244 to any HPAI tip. State for each segment whether reassortment with HPAIV Gs/GD is supported, unsupported or undetermined, and repeat the PA/NA gap.

- [ ] **Step 6: Commit**

```bash
git add scripts/segment_placement.py tests/scripts/test_segment_placement.py docs/superpowers/results/2026-11-05-reassortment-25CR6244.md
git commit -m "feat(scripts): report per-segment placement of 25CR6244 against the HPAI context"
```
