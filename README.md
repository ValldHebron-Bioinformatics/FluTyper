# FluTyper 🧬🐔🐷
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)

FluTyper is a modular, reproducible Nextflow pipeline for genotyping zoonotic influenza viruses (avian and swine) and characterizing relevant mutations. It is designed for genomic surveillance, research, and integration within a **One Health** framework.

---

## ✨ Features

- Automated organization of input samples and extraction of individual segments.
- Subtype detection (H/N typing and pathotype inference) from sequence data.
- Reference dataset selection and download based on detected subtypes.
- Genotyping using Nextclade with per-sample and merged reports.
- Extraction of coding sequences (CDS) and translation to protein sequences.
- Mutation detection and annotation with standardized cross-subtype numbering (optional, configurable).
- Comprehensive error reporting and logging.
- Support for both avian and swine influenza viruses (SWINE protocol in development).
- Modular, reproducible workflow built with Nextflow DSL2.

---

## 🚀 Installation

```bash
git clone https://github.com/ValldHebron-Bioinformatics/FluTyper.git
cd FluTyper
```

---

## 🏃 Usage

### Input Requirements
MultiFASTA headers must use either an underscore (`_`) or a pipe (`|`) as a separator.

**Header Format:**
- `{SequenceID}_{Protein}_{OptionalInformation}`
- `{SequenceID}|{Protein}|{OptionalInformation}`

**Examples:**
- Underscore: `>Sample01_HA_2024_Spain`
- Pipe: `>Sample01|NA|Hebei_SJ27`

### Execution
Run the pipeline with the following command:
```bash
nextflow run nf_pipeline/main.nf \
  --inputFasta <input.fasta> \
  --protocol <AVIAN|SWINE> \
  --outDir <output_directory> \
	--extraMarkers <extra_markers.csv> \
  --threshold <0-1>
```
- Default input: `docs/fastas/prova.fasta`
- Default protocol: `AVIAN` (SWINE is under development)

#### Extra Markers

You can provide additional mutation marker data using the `--extraMarkers` flag. This flag should point to a single CSV file containing all extra markers.

**CSV format:**

- The file must have exactly seven columns, with headers:
	- `MARKER_ID`, `POSITION`, `AA`, `PROTEIN`, `EFFECT`, `FOUND_IN`, `REFERENCE`
- `MARKER_ID` must be an integer and should start at `1000` for custom markers.
- Use the same `MARKER_ID` across multiple rows when a marker is defined by a combination of mutations.
- `FOUND_IN` should indicate the subtype/context where that marker effect was reported.
- Example:
	```csv
	MARKER_ID,POSITION,AA,PROTEIN,EFFECT,FOUND_IN,REFERENCE
	1000,155,K,HA1,Antigenic,H5N1,Reference1
	1001,24,R,M1,CANCER,H1N1,PMID:12345678
	1002,190,E,NA,Resistance,H7N9,Reference2
	1010,30,D,M1,Composite marker example,H5N1,PMID:11111111
	1010,31,N,M1,Composite marker example,H5N1,PMID:11111111
	```

**Valid protein names:**
HA1, HA2, M1, M2, NA, NP, NS-1, NS-2, PA, PB1, PB1-F2, PB2

Each row should specify the protein in the `PROTEIN` column. The pipeline will automatically detect and use all valid markers from this file.
Rows that share the same `MARKER_ID` are linked as a single marker definition, allowing the pipeline to evaluate combined mutation patterns together.
If required columns are missing or the file cannot be parsed, the pipeline prints a warning and skips loading extra markers from that file.

### Threshold parameter for mutation relevance

You can configure the threshold for reporting mutations that are frequent within the same protein using the `--threshold` parameter (default: 0.25).

- This parameter determines the minimum fraction of samples containing a specific protein in which a mutation at a given position must appear to be considered relevant and included in the `relevant_mutations.xlsx` report.
- Because the pipeline adjusts for missing data, the threshold is calculated independently for each protein based only on the number of samples where that protein was successfully identified rather than the total number of samples in the run.
- For example, if the default value is 0.25 and a specific protein is only detected in 40 out of 100 total samples, a mutation must appear in more than 10 of those 40 samples to be reported as relevant.
- You can adjust this value when running the pipeline:

```bash
nextflow run nf_pipeline/main.nf --threshold 0.5
```

This would only report mutations present in more than 50% of all run samples as relevant.

If some samples lack a given protein/segment, they are still part of the denominator in the current behavior.

### Testing
The project uses `nf-test` for verification.
- **All tests:** `nf-test test tests/main.nf.test`
- **Module tests:** `nf-test test tests/modules/*.nf.test`
- **Individual module test:** `nf-test test tests/modules/<module_name>.nf.test`

### 🛠️ Continuous Integration
FluTyper uses GitHub Actions for automated testing and quality assurance:
- **Unit Testing:** Individual Nextflow modules are tested using `nf-test` on every pull request and push to the main branches.
- **Integration Testing:** The complete pipeline is verified against test datasets to ensure end-to-end functionality.

---

## 🔄 Pipeline Overview

![FluTyper pipeline walkthrough](docs/TFM/FluTyper.drawio.png)

The workflow consists of several key stages:

1. **OrganizeBySample**  
	Organizes input sequences by sample, automatically detects correct sequence orientation (handling reverse complements), extracts segments, and creates per-sample directories.
2. **SubtypeDetection**  
	Uses Nextclade minimizer-based subtyping to infer H/N subtypes and pathotypes (H5/H7/H9).
3.  **Database & Dataset Management:**
    - **FluMutDB:** Automatically fetches or updates the latest `flumut_db.sqlite` from the [izsvenezie-virology/FluMutDB](https://github.com/izsvenezie-virology/FluMutDB) repository.
    - **MarkersFiles:** Queries the SQLite database to generate protein-specific marker CSVs for mutation annotation.
    - **GetDatasets:** Downloads/selects Nextclade reference datasets based on detected subtypes.
4. **GenotypingNextclade**  
	Runs Nextclade genotyping for each sample using the selected datasets.
5. **GenotypingResults**  
	Merges genotyping results into a summary report (clade, QC, dataset version).
6. **GetCDS & TranslateToProtein**  
	Extracts coding sequences (CDS) using reference alignments and translates them to proteins.
7. **MutationsFinder**  
	Compares sample proteins to references, applies the designated HA and NA numeration schema for standardized coordinates, annotates mutations, and checks for known markers.
8. **MutationsCompiler**  
	Compiles all mutation data into a single Excel report with one sheet per protein.
9. **CompileErrors**  
	Aggregates and formats error logs for each sample into a final report.

---

## 🔢 Standardized Cross-Subtype Numbering

All mutation markers stored in `flumut_db.sqlite` use a unified reference numbering based on **H5 (for HA proteins) and N1 (for NA protein)**. To enable biologically meaningful cross-subtype comparisons, FluTyper includes two positional correspondence dictionaries that translate residue positions from H5/N1 coordinates into the equivalent positions of any other detected subtype.

### HA Dictionary (`HA_DICT.csv`)

This dictionary covers the hemagglutinin protein and maps H5-based residue positions to the equivalent positions in 17 additional HA subtypes (H1–H18, excluding H5 itself as the reference). It contains 585 alignment positions across the following subtypes:

| Reference | Subtypes covered |
|-----------|-----------------|
| H5 | H1, H2, H3, H4, H6, H7, H8, H9, H10, H11, H12, H13, H14, H15, H16, H17, H18 |

### NA Dictionary (`NA_DICT.csv`)

This dictionary covers the neuraminidase protein and maps N1-based residue positions to the equivalent positions in 10 additional NA subtypes (N2–N11). It contains 493 alignment positions across the following subtypes:

| Reference | Subtypes covered |
|-----------|-----------------|
| N1 | N2, N3, N4, N5, N6, N7, N8, N9, N10, N11 |

### How it works

During the **MutationsFinder** step, once mutations are identified against the H5 or N1 reference, the `MutationsDictionary.py` script is invoked for any sample whose detected subtype differs from the reference. The script performs a lookup in the appropriate dictionary and populates the `POSITION_SUBTYPE` column in the output CSV with the subtype-specific residue number, while the `POSITION` column retains the original H5/N1 coordinate. This dual-coordinate system allows results to be interpreted both in the standardized reference framework and in the subtype-native numbering, facilitating direct comparison with published literature for any influenza subtype.

---

## 📂 Output

![FluTyper output folder organization](docs/TFM/Folderorganization.drawio.png)

- `final_genotyping_results.csv`: Summary of genotyping and QC for all samples.
- `final_mutations_report.xlsx`: All detected mutations, organized by protein.
- `relevant_mutations.xlsx`: Filtered mutation report including markers and mutations above the configured relevance threshold.
- `samples/<sample_id>/`: Per-sample folders with intermediate and final sequence files.
- `pipeline_errors.log`: Aggregated error log for the entire run.

### `final_mutations_report.xlsx` columns

The `All_Proteins` sheet and each protein-specific sheet contain the same columns:

- `SAMPLE_ID`: Sample identifier.
- `SUBTYPE`: Detected subtype for the sample (for example `H3N5` or `H5N1(HPAI)`).
- `PROTEIN`: Protein where the mutation/marker was found.
- `REF_SUBTYPE`: Reference subtype used during the alignment for that protein.
- `POSITION`: Residue position in the target numbering used for reporting.
- `POSITION_REF`: Standardized reference position used internally for marker matching (H5 for HA, N1 for NA).
- `REFERENCE_AA`: Amino acid in the reference sequence at that position.
- `QUERY_AA`: Amino acid in the sample sequence at that position.
- `AA_MUTATION`: Mutation label. For substitutions it is formatted like `N30D`; marker-only entries can appear as `<position><AA>`. The position used is the one found in `POSITION`.
- `MUTATION_TYPE`: One of `Substitution`, `Insertion`, `Deletion`, or `Marker`.
- `MARKER`: `Yes` if the event matches a known marker definition, otherwise `No`.
- `MARKER_ID`: Marker identifier(s). Multiple IDs are separated by ` | `.
- `IS_COMBINATION`: `Yes` when the matched marker ID represents a multi-mutation combination, otherwise `No`.
- `EFFECT`: Marker effect annotation(s), if available.
- `FOUND_IN`: Subtype where the marker was reported in the source data.
- `REFERENCE`: Literature/source reference(s) associated with the marker.

---

## ⚙️ Configuration

- Edit `nf_pipeline/nextflow.config` to customize protocols, reference paths, and default parameters.
- Protocol-specific resources are managed in the `protocols/` directory.

---

### Continuous Integration (CI)
The project uses GitHub Actions (`.github/workflows/ci.yml`) to automate:
- **Module Unit Tests:** Using `nf-test` for individual modules.
- **Full Pipeline Integration Tests:** End-to-end validation of the workflow.
CI runs on `ubuntu-latest` and handles the installation of all necessary bioinformatics tools (Nextclade, seqkit, MAFFT) and Python dependencies.

---

## 🧩 Dependencies

- [Nextflow](https://www.nextflow.io/)
- [Nextclade](https://clades.nextstrain.org/)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Python 3](https://www.python.org/) (with [pandas](https://pandas.pydata.org/docs/index.html), [biopython](https://biopython.org/wiki/Documentation), [openpyxl](https://openpyxl.readthedocs.io/en/stable/), [sqlite3](https://docs.python.org/3/library/sqlite3.html))
- [nf-test](https://www.nf-test.com/)