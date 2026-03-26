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
- Mutation detection and annotation (optional, configurable).
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
  --outDir <output_directory>
	[--extraMarkers <folder_with_marker_csvs>]
```
- Default input: `docs/fastas/prova.fasta`
- Default protocol: `AVIAN` (SWINE is under development)

#### Extra Markers


You can provide additional mutation marker data using the `--extraMarkers` flag. This flag should point to a single CSV file containing all extra markers.

**CSV format:**

- The file must have exactly five columns, with headers:
	- `POSITION`, `AA`, `PROTEIN`, `EFFECT`, `REFERENCE`
- Example:
	```csv
	POSITION,AA,PROTEIN,EFFECT,REFERENCE
	155,K,HA1,Antigenic,Reference1
	24,R,M1,CANCER,PMID: 12345678
	190,E,NA,Resistance,Reference2
	```

**Valid protein names:**
HA1, HA2, NA, NP, M1, M2, NS1, NS2, PA, PB1, PB2

Each row should specify the protein in the `PROTEIN` column. The pipeline will automatically detect and use all valid markers from this file.

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

The workflow consists of several key stages:

1. **OrganizeBySample**  
	Organizes input sequences by sample, extracts segments, and creates per-sample directories.
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
	Compares sample proteins to references, annotates mutations, and checks for known markers.
8. **MutationsCompiler**  
	Compiles all mutation data into a single Excel report with one sheet per protein.
9. **CompileErrors**  
	 Aggregates and formats error logs for each sample into a final report.

---

## 📂 Output

- `final_genotyping_results.csv`: Summary of genotyping and QC for all samples.
- `final_mutations_report.xlsx`: All detected mutations, organized by protein.
- `samples/<sample_id>/`: Per-sample folders with intermediate and final sequence files.
- `pipeline_errors.log`: Aggregated error log for the entire run.

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
