# FluTyper 🧬🐔🐷
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)

FluTyper is a modular, reproducible Nextflow pipeline for genotyping influenza viruses (avian and human) and characterizing relevant mutations. It is designed for genomic surveillance, research, and integration within a **One Health** framework.

---

## ✨ Features

- Automated organization of input samples and extraction of individual segments.
- Subtype detection (H/N typing and pathotype inference) from sequence data.
- Reference dataset selection and download based on detected subtypes.
- Genotyping using Nextclade with per-sample and merged reports.
- Extraction of coding sequences (CDS) and translation to protein sequences.
- Mutation detection and annotation with standardized cross-subtype numbering (optional, configurable).
- Aggregate, per-sample, and time-series HTML mutation reports.
- Comprehensive error reporting and logging.
- Support for avian and human influenza workflows (SWINE protocol is still in development).
- Modular, reproducible workflow built with Nextflow DSL2.

---

## 🚀 Quick Start

### Installation
Clone the repository to your local machine:
```bash
git clone https://github.com/ValldHebron-Bioinformatics/FluTyper.git
cd FluTyper
```

### Basic Execution
Run the pipeline with your sequence data. If no parameters are provided, it defaults to the `AVIAN` protocol using the provided test dataset.
```bash
nextflow run nf_pipeline/main.nf \
    --inputFasta <path/to/input.fasta> \
    --outDir <path/to/output_directory>
```

### Advanced Execution
For a fully customized run utilizing all available features and the colorblind-friendly reporting mode:
```bash
nextflow run nf_pipeline/main.nf \
  --inputFasta <input.fasta> \
  --outDir <output_directory> \
  --protocol HUMAN \
  --extraMarkers <extra_markers.csv> \
  --metadata <metadata.csv> \
  --threshold 0.25 \
  --colorblind true
  --IndividualReports true
```

---

## 📥 Input Data Preparation

### FASTA Headers
MultiFASTA headers must use either an underscore (`_`) or a pipe (`|`) as a separator to ensure the pipeline correctly parses the sequence identity and segment. 

**Format:** `{SequenceID}_{Segment}_{OptionalInformation}` or `{SequenceID}|{Segment}|{OptionalInformation}`
*   **Underscore Example:** `>Sample01_HA_2024_Spain`
*   **Pipe Example:** `>Sample01|NA|Hebei_SJ27`

### Metadata CSV (Optional)
To generate date-based frequency reports per protein, you must provide a metadata CSV file using the `--metadata` flag. The file requires strict headers:
```csv
ID,DATE
Sample_01,YYYY-MM-DD
```

---

## 🛠️ Command-Line Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `inputFasta` | Required | `docs/fastas/prova.fasta` | Path to the input FASTA file containing your sequences. |
| `outDir` | Required | `RESULTS` | Destination directory where all pipeline results and reports will be saved. |
| `protocol` | Optional | `AVIAN` | Defines the viral protocol. Supported options are `AVIAN` and `HUMAN`. |
| `threshold` | Optional | `0.25` | Minimum mutation frequency required to report a non-marker mutation (0.0 to 1.0). |
| `extraMarkers` | Optional | *None* | Path to a custom CSV file containing specific markers to track across samples. |
| `metadata` | Optional | *None* | Path to a metadata CSV to enable time-series frequency tracking. |
| `colorblind` | Optional | `false` | Set to `true` to apply an Okabe-Ito colorblind-friendly palette to all HTML charts. |
| `IndividualReports` | Optional | `false` | Set to `true` to make the Individual genomic barcode for each sample. |

---

## 🔬 Advanced Configuration & Behaviors

### Threshold Parameter Behavior
The threshold parameter establishes a baseline frequency cutoff that impacts data visualization. This threshold value initializes the slider in the `MutationsReport.html` interactive plot, allowing users to dynamically adjust the view without needing to re-execute the pipeline. The frequency denominator used for this calculation is the number of samples containing that specific protein, rather than the total number of samples in the entire run.

### H Subtype Threshold Summary

| H subtype | HA1 threshold (%) | HA2 threshold (%) | Samples (n) |
| :--- | ---: | ---: | ---: |
| 1 | 68.46 | 76.57 | 17 |
| 2 | 78.83 | 83.15 | 21 |
| 3 | 73.89 | 80.53 | 18 |
| 4 | 76.91 | 81.13 | 24 |
| 5 | 73.59 | 79.86 | 49 |
| 6 | 72.40 | 79.86 | 27 |
| 7 | 70.77 | 77.44 | 50 |
| 8 | 79.97 | 83.45 | 8 |
| 9 | 77.71 | 81.53 | 30 |
| 10 | 76.00 | 80.08 | 34 |
| 11 | 75.82 | 77.86 | 29 |
| 12 | 78.63 | 80.76 | 18 |
| 13 | 72.50 | 78.01 | 9 |
| 14 | 85.48 | 87.14 | 6 |
| 15 | 86.22 | 88.18 | 4 |
| 16 | 89.04 | 92.59 | 3 |

### N Subtype Threshold Summary

| N subtype | NA threshold (%) | Samples (n) |
| :--- | ---: | ---: |
| 1 | 75.87 | 41 |
| 2 | 77.36 | 78 |
| 3 | 70.20 | 34 |
| 4 | 80.59 | 19 |
| 5 | 76.57 | 21 |
| 6 | 76.06 | 31 |
| 7 | 72.32 | 30 |
| 8 | 74.36 | 31 |
| 9 | 79.05 | 62 |

### HUMAN Protocol Notes
The HUMAN protocol utilizes dedicated resources located under `protocols/HUMAN/v1` and introduces marker annotations specifically tailored to human seasonal influenza. It supports genotyping for `H1` (using Nextclade dataset `flu_h1n1pdm_ha`) and `H3` (using Nextclade dataset `flu_h3n2_ha`). Marker files are read directly from the protocol's marker directory rather than querying FluMutDB. When metadata is provided, the human protocol fully supports generating time-evolution frequency reports.

### Integrating Extra Markers

You can introduce your own mutation markers into the pipeline by supplying a structured CSV file via the `--extraMarkers` parameter. The file must strictly contain seven columns: `MARKER_ID`, `POSITION`, `AA`, `PROTEIN`, `EFFECT`, `FOUND_IN`, and `REFERENCE`. It is crucial to note that the `POSITION` value must always be specified using H5N1 numbering. Custom `MARKER_ID` values must be integers starting at 1000. Valid inputs for the `PROTEIN` column include HA1, HA2, M1, M2, NA, NP, NS-1, NS-2, PA, PB1, PB1-F2, and PB2.

You may use "X" in the `AA` column as a wildcard, which forces the pipeline to trigger the marker upon any true amino-acid change at that specific position. This is particularly useful for tracking biologically relevant epitopes regardless of the specific resulting mutation. To define a marker that requires a combination of mutations, simply assign the same `MARKER_ID` to multiple rows. For the HUMAN protocol, ensure the `FOUND_IN` column is accurately populated, as the pipeline will selectively check samples against markers that match their specific subtype context.

```csv
MARKER_ID,POSITION,AA,PROTEIN,EFFECT,FOUND_IN,REFERENCE
1000,631,L,PB2,Increased pandemic risk,H5N1,Capalastegui & Goldhill 2025
1001,141,X,HA1,RBD,H5N1 | H7N9, Luczo & Spackman 2024

```

---

## 🔄 Pipeline Architecture

![FluTyper pipeline walkthrough](docs/TFM/FluTyper.drawio.png)

| Step | Process Name | Description |
| :--- | :--- | :--- |
| **1** | **OrganizeBySample** | Organizes input sequences, detects orientation, extracts segments, and builds directories. |
| **2** | **SubtypeDetection** | Uses minimizer-based subtyping via Nextclade to infer H/N subtypes and pathotypes. |
| **3** | **DB & Dataset Prep** | Fetches the latest FluMutDB, generates protein-specific markers, and downloads references. |
| **4** | **GenotypingNextclade** | Runs Nextclade genotyping for each sample against the assigned reference datasets. For H5N1 clade 2.3.4.4b it leverages Genin2 to assign the genotype. |
| **5** | **GenotypingResults** | Aggregates all genotyping outputs into a unified summary report. |
| **6** | **GetCDS** | Maps and extracts the coding sequences (CDS) using the reference alignments. |
| **7** | **TranslateToProtein** | Translates the aligned CDS nucleotides into amino acid sequences. |
| **8** | **MutationsFinder** | Compares samples to references, annotates mutations, and flags known marker hits. |
| **9** | **MutationsCompiler** | Compiles all mutation data into a comprehensive Excel report. |
| **10** | **CompileErrors** | Aggregates and formats all operational error logs into a final text report. |
| **11-15** | **Graphic Reports** | Generates interactive HTML dashboards for clades, overall mutations, markers, and timelines. |

---

## 🔢 Standardized Cross-Subtype Numbering

Mutation markers are matched using unified reference numbering based on H5 for HA proteins and N1 for NA proteins. During the mutation finding step, the pipeline invokes a specialized dictionary script for any sample whose detected subtype differs from the reference. The script performs a lookup in the corresponding dictionary and populates a dual-coordinate system in the final output. The `POSITION_SUBTYPE` column receives the native subtype-specific residue number, while the `POSITION` column retains the standardized H5 or N1 coordinate. This mechanism ensures results can be interpreted natively while remaining directly comparable against published literature across different influenza subtypes.

| Dictionary | Reference | Subtypes Covered 
| :--- | :--- | :--- | 
| **HA_DICT.csv** | H5 | H1, H2, H3, H4, H6, H7, H8, H9, H10, H11, H12, H13, H14, H15, H16, H17, H18 |
| **NA_DICT.csv** | N1 | N2, N3, N4, N5, N6, N7, N8, N9, N10, N11 |

---

## 📂 Outputs

![FluTyper output folder organization](docs/TFM/Folderorganization.drawio.png)

### Core Data Files
| File Name | Description |
| :--- | :--- |
| `final_genotyping_results.csv` | Summary of subtype inference, clade assignments and genotyping for clade 2.3.4.4b.
| `final_mutations_report.xlsx` | Exhaustive record of all detected mutations, organized by protein sheets. |
| `pipeline_errors.log` | Aggregated error and warning log detailing any operational issues during the run. |
| `samples/<sample_id>/` | Individual directories containing intermediate sequences, alignments, and specific data. |

### Interactive HTML Reports
| Report Path | Description |
| :--- | :--- |
| `graphic_reports/CladeGraphicReport.html` | Interactive visualization of subtype and clade distributions. |
| `graphic_reports/MutationsReport.html` | Aggregate per-protein mutation frequencies with dynamic threshold controls. |
| `graphic_reports/MutationsTable.html` | Interactive table detailing marker effects, subtypes, and references. |
| `graphic_reports/FrequencyEvolution/**/*.html` | Time-series plots showing marker frequency over time (requires metadata). |
| `samples/<id>/<id>_MutationsReport.html` | Per-sample mutation barcode plots for rapid visual inspection. |

### Excel Data Schema (`final_mutations_report.xlsx`)
| Column Header | Description |
| :--- | :--- |
| `SAMPLE_ID` | The unique identifier of the sample. |
| `SUBTYPE` | The complete detected subtype (e.g., H5N1(HPAI)). |
| `PROTEIN` | The specific viral protein where the mutation occurs. |
| `REF_SUBTYPE` | The reference subtype used for the baseline alignment. |
| `POSITION` | The standardized reference coordinate used for reporting. |
| `POSITION_REF` | The internal coordinate used for marker matching. |
| `REFERENCE_AA` | The baseline amino acid present in the reference sequence. |
| `QUERY_AA` | The detected amino acid present in the sample sequence. |
| `AA_MUTATION` | Formatted mutation label (e.g., N30D). |
| `MUTATION_TYPE` | Categorized as Substitution, Insertion, Deletion, or Marker. |
| `MARKER` | Boolean flag (Yes/No) indicating a match with a known marker database entry. |
| `MARKER_ID` | The unique identifier(s) of the matched marker(s). |
| `IS_COMBINATION` | Boolean flag (Yes/No) indicating if the marker requires a multi-mutation pattern. |
| `EFFECT` | Biological or phenotypic effect annotation. |
| `FOUND_IN` | The specific viral subtype context where the marker was originally reported. |
| `REFERENCE` | Source literature or database reference validating the marker. |

---

## 🧪 Testing & Continuous Integration

FluTyper is strictly verified using `nf-test`. The repository utilizes GitHub Actions to execute automated Continuous Integration (CI) on `ubuntu-latest` environments for every push and pull request. The CI handles the installation of all necessary bioinformatics dependencies and executes both individual module unit tests and comprehensive end-to-end integration tests against reference datasets.

*   **Run all tests:** `nf-test test tests/main.nf.test`
*   **Run module tests:** `nf-test test tests/modules/*.nf.test`
*   **Run specific module:** `nf-test test tests/modules/<module_name>.nf.test`

---

## 🧩 Dependencies & Acknowledgments

| Software / Library | Usage |
| :--- | :--- |
| **[Nextflow](https://docs.seqera.io/nextflow/?__hstc=247481240.afc94a4be2e71d336bddb8a957545fad.1771512569721.1777885792784.1778573663745.21&__hssc=247481240.1.1778573663745&__hsfp=e02a5757d31090b6dc84c4b8b9f6ddac)** | Core workflow execution and orchestration engine. |
| **[Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/)** | Sequence genotyping, alignment, and clade assignment. |
| **[Genin2](https://izsvenezie-virology.github.io/genin2/)** | Genotype prediction for clade 2.3.4.4b. |
| **[Seqkit](https://bioinf.shenwei.me/seqkit/usage/)** | High-performance sequence parsing and FASTA manipulation. |
| **[MAFFT](https://mafft.cbrc.jp/alignment/software/)** | Multiple sequence alignment for accurate CDS mapping. |
| **Python 3** | Data manipulation and reporting ([`pandas`](https://pandas.pydata.org/docs/user_guide/index.html#user-guide), [`biopython`](https://biopython.org/docs/latest/index.html), [`openpyxl`](https://openpyxl.readthedocs.io/en/stable/), [`sqlite3`](https://docs.python.org/3/library/sqlite3.html), [`plotly`](https://plotly.com/python/)). |
| **[nf-test](https://www.nf-test.com/docs/getting-started/)** | Pipeline testing and validation framework. |

*The minimizer indices used by this pipeline were generated using the methodology and tools developed by the Nextstrain team for the [nextclade_data](https://github.com/nextstrain/nextclade_data.git) repository.*