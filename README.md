
# FluTyper 🧬🐔🐷

FluTyper is a modular, reproducible Nextflow pipeline for genotyping zoonotic influenza viruses (avian and swine) and characterizing relevant mutations. It is designed for genomic surveillance, research, and integration within a **One Health** framework.

---

## ✨ Features

- Automated organization of input samples and extraction of individual segments
- Subtype detection (H/N typing and pathotype inference) from sequence data
- Reference dataset selection and download based on detected subtypes
- Genotyping using Nextclade with per-sample and merged reports
- Extraction of coding sequences (CDS) and translation to protein sequences
- Mutation detection and annotation (optional, configurable)
- Comprehensive error reporting and logging
- Outputs designed for genomic surveillance, research, and downstream analysis
- Support for both avian and swine influenza viruses (SWINE protocol in development)
- Modular, reproducible workflow built with Nextflow DSL2

---

## 🚀 Installation

```bash
git clone https://github.com/ValldHebron-Bioinformatics/FluTyper.git
cd FluTyper
```

---

## 🏃 Usage

**Input requirements:**  
MultiFASTA headers must use either an underscore (_) or a pipe (|) as a separator.

Header Format:
- `{SequenceID}_{Protein}_{OptionalInformation}`
- `{SequenceID}|{Protein}|{OptionalInformation}`

Examples:
- Underscore: `>Sample01_HA_2024_Spain`
- Pipe: `>Sample01|NA|Hebei_SJ27`

**Run the pipeline:**
```bash
nextflow run nf_pipeline/main.nf --inputFasta <your_fasta_file> --protocol AVIAN --outDir <output_folder>
```
- Default input: `docs/fastas/prova.fasta`
- Default protocol: `AVIAN`
- Output: `prova/1/`

---

## 🔄 Pipeline Overview

1. **OrganizeBySample**  
	Organizes input sequences by sample, extracts segments, and creates per-sample directories.

2. **SubtypeDetection**  
	Combines HA and NA segments, runs Nextclade minimizer-based subtyping, and infers H/N subtypes and pathotype (H5/H7/H9).

3. **GetDatasets**  
	Determines which reference datasets to use based on detected subtypes.

4. **GenotypingNextclade**  
	Runs Nextclade genotyping for each sample using the appropriate dataset.

5. **GenotypingResults**  
	Merges genotyping results, extracts clade, QC status, and dataset version, and prepares a final CSV report.

6. **GetCDS**  
	Extracts coding sequences (CDS) for each segment using reference alignments.

7. **TranslateToProtein**  
	Translates CDS FASTA files to protein sequences.

8. **MutationsFinder**  
	Compares sample proteins to references, annotates mutations, and checks for known markers.

9. **MutationsCompiler**  
	Compiles all mutation CSVs into a single Excel report with one sheet per protein.

10. **CompileErrors**  
	 Aggregates and formats error logs for each sample.

---

## 📂 Output

- `final_genotyping_results.csv`: Summary of genotyping and QC for all samples
- `final_mutations_report.xlsx`: All detected mutations, one sheet per protein
- `samples/<sample_id>/`: Per-sample folders with all intermediate and final files
- `pipeline_errors.log`: Aggregated error log for the run

---

## ⚙️ Configuration

- Edit `nf_pipeline/nextflow.config` to customize protocols, reference paths, and default parameters.
- Only the AVIAN protocol is currently supported; SWINE is under development.

---

## 🧩 Dependencies

- [Nextflow](https://www.nextflow.io/)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [Nextclade](https://clades.nextstrain.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Python 3](https://www.python.org/) (with pandas, biopython, openpyxl)

---
