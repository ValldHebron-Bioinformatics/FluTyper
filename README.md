# FluTyper 🧬🐔🐷

FluTyper is a bioinformatics toolkit designed for genotyping zoonotic influenza viruses (avian and swine) and characterizing relevant mutations.
It aims to support genomic surveillance, research, and integration within a **One Health** framework.

---

## ✨ Features

- 🧬 Influenza genotype assignment
- 🔎 Mutation detection and annotation
- 🧪 Support for avian and swine influenza strains
- 📊 Outputs suitable for surveillance and downstream analysis
- 🌍 Designed with One Health integration in mind


## 🚀 Installation

```bash
git clone https://github.com/ValldHebron-Bioinformatics/FluTyper.git
cd FluTyper
```
## USAGE

For the script to parse your data correctly, MultiFASTA headers must follow a specific naming convention using either an underscore (_) or a pipe (|) as a separator.

Header Format
>{SequenceID}_{Protein}_{OptionalInformation}

SequenceID: A unique identifier for the sample (e.g., Sample01, Chicken02).

Extra Info: Any additional metadata (date, location, etc.).

Examples
Underscore: >Sample01_HA_2024_Spain

Pipe: >Sample01|NA|Hebei_SJ27

## 🔄 Nextflow channel flow (sequences_dir)

## 🔄 Nextflow channel flow (main pipeline)

1. **SampleInput_ch**: Created from the input FASTA file(s) using `channel.fromPath(params.inputFasta, checkIfExists: true).splitFasta(...)`. Emits tuples of (sample_id, fasta_file).
2. **OrganizeBySample(SampleInput_ch)**: Organizes input sequences by sample, emitting (sample_id, sample_dir).
3. **SubtypeDetection(OrganizeBySample.out)**: Receives (sample_id, sample_dir), emits (sample_id, tsv_file) with subtyping results.
4. **GenotypingInfo_ch**: Parses subtyping results to extract (sample_id, h_tag, n_tag, pathotype).
5. **GetDatasets(SubtypeMerged_ch)**: Uses merged subtyping results to determine which reference datasets to download.
6. **GenotypingNextcladeInput_ch**: Joins sample HA files, subtyping info, and datasets, filtering by H-type.
7. **GenotypingNextclade(GenotypingNextcladeInput_ch)**: Runs Nextclade genotyping, emits per-sample CSV results.
8. **NextcladeTuple_ch**: Maps Nextclade output files to (sample_id, csv_file).
9. **GenotypingResultsInput_ch**: Joins genotyping info with Nextclade results, emits (sample_id, h_tag, n_tag, pathotype, csv_file).
10. **GenotypingResults(GenotypingResultsInput_ch, GetDatasets.out.collect())**: Prepares the final report.
11. **GetCDS(CDSInput_ch)**: Prepares inputs for sequence extraction.
