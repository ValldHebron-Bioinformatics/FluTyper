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
{SequenceID}_{Protein}_{OptionalInformation}

SequenceID: A unique identifier for the sample (e.g., Sample01, Chicken02).

Extra Info: Any additional metadata (date, location, etc.).

Examples
Underscore: >Sample01_HA_2024_Spain

Pipe: >Sample01|NA|Hebei_SJ27

## 🔄 Nextflow channel flow (sequences_dir)

1. input_ch = channel.of([params.sample, params.dirSample])
2. OrganizeBySpecies(input_ch) creates and emits path("sequences")
3. TranslateToProtein(OrganizeBySpecies.out) passes that emitted path into the process input
4. Inside TranslateToProtein, path(sequences_dir) binds that incoming path to the variable sequences_dir
5. The script uses ${sequences_dir} for --sequences-dir and --output-dir, then emits path("sequences") again
