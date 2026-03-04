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

Protein: Must be either HA or NA.

Extra Info: Any additional metadata (date, location, etc.).

Examples
Underscore: >Sample01_HA_2024_Spain

Pipe: >Sample01|NA|Hebei_SJ27
