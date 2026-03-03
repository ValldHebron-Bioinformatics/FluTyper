## This is a copy of https://github.com/nextstrain/nextclade_data/blob/master/scripts/lib/minimizer.py
import copy
import os

import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import json
import click
import yaml

# Increment this version, if there are changes in the algorithm.
# Backwards compatibility must be ensured to not break client code: all versions must be computed and reffered to
# from the dataset's server 'index.json' file.
MINIMIZER_ALGO_VERSION = "1"

# Increment this version, if there are changes in the layout of the output file.
MINIMIZER_JSON_SCHEMA_VERSION = "3.0.0"

# What is this?
MAGIC_NUMBER_K = 17


# from lh3
def invertible_hash(x):
    m = (1 << 32) - 1
    x = (~x + (x << 21)) & m
    x = x ^ (x >> 24)
    x = (x + (x << 3) + (x << 8)) & m
    x = x ^ (x >> 14)
    x = (x + (x << 2) + (x << 4)) & m
    x = x ^ (x >> 28)
    x = (x + (x << 31)) & m
    return x


# turn a kmer into an integer
def get_hash(kmer, cutoff):
    x = 0
    j = 0
    for i, nuc in enumerate(kmer):
        if i % 3 == 2:
            continue  # skip every third nucleotide to pick up conserved patterns
        if nuc not in "ACGT":
            return cutoff + 1  # break out of loop, return hash above cutoff
        else:  # A=11=3, C=10=2, G=00=0, T=01=1
            if nuc in "AC":
                x += 1 << j
            if nuc in "AT":
                x += 1 << (j + 1)
        j += 2

    return invertible_hash(x)


def get_ref_search_minimizers(seq: SeqRecord, cutoff: int, k=MAGIC_NUMBER_K):
    seq_str = preprocess_seq(seq)
    minimizers = []
    # we know the rough number of minimizers, so we can pre-allocate the array if needed
    for i in range(len(seq_str) - k):
        kmer = seq_str[i : i + k]
        mhash = get_hash(kmer, cutoff)
        if (
            mhash < cutoff
        ):  # accept only hashes below cutoff --> reduces the size of the index and the number of look-ups
            minimizers.append(mhash)
    return np.unique(minimizers)


def make_ref_search_index(refs, cutoff: int):
    # collect minimizers for each reference sequence first
    minimizers_by_reference = list()
    for name, ref in sorted(refs.items()):
        minimizers = get_ref_search_minimizers(ref, cutoff)
        minimizers_by_reference.append(
            {
                "minimizers": minimizers,
                "meta": {
                    "name": name,
                    "length": len(ref.seq),
                    "nMinimizers": len(minimizers),
                },
            }
        )

    # construct an index where each minimizer maps to the references it contains via a bit set (here boolean np array)
    index = {"minimizers": {}, "references": []}
    for ri, minimizer_set in enumerate(minimizers_by_reference):
        for m in minimizer_set["minimizers"]:
            if m not in index["minimizers"]:
                index["minimizers"][m] = []
            index["minimizers"][m].append(ri)

        # reference will be a list in same order as the bit set
        index["references"].append(minimizer_set["meta"])

    normalization = np.array(
        [x["length"] / x["nMinimizers"] for x in index["references"]]
    )

    return {
        "schemaVersion": MINIMIZER_JSON_SCHEMA_VERSION,
        "version": MINIMIZER_ALGO_VERSION,
        "params": {
            "k": MAGIC_NUMBER_K,
            "cutoff": cutoff,
        },
        **index,
        "normalization": normalization,
    }


def preprocess_seq(seq: SeqRecord) -> str:
    return str(seq.seq).upper().replace("-", "")


def serialize_ref_search_index(index):
    index = copy.deepcopy(index)
    index["minimizers"] = {str(k): v for k, v in index["minimizers"].items()}
    index["normalization"] = index["normalization"].tolist()
    return index


def to_bitstring(arr) -> str:
    return "".join([str(x) for x in arr])


@click.command(help="Parse fasta header, only keep if fits regex filter_fasta_headers")
@click.option("--references-fasta", required=True, type=click.Path(exists=True))
@click.option("--minimizer-json", required=True, type=click.Path())
@click.option("--config-file", required=True, type=click.Path(exists=True))
def main(references_fasta: str, minimizer_json: str, config_file: str) -> None:
    with open(config_file, encoding="utf-8") as file:
        full_config = yaml.safe_load(file)
        cutoff_threshold = full_config.get("threshold", 28)
        if cutoff_threshold < 1 or cutoff_threshold > 32:
            raise ValueError(
                f"Invalid threshold value {cutoff_threshold}. Must be between 1 and 31."
            )
    # this is equivalent to 2^(cutoff_threshold) or 1 followed by (cutoff_threshold) zeros in binary
    cutoff = 1 << cutoff_threshold
    refs = {rec.id: rec for rec in SeqIO.parse(references_fasta, "fasta")}

    index = make_ref_search_index(refs, cutoff)

    serialized = serialize_ref_search_index(index)
    os.makedirs(os.path.dirname(minimizer_json), exist_ok=True)

    with open(minimizer_json, "w") as f:
        json.dump(serialized, f)


if __name__ == "__main__":
    main()