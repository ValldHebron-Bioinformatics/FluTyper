#!/usr/bin/env python3
import sys


def rewrite_header(line: str) -> str:
    line = line.rstrip("\n")
    if not line.startswith(">"):
        return line
    parts = line[1:].split("|")

    # Format antic: ...|SEGMENT|...|species_id|clade
    if len(parts) >= 5:
        protein = parts[1]
        species_id = parts[3].replace("_", "")
        clade = parts[4]
    # Format habitual: EPI_ISL_xxx|SEGMENT|strain|clade
    elif len(parts) >= 4:
        species_id = parts[0].replace("_", "")
        protein = parts[1]
        clade = parts[3]
    else:
        return line

    return f">{species_id}_{protein}_{clade}"

def main():
    if len(sys.argv) != 3:
        print("Usage: rewrite_headers.py <input.fasta> <output.fasta>", file=sys.stderr)
        sys.exit(1)

    inp, out = sys.argv[1], sys.argv[2]
    with open(inp, encoding="utf-8") as fin, open(out, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(rewrite_header(line) + "\n")
            else:
                fout.write(line)

if __name__ == "__main__":
    main()