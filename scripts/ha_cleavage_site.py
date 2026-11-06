#!/usr/bin/env python3
"""Read the HA cleavage site motif from a sample's HA protein and call the pathotype from it.

This is a reading of the sequence, not a nearest-reference label. FluTyper's GetCDS and
TranslateToProtein already map each sample's HA onto the protocol's reference CDS and split it
into HA1 and HA2, so the cleavage site is exactly the junction between the two: the motif is
the C-terminal residues of HA1, and HA2 begins with the conserved GLF of the fusion peptide.

The call follows the WOAH/EFSA criterion for H5 and H7: multiple basic residues at the
cleavage site make a virus highly pathogenic. In practice, for an H5:

  HPAI   >= 3 basic residues (K or R) in P1-P5, or a furin consensus R-X-K/R-R at P4-P1
         e.g. ...SPLREKRRKR / GLF   (clade 2.3.4.4b)
  LPAI   a monobasic site, i.e. a basic P1 without that polybasic stretch
         e.g. ...LRNVPQRETR / GLF

Anything that cannot be read is 'indeterminate' and never LPAI: no HA1, no conserved GLF
anchor at the start of HA2, an unresolved residue (X) or a gap inside the motif, an HA1 too
short to hold one, or a P1 that is not basic (which means the junction is not where the
alignment put it).

Usage: ha_cleavage_site.py <results_dir> <out.csv>
       reads <results_dir>/samples/<id>/proteins/<id>_HA1_PROT.fasta (and _HA2_PROT.fasta)
       and writes ID,cleavage_motif,pathotype_cleavage,basic_p1_p5,cleavage_reason
"""
import csv
import glob
import os
import sys

BASIC = set("KR")
ANCHOR = "GLF"          # start of HA2, the fusion peptide
MOTIF_FLANK = 9         # HA1 residues reported before the junction
UNREADABLE = set("X-*?")
OUT_COLUMNS = ["ID", "cleavage_motif", "pathotype_cleavage", "basic_p1_p5", "cleavage_reason"]


def _clean_tail(seq):
    return (seq or "").strip().upper().rstrip("*-")


def read_cleavage(ha1_protein, ha2_protein, flank=MOTIF_FLANK, anchor=ANCHOR):
    """Motif and pathotype call at the HA1/HA2 junction. Never guesses: see module docstring."""
    result = {"motif": "", "call": "indeterminate", "basic_p1_p5": None, "reason": ""}

    ha1 = _clean_tail(ha1_protein)
    ha2 = (ha2_protein or "").strip().upper().lstrip("-")
    if not ha1:
        result["reason"] = "no HA1 protein"
        return result
    if not ha2.startswith(anchor):
        result["reason"] = f"HA2 does not start with the conserved {anchor} anchor"
        return result

    p1_p5 = ha1[-5:]
    if len(ha1) < 5:
        result["reason"] = f"HA1 is only {len(ha1)} residues: too short to read a cleavage site"
        return result

    result["motif"] = f"{ha1[-flank:]}/{ha2[:len(anchor)]}"
    bad = sorted(set(p1_p5) & UNREADABLE)
    if bad:
        result["reason"] = f"unreadable residue(s) {''.join(bad)} in the cleavage motif"
        return result

    basic = sum(1 for aa in p1_p5 if aa in BASIC)
    result["basic_p1_p5"] = basic
    if p1_p5[-1] not in BASIC:
        result["reason"] = f"P1 residue {p1_p5[-1]} is not basic: the HA1/HA2 junction does not look aligned"
        return result

    p4, p2, p1 = p1_p5[-4], p1_p5[-2], p1_p5[-1]
    furin = p4 in BASIC and p2 in BASIC and p1 in BASIC
    if basic >= 3 or furin:
        result["call"] = "HPAI"
        result["reason"] = (f"{basic} basic residues in P1-P5"
                            + (" and a furin consensus R-X-K/R-R at P4-P1" if furin else ""))
    else:
        result["call"] = "LPAI"
        result["reason"] = f"monobasic cleavage site ({basic} basic residues in P1-P5)"
    return result


def _read_protein(path):
    if not path or not os.path.exists(path):
        return ""
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def _one(sample_dir):
    sid = os.path.basename(sample_dir.rstrip("/"))
    ha1 = glob.glob(os.path.join(sample_dir, "proteins", "*_HA1_PROT.fasta"))
    ha2 = glob.glob(os.path.join(sample_dir, "proteins", "*_HA2_PROT.fasta"))
    out = read_cleavage(_read_protein(ha1[0] if ha1 else ""),
                        _read_protein(ha2[0] if ha2 else ""))
    return {"ID": sid, "cleavage_motif": out["motif"], "pathotype_cleavage": out["call"],
            "basic_p1_p5": "" if out["basic_p1_p5"] is None else out["basic_p1_p5"],
            "cleavage_reason": out["reason"]}


def main(results_dir, out_csv):
    sample_dirs = sorted(glob.glob(os.path.join(results_dir, "samples", "*")))
    sample_dirs = [d for d in sample_dirs if os.path.isdir(d)]
    if not sample_dirs:
        sys.exit(f"{results_dir}: no samples/ subdirectories found.")
    rows = [_one(d) for d in sample_dirs]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    counts = {}
    for r in rows:
        counts[r["pathotype_cleavage"]] = counts.get(r["pathotype_cleavage"], 0) + 1
    print(f"{len(rows)} samples -> {out_csv}")
    for call in ("HPAI", "LPAI", "indeterminate"):
        if call in counts:
            print(f"  {call}: {counts[call]}")
    return rows


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
