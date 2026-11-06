#!/usr/bin/env python3
"""Where does one sample's segment sit, relative to the HPAI tips of the same tree?

For each segment tree this reports, for the focus sample:

  sister_tips         the tips of its sister clade (its closest relatives in that tree)
  sister_pathotypes   their pathotype composition; a tip with no pathotype is UNKNOWN, never
                      assumed to be one or the other
  sister_genotypes    their genotype composition, when genotypes are supplied
  hpai_node_support   the bootstrap support of the SMALLEST clade that contains the focus tip
                      and at least one HPAI tip, or None when the tree has no HPAI tip at all
                      or carries no support values

A focus segment whose sister clade is HPAI, under a well-supported node, is what reassortment
with an HPAI virus looks like in a single-segment tree. Nothing here decides that: it reports
the placement and the support, and the interpretation is written up separately.

Usage: segment_placement.py <pathotype.csv> <focus_id> <tree.nwk> [<tree.nwk> ...]
       The pathotype CSV needs an ID column (ID/Sample_ID/SampleID) and a pathotype column;
       CLADE and GENOTYPE columns are used for the composition when present.
"""
import os
import sys
from collections import Counter

from Bio import Phylo


def _support(clade):
    """Bootstrap support of an internal node, or None when the tree carries none.

    Newick internal labels land in confidence or in name depending on the file, and a label
    that is not a number is a name, not a support value."""
    for value in (clade.confidence, clade.name):
        if value is None:
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return None


def placement(tree_path, focus, pathotype_by_tip, genotype_by_tip=None):
    genotype_by_tip = genotype_by_tip or {}
    tree = Phylo.read(tree_path, "newick")

    focus_clade = next((t for t in tree.get_terminals() if t.name == focus), None)
    result = {
        "tree": os.path.basename(os.path.dirname(os.path.abspath(tree_path))) or tree_path,
        "focus_present": focus_clade is not None,
        "n_tips": len(tree.get_terminals()),
        "sister_tips": [],
        "sister_pathotypes": {},
        "sister_genotypes": {},
        "sister_support": None,
        "hpai_node_support": None,
        "hpai_node_tips": None,
    }
    if focus_clade is None:
        return result

    # Ancestors of the focus tip, innermost first: its parent, then outwards to the root
    path = tree.get_path(focus_clade)
    ancestors = list(reversed(path[:-1])) + [tree.root]
    parent = ancestors[0]

    sister_clades = [c for c in parent.clades if c is not focus_clade]
    sister_tips = [t.name for c in sister_clades for t in c.get_terminals() if t.name]
    result["sister_tips"] = sister_tips
    result["sister_pathotypes"] = dict(Counter(
        (pathotype_by_tip.get(name) or "UNKNOWN").upper() for name in sister_tips))
    if genotype_by_tip:
        result["sister_genotypes"] = dict(Counter(
            genotype_by_tip.get(name, "UNKNOWN") for name in sister_tips))
    result["sister_support"] = _support(parent)

    # Smallest clade holding the focus tip and at least one HPAI tip
    for node in ancestors:
        tips = [t.name for t in node.get_terminals() if t.name and t.name != focus]
        if any((pathotype_by_tip.get(name) or "").upper() == "HPAI" for name in tips):
            result["hpai_node_support"] = _support(node)
            result["hpai_node_tips"] = len(tips) + 1
            break
    return result


def _read_pathotypes(path):
    import pandas as pd
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    df.columns = [str(c).strip() for c in df.columns]
    id_col = next((c for c in df.columns
                   if c.lower().replace("_", "").replace(" ", "") in ("id", "sampleid")), None)
    if id_col is None:
        sys.exit(f"{path}: no sample ID column. Columns: {list(df.columns)}")
    patho = next((c for c in df.columns if c.lower() in ("pathotype", "pathogenicity")), None)
    if patho is None:
        sys.exit(f"{path}: no pathotype/PATHOGENICITY column. Columns: {list(df.columns)}")
    geno = next((c for c in df.columns if c.lower() == "genotype"), None)
    ids = df[id_col].str.strip()
    return (dict(zip(ids, df[patho].str.strip())),
            dict(zip(ids, df[geno].str.strip())) if geno else {})


def main(pathotype_csv, focus, tree_paths):
    pathotypes, genotypes = _read_pathotypes(pathotype_csv)
    print(f"focus: {focus}\n")
    for tree_path in tree_paths:
        r = placement(tree_path, focus, pathotypes, genotypes)
        print(f"== {r['tree']}  ({r['n_tips']} tips)")
        if not r["focus_present"]:
            print(f"   {focus} is not a tip of this tree\n")
            continue
        support = "no support values" if r["sister_support"] is None \
            else f"support {r['sister_support']:g}"
        print(f"   sister clade: {len(r['sister_tips'])} tip(s), {support}")
        print(f"   sister pathotypes: {r['sister_pathotypes'] or 'none'}")
        if r["sister_genotypes"]:
            print(f"   sister genotypes:  {r['sister_genotypes']}")
        if r["hpai_node_support"] is None and r["hpai_node_tips"] is None:
            print("   no HPAI tip anywhere in this tree")
        else:
            node_support = "no support value" if r["hpai_node_support"] is None \
                else f"support {r['hpai_node_support']:g}"
            print(f"   smallest clade with an HPAI tip: {r['hpai_node_tips']} tips, {node_support}")
        print()


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
