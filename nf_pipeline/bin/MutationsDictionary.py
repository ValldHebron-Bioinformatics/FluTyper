#!/usr/bin/env python3
import argparse, csv, re
from pathlib import Path
import pandas as pd

COL = {"H3":"reference_site(H3_numbering)","H1":"reference_H1_site(H1_numbering)","H5":"mature_H5_site(no_signal_peptide)","H7":"H7_NUMBERING","H9":"H9_NUMBERING"}
RX = re.compile(r"^\s*(-?\d+[A-Za-z]?)\s*([A-Za-z])\s*$")

def parse_args() -> argparse.Namespace:
    """Llegeix els arguments de comanda (subtipus destí i sortida opcional)."""
    parser = argparse.ArgumentParser(description="Translate HA markers (H5 mature numbering) to a target subtype")
    parser.add_argument("--subtype", required=True, choices=sorted(COL))
    parser.add_argument("--output", default=None, help="Optional output CSV path")
    return parser.parse_args()


def build_mapping(dictionary_path: Path, subtype: str) -> dict:
    """Construeix el mapa de posicions H5 madures cap al subtipus objectiu."""
    h5_to_target = {}
    with dictionary_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            h5 = (row.get(COL["H5"]) or "").strip()
            target = (row.get(COL[subtype]) or "").strip()
            if h5 and target:
                h5_to_target[h5] = target
    return h5_to_target


def translate_markers(markers_path: Path, subtype: str, h5_to_target: dict) -> pd.DataFrame:
    """Llegeix MARKERS.xlsx, filtra HA i retorna només les columnes H5 i subtipus."""
    df = pd.read_excel(markers_path, sheet_name="Markers", header=None)
    translated = []
    for col_idx in range(df.shape[1]):
        if str(df.iat[1, col_idx]).strip().upper() != "HA":
            continue
        for row_idx in range(2, df.shape[0]):
            match = RX.match(str(df.iat[row_idx, col_idx]).strip())
            if not match:
                continue
            h5_pos, aa = match.group(1), match.group(2).upper()
            target_pos = h5_to_target.get(h5_pos, "")
            translated.append({"H5": f"{h5_pos}{aa}", subtype: f"{target_pos}{aa}" if target_pos else "-"})
    return pd.DataFrame(translated, columns=["H5", subtype])


def main() -> None:
    """Orquestra la càrrega de fitxers, la traducció i l'escriptura de resultats."""
    args = parse_args()
    markers = Path("docs/mutations/MARKERS.xlsx")
    dictionary = Path("docs/mutations/AA_dictionary_proposal/AA_Sites.csv")
    h5_to_target = build_mapping(dictionary, args.subtype)
    result = translate_markers(markers, args.subtype, h5_to_target)

    if args.output:
        output_path = Path(args.output).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_path, index=False)
    else:
        print(result.to_csv(index=False), end="")


if __name__ == "__main__":
    main()
