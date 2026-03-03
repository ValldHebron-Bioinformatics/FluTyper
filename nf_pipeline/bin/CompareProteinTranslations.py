#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path
from statistics import mean
from typing import Optional, Tuple

from Bio import SeqIO


ALIASES = {
    'M1': 'MP',
    'M2': 'MP2',
    'NS1': 'NS',
    'NEP': 'NS2',
    'NS2': 'NS2',
}


def normalize_token(value: str) -> str:
    return re.sub(r'[^A-Z0-9]', '', value.upper())


def parse_reference(reference_fasta: Path, selected_epis: set[str]) -> dict[str, dict[str, str]]:
    reference: dict[str, dict[str, str]] = {}
    for rec in SeqIO.parse(str(reference_fasta), 'fasta'):
        header = rec.description.strip()
        if '|' not in header:
            continue
        epi, protein = header.split('|', 1)
        epi = epi.strip().upper()
        if selected_epis and epi not in selected_epis:
            continue
        protein = protein.strip().upper()
        seq = str(rec.seq).upper().replace('*', '')
        reference.setdefault(epi, {})[protein] = seq
    return reference


def load_observed_proteins(sequences_root: Path, epi_id: str) -> Tuple[Optional[str], dict[str, str]]:
    norm_epi = normalize_token(epi_id)
    best_dir = None
    for sample_dir in sequences_root.iterdir():
        if not sample_dir.is_dir():
            continue
        if normalize_token(sample_dir.name).endswith(norm_epi) or norm_epi.endswith(normalize_token(sample_dir.name)):
            best_dir = sample_dir
            break
    if best_dir is None:
        return None, {}

    proteins_dir = best_dir / 'proteins'
    if not proteins_dir.exists():
        return best_dir.name, {}

    observed = {}
    for fp in proteins_dir.glob('*.fasta'):
        recs = list(SeqIO.parse(str(fp), 'fasta'))
        if recs:
            observed[fp.stem.upper()] = str(recs[0].seq).upper().replace('*', '')
    return best_dir.name, observed


def resolve_observed_key(ref_key: str, observed: dict[str, str]) -> Optional[str]:
    if ref_key in observed:
        return ref_key
    alt = ALIASES.get(ref_key)
    if alt and alt in observed:
        return alt
    norm_ref = normalize_token(ref_key)
    for key in observed:
        if normalize_token(key) == norm_ref:
            return key
    return None


def compare_sequences(ref_seq: str, obs_seq: str) -> tuple[int, int, float]:
    min_len = min(len(ref_seq), len(obs_seq))
    mismatches = sum(1 for a, b in zip(ref_seq[:min_len], obs_seq[:min_len]) if a != b)
    extra = len(obs_seq) - len(ref_seq)
    identity = ((min_len - mismatches) / len(ref_seq)) if ref_seq else 0.0
    return mismatches, extra, identity


def build_summary_rows(rows: list[list]) -> list[list]:
    grouped: dict[str, list[list]] = {}
    for row in rows:
        epi_id = str(row[0])
        grouped.setdefault(epi_id, []).append(row)

    summary_rows = []
    for epi_id, epi_rows in sorted(grouped.items()):
        sample_dir = str(epi_rows[0][1]) if epi_rows else '-'
        total = len(epi_rows)
        ok_count = sum(1 for r in epi_rows if r[6] == 'OK')
        diff_count = sum(1 for r in epi_rows if r[6] == 'DIFF')
        missing_count = sum(1 for r in epi_rows if r[6] == 'MISSING')

        valid_identities = []
        for r in epi_rows:
            try:
                valid_identities.append(float(r[9]))
            except (ValueError, TypeError):
                continue

        mean_identity = f"{mean(valid_identities):.4f}" if valid_identities else '-'
        min_identity = f"{min(valid_identities):.4f}" if valid_identities else '-'

        diff_proteins = [str(r[2]) for r in epi_rows if r[6] == 'DIFF']
        missing_proteins = [str(r[2]) for r in epi_rows if r[6] == 'MISSING']
        status = 'OK' if diff_count == 0 and missing_count == 0 else ('MISSING' if missing_count > 0 else 'DIFF')

        summary_rows.append([
            epi_id,
            sample_dir,
            status,
            total,
            ok_count,
            diff_count,
            missing_count,
            mean_identity,
            min_identity,
            ';'.join(diff_proteins) if diff_proteins else '-',
            ';'.join(missing_proteins) if missing_proteins else '-',
        ])

    return summary_rows


def main() -> None:
    parser = argparse.ArgumentParser(description='Compare translated proteins vs reference protein FASTA for EPI samples.')
    parser.add_argument('--reference-fasta', required=True, help='Path to reference protein FASTA (headers like EPI_ISL_xxx|PROTEIN).')
    parser.add_argument('--sequences-root', default='results/converted.fasta/sequences', help='Path to pipeline sequences root containing sample/proteins folders.')
    parser.add_argument('--epi', action='append', default=[], help='Optional EPI_ISL identifier. Can be passed multiple times.')
    parser.add_argument('--output-csv', default=None, help='Optional output CSV with detailed per-protein rows.')
    parser.add_argument('--output-summary-csv', default='comparison_summary.csv', help='Output CSV with one row per EPI (easy to evaluate).')
    args = parser.parse_args()

    reference_fasta = Path(args.reference_fasta).expanduser().resolve()
    sequences_root = Path(args.sequences_root).expanduser().resolve()
    selected_epis = {x.strip().upper() for x in args.epi if x.strip()}

    reference = parse_reference(reference_fasta, selected_epis)
    rows = []

    for epi_id in sorted(reference.keys()):
        sample_dir_name, observed = load_observed_proteins(sequences_root, epi_id)
        for protein, ref_seq in sorted(reference[epi_id].items()):
            obs_key = resolve_observed_key(protein, observed)
            if not obs_key:
                rows.append([epi_id, sample_dir_name or '-', protein, '-', len(ref_seq), '-', 'MISSING', '-', '-', '-'])
                continue

            obs_seq = observed[obs_key]
            mismatches, extra, identity = compare_sequences(ref_seq, obs_seq)
            status = 'OK' if mismatches == 0 and extra == 0 and len(obs_seq) == len(ref_seq) else 'DIFF'
            rows.append([
                epi_id,
                sample_dir_name or '-',
                protein,
                obs_key,
                len(ref_seq),
                len(obs_seq),
                status,
                mismatches,
                extra,
                f'{identity:.4f}',
            ])

    header = ['epi_id', 'sample_dir', 'reference_protein', 'observed_protein', 'ref_len', 'obs_len', 'status', 'mismatches', 'extra_len', 'identity']
    summary_header = ['epi_id', 'sample_dir', 'status', 'total_proteins', 'ok', 'diff', 'missing', 'mean_identity', 'min_identity', 'diff_proteins', 'missing_proteins']

    print('\t'.join(header))
    for row in rows:
        print('\t'.join(map(str, row)))

    summary_rows = build_summary_rows(rows)

    summary_path = Path(args.output_summary_csv).expanduser().resolve()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(summary_header)
        writer.writerows(summary_rows)

    if args.output_csv:
        output_csv = Path(args.output_csv).expanduser().resolve()
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    print(f"\nSummary CSV written to: {summary_path}")


if __name__ == '__main__':
    main()