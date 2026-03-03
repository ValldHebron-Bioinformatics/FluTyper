#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Optional
from Bio import SeqIO
from Bio.Seq import Seq

# Coordenades d'anotació per proteïna (1-based, extretes de FluMut)
ANNOTATIONS = {
    'PB2': {
        'PB2': [(1, 2280)],
    },
    'PB1': {
        'PB1': [(1, 2274)],
        'PB1-F2': [(95, 367)],
    },
    'PA': {
        'PA': [(1, 2151)],
        'PA-X': [(1, 570), (572, 760)],
    },
    'HA': {
        'HA1': [(49, 1035)],
        'HA2': [(1036, 1704)],
    },
    'NP': {
        'NP': [(1, 1497)],
    },
    'NA': {
        'NA': [(1, 1410)],
    },
    'MP': {
        'M1': [(1, 759)],
        'M2': [(1, 26), (715, 982)],
    },
    'NS': {
        'NS-1': [(1, 693)],
        'NS-2': [(1, 30), (503, 838)],
    },
}

OUTPUT_NAME = {
    # Normalitza noms interns a noms de sortida esperats
    'NS-1': 'NS1',
    'NS-2': 'NEP',
}

# Longituds de referència aproximades per detectar traduccions massa curtes
EXPECTED_LENGTHS = {
    'PB2': 759,
    'PB1': 757,
    'PB1-F2': 90,
    'PA': 716,
    'PA-X': 252,
    'HA': 551,
    'HA1': 329,
    'HA2': 223,
    'NP': 498,
    'NA': 469,
    'M1': 252,
    'M2': 97,
    'NS1': 230,
    'NEP': 121,
}


def assemble_nt(nt_seq: str, ranges: list[tuple[int, int]]) -> str:
    # Uneix trams nucleotídics segons les coordenades anotades
    seq = str(nt_seq).upper().replace('U', 'T')
    seq_len = len(seq)
    chunks = []
    for start, end in ranges:
        if start > seq_len:
            continue
        safe_start = max(1, start)
        safe_end = min(end, seq_len)
        if safe_start > safe_end:
            continue
        chunks.append(seq[safe_start - 1:safe_end])
    return ''.join(chunks)


def translate_nt(nt_seq: str) -> str:
    # Tradueix fins al primer codó stop; retalla cua no múltiple de 3
    if not nt_seq:
        return ''
    nt_seq = nt_seq.upper().replace('U', 'T')
    usable_len = (len(nt_seq) // 3) * 3
    if usable_len == 0:
        return ''
    return str(Seq(nt_seq[:usable_len]).translate(to_stop=True))


def select_best_translation(
    nt_seq: str,
    expected_len: Optional[int] = None,
    required_start: Optional[str] = None,
) -> str:
    nt_seq = nt_seq.upper().replace('U', 'T')
    candidates = []

    # Prova els tres marcs
    for frame in (0, 1, 2):
        aa = translate_nt(nt_seq[frame:])
        if aa:
            candidates.append(aa)

    # Prova també inicis ATG propers al principi
    search_limit = min(180, max(0, len(nt_seq) - 2))
    for idx in range(search_limit):
        if nt_seq[idx:idx + 3] == 'ATG':
            aa = translate_nt(nt_seq[idx:])
            if aa:
                candidates.append(aa)

    if not candidates:
        return ''

    # Elimina duplicats mantenint ordre
    candidates = list(dict.fromkeys(candidates))

    # Si cal un inici específic (HA amb D), prioritza aquests candidats
    if required_start:
        preferred = [aa for aa in candidates if aa.startswith(required_start)]
        if preferred:
            candidates = preferred
        else:
            shifted = []
            for aa in candidates:
                pos = aa.find(required_start)
                if pos != -1:
                    shifted.append(aa[pos:])
            if shifted:
                candidates = list(dict.fromkeys(shifted))

    if expected_len:
        return min(candidates, key=lambda aa: abs(len(aa) - expected_len))

    return max(candidates, key=len)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequences-dir', '--input', dest='sequences_dir', required=True)
    parser.add_argument('--output-dir', '--output', dest='output_dir', required=True)
    parser.add_argument('--protocol', default='')
    args = parser.parse_args()
    protocol = args.protocol.strip().upper()

    sequences_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir)

    for fasta in sequences_dir.rglob('*.fasta'):
        # Processa només FASTA de segments (evita rellegir proteins/*.fasta)
        if 'proteins' in fasta.parts or fasta.parent.name != 'segments':
            continue

        # El segment s'infereix del nom del fitxer (p. ex. HA.fasta)
        segment = fasta.stem.upper()
        if segment not in ANNOTATIONS:
            continue

        species = fasta.parent.parent.name
        out_dir = output_dir / species / 'proteins'
        out_dir.mkdir(parents=True, exist_ok=True)

        for record in SeqIO.parse(fasta, 'fasta'):
            nt = str(record.seq)

            rec_id = record.id.split('|')[0] if '|' in record.id else record.id

            proteins = {}
            # Traducció general per anotacions
            for protein_name, ranges in ANNOTATIONS[segment].items():
                cds_nt = assemble_nt(nt, ranges)
                out_name = OUTPUT_NAME.get(protein_name, protein_name)
                aa_seq = select_best_translation(cds_nt, EXPECTED_LENGTHS.get(out_name))
                if not aa_seq:
                    continue
                proteins[out_name] = aa_seq

            if segment == 'HA':
                # HA completa des de la regió madura anotada (sense pèptid senyal)
                ha_nt = assemble_nt(nt, [(49, 1704)])
                ha_aa = select_best_translation(
                    ha_nt,
                    EXPECTED_LENGTHS.get('HA'),
                    required_start='D',
                )
                if ha_aa:
                    proteins['HA'] = ha_aa

            # Escriu cada proteïna en un FASTA independent
            for out_name, aa_seq in proteins.items():
                with open(out_dir / f'{out_name}.fasta', 'w') as handle:
                    handle.write(f'>{rec_id}|{out_name}\n{aa_seq}\n')


if __name__ == '__main__':
    main()
