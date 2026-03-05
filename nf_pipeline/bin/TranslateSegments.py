#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Optional
from Bio import SeqIO
from Bio.Seq import Seq

# Aquí guardem el "mapa" de cada segment: on comença i acaba cada proteïna.
ANNOTATIONS = {
    'PB2': {'PB2': [(1, 2280)]},
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
    'NP': {'NP': [(1, 1497)]},
    'NA': {'NA': [(1, 1410)]},
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
    'NS-1': 'NS1',
    'NS-2': 'NEP',
}

# Longituds normals de cada proteïna per saber si la traducció ha anat bé.
EXPECTED_LENGTHS = {
    'PB2': 759, 'PB1': 757, 'PB1-F2': 90, 'PA': 716, 'PA-X': 252,
    'HA': 551, 'HA1': 329, 'HA2': 223, 'NP': 498, 'NA': 469,
    'M1': 252, 'M2': 97, 'NS1': 230, 'NEP': 121,
}

def assemble_nt(nt_seq: str, ranges: list[tuple[int, int]]) -> str:
    """
    Retalla i uneix els trossos de la seqüència segons el mapa de coordenades.
    Si una proteïna està partida en dos trossos, aquí els enganxem.
    """
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
    """
    Passa de lletres de genoma (nucleòtids) a lletres de proteïna (aminoàcids).
    Es para quan troba un senyal de "STOP" i ignora les lletres sobrants del final.
    """
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
    """
    Com que les seqüències a vegades fallen, l'script:
    1. Prova de traduir començant des de la posició 1, 2 i 3. Tria el millor frame
    2. Busca si hi ha un inici oficial (ATG) una mica més endavant.
    3. Si sabem que ha de començar per una lletra concreta (com la D a la HA), la prioritza.
    4. Al final, es queda amb la versió que té la mida més semblant a la real.
    """
    nt_seq = nt_seq.upper().replace('U', 'T')
    candidates = []

    for frame in (0, 1, 2):
        aa = translate_nt(nt_seq[frame:])
        if aa:
            candidates.append(aa)

    search_limit = min(180, max(0, len(nt_seq) - 2))
    for idx in range(search_limit):
        if nt_seq[idx:idx + 3] == 'ATG':
            aa = translate_nt(nt_seq[idx:])
            if aa:
                candidates.append(aa)

    if not candidates:
        return ''

    candidates = list(dict.fromkeys(candidates))

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
    """
    1. Busca els fitxers a la carpeta 'segments'.
    2. Mira de quin segment es tracta.
    3. Tradueix cada proteïna que toca segons el mapa inicial.
    4. Guarda els resultats en fitxers nous dins d'una carpeta 'proteins'.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequences-dir', '--input', dest='sequences_dir', required=True)
    parser.add_argument('--output-dir', '--output', dest='output_dir', required=True)
    args = parser.parse_args()

    sequences_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir)

    for fasta in sequences_dir.rglob('*.fasta'):
        if 'proteins' in fasta.parts or fasta.parent.name != 'segments':
            continue

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

            for protein_name, ranges in ANNOTATIONS[segment].items():
                cds_nt = assemble_nt(nt, ranges)
                out_name = OUTPUT_NAME.get(protein_name, protein_name)
                aa_seq = select_best_translation(cds_nt, EXPECTED_LENGTHS.get(out_name))
                if aa_seq:
                    proteins[out_name] = aa_seq

            # Cas especial: la proteïna HA madura comença per la lletra 'D'.
            if segment == 'HA':
                ha_nt = assemble_nt(nt, [(49, 1704)])
                ha_aa = select_best_translation(ha_nt, EXPECTED_LENGTHS.get('HA'), required_start='D')
                if ha_aa:
                    proteins['HA'] = ha_aa

            for out_name, aa_seq in proteins.items():
                with open(out_dir / f'{out_name}.fasta', 'w') as handle:
                    handle.write(f'>{rec_id}|{out_name}\n{aa_seq}\n')

if __name__ == '__main__':
    main()