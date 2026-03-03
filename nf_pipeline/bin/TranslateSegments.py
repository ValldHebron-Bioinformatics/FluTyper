#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from typing import Optional
from Bio import SeqIO
from Bio.Seq import Seq

# Codons de parada estàndard.
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# Regles d'ORFs alternatius per segment.
ALT_ORF_TARGETS = {
    'PB1': {
        'PB1-F2': {
            'target_lengths': [52, 57, 87, 90],
            'expected_start_nt': 95,
            'start_window_nt': 70,
            'max_len_aa': 130,
        },
    },
}

# Regles d'ORFs per frameshift (PA-X).
FRAMESHIFT_ALT_RULES = {
    'PA': {
        'orf_name': 'PA-X',
        'target_len': 252,
        'expected_break_nt': 573,
        'break_window': 18,
        'preferred_motif': 'TCCTTTCGT',
        'motif_break_offset': 6,
        'min_len_aa': 200,
    },
}

# Regles d'ORFs alternatius per splicing.
SPLICE_ALT_RULES = {
    'MP': {
        'orf_name': 'MP2',
        'target_len': 97,
        'prefix_aa': 'MSLL',
        'donor_range': (25, 95),
        'acceptor_range': (650, 900),
        'expected_donor': 51,
        'expected_acceptor': 740,
        'min_intron': 30,
    },
    'NS': {
        'orf_name': 'NS2',
        'target_len': 121,
        'prefix_aa': 'MDSN',
        'donor_range': (25, 95),
        'acceptor_range': (430, 700),
        'expected_donor': 55,
        'expected_acceptor': 529,
        'min_intron': 30,
    },
}

# Límit de candidats de splice per mantenir rendiment.
MAX_SPLICE_DONOR_CANDIDATES = 12
MAX_SPLICE_ACCEPTOR_CANDIDATES = 16


def find_orf_candidates(nt_seq: str) -> list:
    # Troba ORFs iniciades en ATG als 3 frames i en retorna candidates.
    nt_seq = nt_seq.upper().replace('U', 'T')
    candidates = []
    for frame in (0, 1, 2):
        for pos in range(frame, len(nt_seq) - 2, 3):
            if nt_seq[pos:pos + 3] != 'ATG':
                continue
            end = pos
            while end + 3 <= len(nt_seq):
                codon = nt_seq[end:end + 3]
                if end > pos and codon in STOP_CODONS:
                    break
                end += 3
            coding_nt = nt_seq[pos:end]
            if len(coding_nt) < 3:
                continue
            aa = str(Seq(coding_nt).translate(to_stop=False)).replace('*', '')
            if not aa:
                continue
            candidates.append({
                'start_nt': pos + 1,
                'frame': frame + 1,
                'aa': aa,
                'length_aa': len(aa),
            })
    return candidates


def translate_from_first_atg_to_stop(nt_seq: str) -> str:
    # Tradueix des del primer ATG fins al primer codó stop.
    nt_seq = nt_seq.upper().replace('U', 'T')
    first_atg = nt_seq.find('ATG')
    if first_atg == -1:
        return ''
    coding = nt_seq[first_atg:]
    aa = []
    for i in range(0, len(coding) - 2, 3):
        codon = coding[i:i + 3]
        residue = str(Seq(codon).translate())
        if residue == '*':
            break
        aa.append(residue)
    return ''.join(aa)


def translate_from_start_to_stop(nt_seq: str) -> str:
    # Tradueix des de l'inici de seqüència (frame actual) fins a stop.
    nt_seq = nt_seq.upper().replace('U', 'T')
    nt_seq = nt_seq[:len(nt_seq) - (len(nt_seq) % 3)]
    aa = []
    for i in range(0, len(nt_seq), 3):
        residue = str(Seq(nt_seq[i:i + 3]).translate())
        if residue == '*':
            break
        aa.append(residue)
    return ''.join(aa)


def translate_best_frame(nt_seq: str) -> Seq:
    # Selecciona la millor ORF principal (més llarga, inici més primerenc).
    candidates = find_orf_candidates(nt_seq)
    if not candidates:
        return Seq("")
    primary = min(candidates, key=lambda c: (-c['length_aa'], c['start_nt']))
    return Seq(primary['aa'])


def write_position_csv(aa_seq: str, csv_path: Path) -> None:
    # Escriu CSV de posició aminoacídica (1-based).
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['position', 'aa'])
        for idx, aa in enumerate(aa_seq, start=1):
            writer.writerow([idx, aa])


def mature_from_first_d(aa_seq: str) -> str:
    # Retalla la seqüència des del primer D (numeració madura HA).
    first_d = aa_seq.find('D')
    return aa_seq[first_d:] if first_d != -1 else aa_seq


def write_orf_metadata_csv(metadata_rows: list, csv_path: Path) -> None:
    # Guarda metadades d'ORFs detectades per traçabilitat.
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['orf_name', 'start_nt', 'frame', 'length_aa', 'confidence'])
        for row in metadata_rows:
            writer.writerow([row['orf_name'], row['start_nt'], row['frame'], row['length_aa'], row['confidence']])


def find_donor_sites(nt_seq: str) -> list:
    # Donor canònic: límit d'exó abans del motiu intrònic GT.
    return [i for i in range(len(nt_seq) - 1) if nt_seq[i:i + 2] == 'GT']


def find_acceptor_sites(nt_seq: str) -> list:
    # Acceptor canònic: motiu AG just abans de l'inici del següent exó.
    return [i for i in range(2, len(nt_seq) + 1) if nt_seq[i - 2:i] == 'AG']


def select_ranked_sites(sites: list, start: int, end: int, expected: Optional[int], max_sites: int) -> list:
    # Prioritza posicions properes a l'esperada i limita el nombre final.
    in_range = [x for x in sites if start <= x <= end]
    if expected is not None and start <= expected <= end:
        in_range.append(expected)
    ranked = sorted(set(in_range), key=lambda x: (abs(x - expected), x) if expected is not None else (x,))
    return ranked[:max_sites]


def detect_spliced_alt_orf(nt_seq: str, segment_name: str) -> Optional[dict]:
    # Cerca ORF alternativa basada en splice (MP2/NS2).
    rule = SPLICE_ALT_RULES.get(segment_name.upper())
    if not rule:
        return None

    nt_seq = nt_seq.upper().replace('U', 'T')
    donor_start, donor_end = rule['donor_range']
    acceptor_start, acceptor_end = rule['acceptor_range']
    min_intron = rule['min_intron']
    target_len = rule['target_len']
    prefix = rule.get('prefix_aa', '')
    expected_donor = rule.get('expected_donor')
    expected_acceptor = rule.get('expected_acceptor')

    motif_donors = find_donor_sites(nt_seq)
    motif_acceptors = find_acceptor_sites(nt_seq)
    donor_candidates = select_ranked_sites(
        motif_donors,
        donor_start,
        donor_end,
        expected_donor,
        MAX_SPLICE_DONOR_CANDIDATES,
    )
    acceptor_candidates = select_ranked_sites(
        motif_acceptors,
        acceptor_start,
        acceptor_end,
        expected_acceptor,
        MAX_SPLICE_ACCEPTOR_CANDIDATES,
    )

    best = None
    for donor in donor_candidates:
        for acceptor in acceptor_candidates:
            # Evita introns massa curts.
            if acceptor <= donor + min_intron:
                continue
            spliced_nt = nt_seq[:donor] + nt_seq[acceptor:]
            aa = translate_from_first_atg_to_stop(spliced_nt)
            if len(aa) < 10:
                continue

            prefix_mismatches = 0
            if prefix:
                prefix_mismatches = sum(1 for x, y in zip(aa[:len(prefix)], prefix) if x != y)
            donor_has_motif = nt_seq[donor:donor + 2] == 'GT'
            acceptor_has_motif = nt_seq[max(0, acceptor - 2):acceptor] == 'AG'
            motif_penalty = (0 if donor_has_motif else 1) + (0 if acceptor_has_motif else 1)
            donor_dist = abs(donor - expected_donor) if expected_donor is not None else 0
            acceptor_dist = abs(acceptor - expected_acceptor) if expected_acceptor is not None else 0
            # Ordena per longitud objectiu, qualitat de motiu i proximitat.
            score = (abs(len(aa) - target_len), motif_penalty, prefix_mismatches, donor_dist + acceptor_dist)
            if best is None or score < best['score']:
                confidence = 'exact_splice_match' if score == (0, 0, 0, 0) else 'closest_splice_match'
                best = {
                    'score': score,
                    'orf_name': rule['orf_name'],
                    'start_nt': 1,
                    'frame': 1,
                    'length_aa': len(aa),
                    'confidence': confidence,
                    'aa': aa,
                    'donor': donor,
                    'acceptor': acceptor,
                }

    return best


def detect_frameshift_alt_orf(nt_seq: str, segment_name: str, primary_aa: str) -> Optional[dict]:
    # Cerca ORF alternativa per frameshift (PA-X).
    rule = FRAMESHIFT_ALT_RULES.get(segment_name.upper())
    if not rule:
        return None

    nt_seq = nt_seq.upper().replace('U', 'T')
    target_len = rule['target_len']
    center = rule['expected_break_nt']
    window = rule['break_window']
    preferred_motif = rule.get('preferred_motif', '')
    motif_break_offset = rule.get('motif_break_offset', 0)
    min_len_aa = rule.get('min_len_aa', 20)

    # Finestra base al voltant del breakpoint esperat.
    start = max(3, center - window)
    end = min(len(nt_seq) - 3, center + window)

    candidate_breakpoints = set(range(start, end + 1))
    motif_breakpoints = set()
    if preferred_motif:
        # Afegeix candidats guiats pel motiu canònic de frameshift.
        motif_scan_start = max(0, center - 90)
        motif_scan_end = min(len(nt_seq) - len(preferred_motif) + 1, center + 90)
        for i in range(motif_scan_start, motif_scan_end):
            if nt_seq[i:i + len(preferred_motif)] == preferred_motif:
                bp = i + motif_break_offset
                candidate_breakpoints.add(bp)
                motif_breakpoints.add(bp)

    # Candidats finals vàlids.
    candidate_breakpoints = sorted(b for b in candidate_breakpoints if 3 <= b <= len(nt_seq) - 3)

    best = None
    for breakpoint in candidate_breakpoints:
        left_nt = nt_seq[:breakpoint]
        right_nt = nt_seq[breakpoint + 1:]

        aa_left = translate_from_start_to_stop(left_nt)
        aa_right = translate_from_start_to_stop(right_nt)
        if not aa_left or not aa_right:
            continue

        aa = aa_left + aa_right
        if len(aa) < min_len_aa:
            continue
        shared_prefix = 0
        for a, b in zip(aa, primary_aa):
            if a != b:
                break
            shared_prefix += 1

        motif_penalty = 0 if breakpoint in motif_breakpoints else 1

        # Ordena per motiu, longitud objectiu, prefix compartit i proximitat.
        score = (
            motif_penalty,
            abs(len(aa) - target_len),
            -shared_prefix,
            abs(breakpoint - center),
        )
        if best is None or score < best['score']:
            confidence = 'exact_frameshift_match' if score[0] == 0 else 'closest_frameshift_match'
            best = {
                'score': score,
                'orf_name': rule['orf_name'],
                'start_nt': 1,
                'frame': 1,
                'length_aa': len(aa),
                'confidence': confidence,
                'aa': aa,
                'breakpoint': breakpoint,
            }

    return best


def select_alternative_orfs(candidates: list, primary_orf: dict, segment_name: str) -> list:
    # Selecciona ORFs alternatives no-splice/no-frameshift (p.ex. PB1-F2).
    alternatives = []
    targets = ALT_ORF_TARGETS.get(segment_name.upper(), {})
    for orf_name, target_rule in targets.items():
        target_lengths = target_rule['target_lengths'] if isinstance(target_rule, dict) else [target_rule]
        expected_start = target_rule.get('expected_start_nt') if isinstance(target_rule, dict) else None
        start_window = target_rule.get('start_window_nt') if isinstance(target_rule, dict) else None
        max_len = target_rule.get('max_len_aa') if isinstance(target_rule, dict) else None
        eligible = [
            c for c in candidates
            if c['start_nt'] != primary_orf['start_nt'] and c['frame'] != primary_orf['frame'] and c['length_aa'] >= 20
        ]
        base_eligible = eligible[:]
        # Filtre per inici esperat.
        if expected_start is not None and start_window is not None:
            eligible = [c for c in eligible if abs(c['start_nt'] - expected_start) <= start_window]
        # Filtre per longitud màxima.
        if max_len is not None:
            eligible = [c for c in eligible if c['length_aa'] <= max_len]
        # Fallback si el filtre d'inici deixa buit.
        if not eligible:
            eligible = [c for c in base_eligible if max_len is None or c['length_aa'] <= max_len]
        if not eligible:
            continue
        best = min(
            eligible,
            key=lambda c: (
                min(abs(c['length_aa'] - t) for t in target_lengths),
                abs(c['start_nt'] - expected_start) if expected_start is not None else 0,
                c['start_nt'],
            )
        )
        confidence = 'exact_length_match' if best['length_aa'] in target_lengths else 'closest_match'
        alternatives.append({
            'orf_name': orf_name,
            'start_nt': best['start_nt'],
            'frame': best['frame'],
            'length_aa': best['length_aa'],
            'confidence': confidence,
            'aa': best['aa'],
        })

    return alternatives


def convert_fasta(in_fasta: Path, out_fasta: Path, out_csv: Path, use_mature_numbering: bool, segment_name: str) -> None:
    # Converteix una FASTA de segment a proteïna principal + ORFs alternatives.
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    metadata_rows = []

    for record in SeqIO.parse(str(in_fasta), 'fasta'):
        # 1) ORF principal.
        candidates = find_orf_candidates(str(record.seq))
        if not candidates:
            continue

        primary = min(candidates, key=lambda c: (-c['length_aa'], c['start_nt']))
        primary_aa = primary['aa']
        primary_record = record[:]
        primary_record.seq = Seq(primary_aa)
        SeqIO.write([primary_record], str(out_fasta), 'fasta')

        aa_seq = primary_aa
        csv_seq = mature_from_first_d(aa_seq) if use_mature_numbering else aa_seq
        write_position_csv(csv_seq, out_csv)

        metadata_rows.append({
            'orf_name': segment_name,
            'start_nt': primary['start_nt'],
            'frame': primary['frame'],
            'length_aa': primary['length_aa'],
            'confidence': 'primary_orf',
        })

        # 2) ORF per frameshift (PA-X).
        frameshift_alt = detect_frameshift_alt_orf(str(record.seq), segment_name, primary_aa)
        if frameshift_alt:
            alt_record = record[:]
            alt_record.id = f"{record.id}|{frameshift_alt['orf_name']}"
            alt_record.description = ''
            alt_record.seq = Seq(frameshift_alt['aa'])
            alt_fasta = out_fasta.with_name(f"{frameshift_alt['orf_name']}.fasta")
            alt_csv = out_csv.with_name(f"{frameshift_alt['orf_name']}.csv")
            SeqIO.write([alt_record], str(alt_fasta), 'fasta')
            write_position_csv(frameshift_alt['aa'], alt_csv)
            metadata_rows.append({
                'orf_name': frameshift_alt['orf_name'],
                'start_nt': frameshift_alt['start_nt'],
                'frame': frameshift_alt['frame'],
                'length_aa': frameshift_alt['length_aa'],
                'confidence': frameshift_alt['confidence'],
            })

        # 3) ORF per splice (MP2/NS2).
        splice_alt = detect_spliced_alt_orf(str(record.seq), segment_name)
        if splice_alt:
            alt_record = record[:]
            alt_record.id = f"{record.id}|{splice_alt['orf_name']}"
            alt_record.description = ''
            alt_record.seq = Seq(splice_alt['aa'])
            alt_fasta = out_fasta.with_name(f"{splice_alt['orf_name']}.fasta")
            alt_csv = out_csv.with_name(f"{splice_alt['orf_name']}.csv")
            SeqIO.write([alt_record], str(alt_fasta), 'fasta')
            write_position_csv(splice_alt['aa'], alt_csv)
            metadata_rows.append({
                'orf_name': splice_alt['orf_name'],
                'start_nt': splice_alt['start_nt'],
                'frame': splice_alt['frame'],
                'length_aa': splice_alt['length_aa'],
                'confidence': splice_alt['confidence'],
            })

        # 4) Altres ORFs alternatives (PB1-F2).
        for alt in select_alternative_orfs(candidates, primary, segment_name):
            alt_record = record[:]
            alt_record.id = f"{record.id}|{alt['orf_name']}"
            alt_record.description = ''
            alt_record.seq = Seq(alt['aa'])
            alt_fasta = out_fasta.with_name(f"{alt['orf_name']}.fasta")
            alt_csv = out_csv.with_name(f"{alt['orf_name']}.csv")
            SeqIO.write([alt_record], str(alt_fasta), 'fasta')
            write_position_csv(alt['aa'], alt_csv)
            metadata_rows.append({
                'orf_name': alt['orf_name'],
                'start_nt': alt['start_nt'],
                'frame': alt['frame'],
                'length_aa': alt['length_aa'],
                'confidence': alt['confidence'],
            })

    if metadata_rows:
        write_orf_metadata_csv(metadata_rows, out_csv.with_name(f"{segment_name}_orfs.csv"))


def main() -> None:
    # Punt d'entrada CLI.
    parser = argparse.ArgumentParser(description='Translate segment FASTA files to amino-acid FASTA files.')
    parser.add_argument('--sequences-dir', required=True, help='Path to sequences directory generated by OrganizeBySpecies')
    parser.add_argument('--output-dir', required=True, help='Output base directory where sequences/<species>/proteins/*.fasta will be written')
    parser.add_argument('--protocol', default='', help='Protocol name (e.g. H5, H7). Mature numbering from first D is applied only for H5/H7')
    args = parser.parse_args()

    sequences_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir)
    protocol = args.protocol.strip().upper()
    # Numeració madura només per protocols H5/H7.
    use_mature_numbering = protocol in {'H5', 'H7'}

    for species_dir in sequences_dir.iterdir():
        if not species_dir.is_dir():
            continue
        segment_dirs = [species_dir / 'segment', species_dir / 'segments']
        for segments_dir in segment_dirs:
            if not segments_dir.exists():
                continue
            for seg_fasta in segments_dir.glob('*.fasta'):
                segment_name = seg_fasta.stem.strip().upper()
                out_fasta = output_dir / species_dir.name / 'proteins' / f"{segment_name}.fasta"
                out_csv = output_dir / species_dir.name / 'proteins' / f"{segment_name}.csv"
                # Aplicar numeració madura només al segment HA.
                use_mature_for_segment = use_mature_numbering and segment_name == 'HA'
                convert_fasta(seg_fasta, out_fasta, out_csv, use_mature_for_segment, segment_name)


if __name__ == '__main__':
    main()