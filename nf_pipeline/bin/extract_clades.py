import pandas as pd
import sys
import re


def infer_real_clade(seq_id: str) -> str:
    pipe_parts = seq_id.split('|')
    if pipe_parts:
        candidate = pipe_parts[-1].strip()
        if candidate.startswith('2.3.4.4'):
            return candidate

    match = re.search(r'2\.3\.4\.4[\w.-]*', seq_id)
    if match:
        return match.group(0)

    return 'unassigned'

def compare_clades(csv_file):
    """
    Llegeix un CSV de Nextclade i calcula la precisió comparant el clade 
    real (últim element després de | o _) amb el clade predit.
    """
    try:
        # Carreguem les dades especificant el delimitador de punt i coma
        df = pd.read_csv(csv_file, sep=';')

        total_sequences = 0
        correct_predictions = 0

        # Definim les capçaleres de la taula per a una visualització clara
        print(f"{'seqID':<60}\t{'real_clade':<15}\t{'predicted_clade':<15}\t{'match'}")

        for index, row in df.iterrows():
            seq_id = str(row.get('seqName', 'N/A'))
            predicted_clade = str(row.get('clade', 'N/A'))
            qc_status = str(row.get('qc.overallStatus', 'N/A'))

            # Només analitzem seqüències amb bona qualitat
            if qc_status != 'good':
                continue

            real_clade = infer_real_clade(seq_id)

            total_sequences += 1

            is_2344_like_match = real_clade.startswith('2.3.4.4') and predicted_clade == '2.3.4.4-like'
            match = (real_clade == predicted_clade) or is_2344_like_match
            if match:
                correct_predictions += 1

            print(f"{seq_id:<60}\t{real_clade:<15}\t{predicted_clade:<15}\t{match}")

        # Calculem el percentatge d'encert final
        accuracy = (correct_predictions / total_sequences) * 100 if total_sequences > 0 else 0
        print(f"\nAccuracy: {accuracy:.2f}% ({correct_predictions}/{total_sequences} correct)")

    except FileNotFoundError:
        print(f"Error: No s'ha pogut trobar el fitxer '{csv_file}'.")
    except Exception as e:
        print(f"S'ha produït un error durant el processament: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Ús: python compare_clades.py <nextclade_results.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    compare_clades(input_file)