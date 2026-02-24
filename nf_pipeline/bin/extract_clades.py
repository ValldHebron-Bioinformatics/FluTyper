import pandas as pd
import sys
import re

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

            # Només analitzem seqüències amb bona qualitat i amb clade assignat
            if qc_status != 'good' or predicted_clade == 'unassigned':
                continue
            
            total_sequences += 1

            # Extraiem l'últim segment de l'ID usant | o _ com a delimitadors
            # La funció re.split permet separar per múltiples caràcters alhora
            segments = re.split(r'[|_]', seq_id)
            real_clade = segments[-1] if segments else 'N/A'
        
            match = (real_clade == predicted_clade)
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