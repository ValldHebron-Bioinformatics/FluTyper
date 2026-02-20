import json
import sys

def compare_clades(json_file):
    """
    Llegeix un JSON de Nextclade, filtra les seqüències per qualitat i 
    calcula la precisió comparant el clade real amb el predit.
    Si la seqüència no té una qualitat 'good' o el clade és 'unassigned', 
    es descarta de l'anàlisi.
    """
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)

        total_sequences = 0
        correct_predictions = 0

        print("seqID\treal_clade\tpredicted_clade\tmatch")

        for seq in data.get('results', []):
            seq_id = seq.get('seqName', 'N/A')
            predicted_clade = seq.get('clade', 'N/A')
            
            qc_status = seq.get('qc', {}).get('overallStatus', 'N/A')

            if qc_status != 'good':
                continue
            if predicted_clade == 'unassigned':
                continue    
            total_sequences += 1

            real_clade = seq_id.split('|')[-1] if '|' in seq_id else 'N/A'
        
            match = real_clade == predicted_clade
            if match:
                correct_predictions += 1

            print(f"{seq_id}\t{real_clade}\t{predicted_clade}\t{match}")

        accuracy = (correct_predictions / total_sequences) * 100 if total_sequences > 0 else 0
        print(f"\nAccuracy: {accuracy:.2f}% ({correct_predictions}/{total_sequences} correct)")

    except FileNotFoundError:
        print(f"Error: File '{json_file}' not found.")
    except json.JSONDecodeError:
        print(f"Error: File '{json_file}' is not a valid JSON file.")
    except KeyError as e:
        print(f"Error: Missing expected key in JSON data: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python compare_clades.py <nextclade_results.json>")
        sys.exit(1)

    json_file = sys.argv[1]
    compare_clades(json_file)