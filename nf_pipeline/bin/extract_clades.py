import json
import sys

def compare_clades(json_file):
    try:
        # Load the JSON file
        with open(json_file, 'r') as file:
            data = json.load(file)

        total_sequences = 0
        correct_predictions = 0

        print("seqID\treal_clade\tpredicted_clade\tmatch")

        # Compare real clade and predicted clade
        for seq in data['results']:
            seq_id = seq.get('seqName', 'N/A')  # Default to 'N/A' if seqName is missing
            predicted_clade = seq.get('clade', 'N/A')  # Default to 'N/A' if clade is missing
            qc_status = seq.get('status', 'N/A')  # Default to 'N/A' if qc.overallStatus is missing

            # Extract the real clade from the seqID (after the last '/')
            real_clade = seq_id.split('|')[-1] if '|' in seq_id else 'N/A'

            # Check if the predicted clade matches the real clade
            match = real_clade == predicted_clade
            if match:
                correct_predictions += 1
            total_sequences += 1

            # Print the comparison
            print(f"{seq_id}\t{real_clade}\t{predicted_clade}\t{match}")

        # Calculate and print accuracy
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