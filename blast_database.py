import pandas as pd
from Bio import Entrez, SeqIO
import time
import re

Entrez.email = "marcperez02@gmail.com"

# --- 1. Load and Clean Data ---
df = pd.read_csv('CLADECLASSIFICATION.csv')
df = df.dropna(subset=['link'])

def extract_accession(url):
    match = re.search(r'nuccore/([A-Z0-9_.]+)', str(url).strip())
    if match:
        return match.group(1)
    return None

df['accession'] = df['link'].apply(extract_accession)
df = df.dropna(subset=['accession'])

all_accessions = df['accession'].tolist()
print(f"Found {len(all_accessions)} sequences.")

# --- 2. Batch Download to Memory ---
batch_size = 50
collected_sequences = [] # We store data here instead of writing immediately

print("Starting batch download...")

for i in range(0, len(all_accessions), batch_size):
    batch = all_accessions[i:i + batch_size]
    print(f"Downloading batch {i} to {i + len(batch)}...")
    
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=",".join(batch),
            rettype="fasta",
            retmode="text"
        )
        # Read data into a string variable
        batch_data = handle.read()
        handle.close()
        
        # Store in our list
        collected_sequences.append(batch_data)
        
        time.sleep(1) 
        
    except Exception as e:
        print(f"Error downloading batch starting at {i}: {e}")

# --- 3. Write File Once (Safe for Network Mounts) ---
print(f"Download complete. Writing {len(collected_sequences)} batches to disk...")

try:
    with open("blast_database.fasta", "w") as output_file:
        for seq_batch in collected_sequences:
            output_file.write(seq_batch)
    print("Mission accomplished! blast_database.fasta is ready.")
    
except OSError as e:
    print(f"File system error: {e}")
    print("Try saving to a local folder (like /tmp/) and moving the file manually.")

# --- 4. Add Clade Information to FASTA File ---
def add_clade_to_fasta(fasta_file, csv_file, output_file):
    # Load clade information from CSV
    clade_df = pd.read_csv(csv_file)

    # Normalize the 'accession' column by removing version numbers
    clade_df['accession'] = clade_df['accession'].str.split('.').str[0]  # Remove version numbers

    # Create a mapping of accession to clade
    clade_mapping = dict(zip(clade_df['accession'], clade_df['CLADE']))

    # Read and modify FASTA file
    updated_sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Normalize accession number by removing version numbers
        accession = record.id.split("|")[0].split('.')[0]  # Remove version number
        clade = clade_mapping[accession]  # Get clade (all accessions are guaranteed to have clades)
        record.description += f" | Clade: {clade}"  # Append clade to description
        updated_sequences.append(record)

    # Write updated sequences to a new FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(updated_sequences, output_handle, "fasta")

    print(f"Updated FASTA file written to {output_file}")

# Call the function to update the FASTA file
add_clade_to_fasta("blast_database.fasta", "CLADECLASSIFICATION.csv", "blast_database_with_clades.fasta")