import pandas as pd
import sys

def extract_clade_data(input_csv, output_csv):
    """
    Llegeix el CSV de Nextclade i exporta un nou CSV amb les columnes:
    seqID, predicted_clade i qc.overallStatus.
    """
    try:
        # Carreguem el fitxer original (separat per ;)
        df = pd.read_csv(input_csv, sep=';')

        # Seleccionem i reanomenem les columnes per ajustar-les a la teva petició
        # Nota: 'clade' passa a ser 'predicted_clade' en la sortida
        columns_to_extract = {
            'seqName': 'seqID',
            'clade': 'predicted_clade',
            'qc.overallStatus': 'qc.overallStatus',
            'qc.overallScore': 'qc.overallScore'
        }

        # Verifiquem que les columnes existeixin al fitxer original
        available_cols = [col for col in columns_to_extract.keys() if col in df.columns]
        
        # Filtrem el dataframe
        filtered_df = df[available_cols].rename(columns=columns_to_extract)

        # Exportem a un nou CSV (separat per , per defecte)
        filtered_df.to_csv(output_csv, index=False)
        

    except FileNotFoundError:
        print(f"Error: No s'ha trobat el fitxer '{input_csv}'.")
    except Exception as e:
        print(f"S'ha produït un error: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Ús: python extract_clades.py <fitxer_entrada.csv> [fitxer_sortida.csv]")
        sys.exit(1)

    input_file = sys.argv[1]
    # Si no es dóna nom de sortida, usem un per defecte
    output_file = sys.argv[2] if len(sys.argv) > 2 else "filtered_results.csv"

    extract_clade_data(input_file, output_file)