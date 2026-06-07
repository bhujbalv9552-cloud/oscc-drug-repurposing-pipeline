import pandas as pd
import os

def process_ilincs_data(input_filepath, output_filepath):
    """
    Reads the iLINCS export CSV, filters for negative concordance (reversers),
    and formats the output so the downstream visualization scripts can read it.
    """
    print("Loading iLINCS connectivity results...")
    
    try:
        # Read the CSV downloaded from iLINCS
        df = pd.read_csv(input_filepath)
    except FileNotFoundError:
        print(f"Error: Could not find {input_filepath}. Ensure the file is saved in the correct location.")
        return

    # Check if this is actually an iLINCS file by verifying key columns
    expected_columns = ['Perturbagen', 'Concordance', 'Perturbagen targets', 'pValue']
    for col in expected_columns:
        if col not in df.columns:
            print(f"Error: Missing expected column '{col}'. Please ensure you downloaded the correct CSV from iLINCS.")
            return

    # Filter for therapeutic candidates (Concordance < 0 indicates reversal)
    reversers = df[df['Concordance'] < 0].copy()

    # Sort the dataframe so the most negative scores (best candidates) are at the top
    reversers = reversers.sort_values(by='Concordance', ascending=True)

    # Rename the iLINCS columns to match what the original CLUE.io parser outputted,
    # ensuring 05_visualize.py does not break.
    formatted_df = reversers.rename(columns={
        'Perturbagen': 'Name',
        'Concordance': 'Score',
        'Perturbagen targets': 'Target',
        'pValue': 'p_value',
        'Cell Line': 'Cell_Line'
    })

    # Select only the most relevant columns for the final output
    cols_to_keep = ['Name', 'Score', 'Target', 'p_value', 'Cell_Line']
    final_cols = [c for c in cols_to_keep if c in formatted_df.columns]
    final_df = formatted_df[final_cols]

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    # Save the processed results
    final_df.to_csv(output_filepath, index=False)
    
    print(f"Success! Processed {len(final_df)} drug candidates.")
    if not final_df.empty:
        top_drug = final_df.iloc[0]['Name']
        top_score = final_df.iloc[0]['Score']
        print(f"Top reversal candidate: {top_drug} (Score: {top_score})")

if __name__ == "__main__":
    # The pipeline expects to find the raw downloaded iLINCS file saved here
    INPUT_CSV = "data/clue_output/clue_results.csv" 
    
    # The pipeline expects the final parsed list to be saved here for the visualizer
    OUTPUT_CSV = "results/top_drug_candidates.csv" 
    
    process_ilincs_data(INPUT_CSV, OUTPUT_CSV)
