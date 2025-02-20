import pandas as pd
import subprocess
import argparse
import os

def main(motif_list_sorted, variant_file, output_file, motif_file_dir):
    with open(motif_list_sorted, 'r') as f:
        motif_files = [f"'{line.strip()}'" for line in f if line.strip()]
    
    #results_df = pd.DataFrame()
    result_columns = []
    temp_files = []
    
    # Loop through each motif file
    for motif_file in motif_files:
        motif_name = motif_file.replace('_sorted.bed', '').replace("'", "")
        temp_output = f"temp_distance_{motif_name}.bed"
        temp_files.append(temp_output)
        temp_output_with_quotes = f"'{temp_output}'"
        motif_file_path = os.path.join(motif_file_dir, motif_file)
        command = f"bedtools closest -a {variant_file} -b {motif_file_path} -d -t first > {temp_output_with_quotes}"
        subprocess.run(command, shell=True, check=True)
        temp_df = pd.read_csv(temp_output, sep='\t', header=None)
        distance_column = temp_df.iloc[:, -1]
        #results_df[f"distance_{motif_name}"] = distance_column
        result_columns.append(pd.DataFrame({f"Distance_{motif_name}": distance_column}))

    variants_df = pd.read_csv(variant_file, sep='\t', header=None)
    result_columns.insert(0, pd.DataFrame({"variant_id": variants_df[3]}))
    results_df = pd.concat(result_columns, axis=1)
    #results_df.insert(0, "variant_id", variants_df[3])
    results_df.to_csv(output_file, index=False, sep='\t')
    print(f"Results saved to {output_file}")
    # Delete all the temp files
    for temp_file in temp_files:
        try:
            os.remove(temp_file)
            #print(f"Deleted temporary file: {temp_file}")
        except FileNotFoundError:
            print(f"Temporary file not found (already deleted): {temp_file}")
        except Exception as e:
            print(f"Error deleting temporary file {temp_file}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate distances from variants to motifs and aggregate results.")
    parser.add_argument(
        "--motif_list_sorted", 
        type=str, 
        required=True, 
        help="Path to the file containing sorted motif BED file paths."
    )
    parser.add_argument(
        "--variant_file", 
        type=str, 
        required=True, 
        help="Path to the variants BED file."
    )
    parser.add_argument(
        "--output_file", 
        type=str, 
        required=True, 
        help="Path to the output CSV file."
    )
    parser.add_argument(
        "--motif_file_dir", 
        type=str, 
        required=True, 
        help="Path to the motif file folder path."
    )
    args = parser.parse_args()
    main(args.motif_list_sorted, args.variant_file, args.output_file, args.motif_file_dir)
