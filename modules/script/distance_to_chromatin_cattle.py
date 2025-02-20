import os
import argparse
import subprocess
import pandas as pd

def calculate_distances(input_bed, chromatin_dir, output_file):
    TISSUES = ["Adipose", "Cerebellum", "Cortex", "Hypothalamus", "Liver", "Lung", "Muscle", "Spleen"]
    MARKERS = ["ATAC", "CTCF", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3"]

    base_df = pd.read_csv(input_bed, sep='\t', header=None, names=["chrom", "start", "end", "variant_id"])

    for tissue in TISSUES:
        for marker in MARKERS:    
            # Find BED files for the current marker in the tissue folder
            tissue_path = os.path.join(chromatin_dir, tissue)
            if not os.path.exists(tissue_path):
                print(f"  Warning: Tissue folder {tissue_path} does not exist. Skipping.")
                continue
            
            for filename in os.listdir(tissue_path):
                if filename.startswith(marker) and filename.endswith("Peaks.bed"):
                    # Extract sample name from the filename (M08 or M22)
                    parts = filename.split('_')
                    sample = parts[2]
                    column_name = f"distance_to_{marker}_{tissue}_{sample}"
                    bed_file_path = os.path.join(tissue_path, filename)
                    try:
                        result = subprocess.run(
                            ["bedtools", "closest", "-a", input_bed, "-b", bed_file_path, "-D", "b", "-t", "first"],
                            capture_output=True, text=True, check=True
                        )
                    except subprocess.CalledProcessError as e:
                        print(f"Error running bedtools for file {bed_file_path}: {e}")
                        continue
                    
                    distances = [line.split('\t')[-1].strip() for line in result.stdout.strip().split('\n')]
                    base_df[column_name] = distances
    base_df.to_csv(output_file, sep='\t', index=False)
    print(base_df)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate distances from BED file to chromatin data.")
    parser.add_argument("--input_bed", required=True, help="Path to the input BED file.")
    parser.add_argument("--chromatin_dir", required=True, help="Path to the chromatin data directory.")
    parser.add_argument("--output_file", required=True, help="Path to the final output file.")
    args = parser.parse_args()
    calculate_distances(args.input_bed, args.chromatin_dir, args.output_file)
