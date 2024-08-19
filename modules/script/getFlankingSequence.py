#!/usr/bin/env python

import pandas as pd
import subprocess
import sys
import math

def extract_flanking_sequence(chrom, pos, flanking_len, reference_file):
    """
    Extract the flanking sequence around a given position using samtools.
    """
    coord = f"{chrom}:{pos - flanking_len}-{pos + flanking_len}"
    result = subprocess.run(
        ['samtools', 'faidx', reference_file, coord],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    # Split the result to get the sequence, handle potential decoding issues
    sequence_lines = result.stdout.decode('utf-8').splitlines()
    if len(sequence_lines) >= 2:
        return sequence_lines[1].upper()
    else:
        return None

def process_variants(data_path, flanking_len, reference_file, output_path):
    """
    Process the variant file to extract flanking sequences for each variant.
    """
    # Load the input data
    try:
        df = pd.read_table(data_path, header=None, names=['chrom', 'start', 'end', 'variant_id'])
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")

    flanking_sequences = []

    # Iterate through the DataFrame and process each variant
    for index, row in df.iterrows():
        if pd.isna(row['start']):
            flanking_sequences.append(float('nan'))
        else:
            try:
                flanking_sequence = extract_flanking_sequence(row['chrom'], row['end'], flanking_len, reference_file)
                flanking_sequences.append(flanking_sequence)
            except subprocess.CalledProcessError as e:
                print(f"Error processing variant {row['variant_id']} at {row['chrom']}:{row['end']}: {e}", file=sys.stderr)
                flanking_sequences.append(None)
    
    # Add the flanking sequences to the DataFrame
    df['flanking_sequence'] = flanking_sequences

    # Save the annotated data to the output file
    try:
        df.to_csv(output_path, header=True, index=False, sep='\t')
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")

def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: script.py <data_path> <flanking_len> <reference_file> <output_path>")

    data_path = sys.argv[1]
    flanking_len = int(sys.argv[2])
    reference_file = sys.argv[3]
    output_path = sys.argv[4]

    process_variants(data_path, flanking_len, reference_file, output_path)

if __name__ == "__main__":
    main()
