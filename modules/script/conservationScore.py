#!/usr/bin/env python

import pandas as pd
import pyBigWig
import sys
import argparse
import os

def get_conservation_score(variants_path, conservation_file_path, conservation_type):
    """
    Extract conservation scores for each variant from a given bigWig file.
    :param variants_path: Path to the file containing variant positions.
    :param conservation_file_path: Path to the conservation bigWig file.
    :param conservation_type: The type of conservation score being extracted.
    :return: A DataFrame with conservation scores.
    """
    bw = pyBigWig.open(conservation_file_path)
    variants = pd.read_table(variants_path, header=None, names=None)
    conservation = pd.DataFrame(columns=[conservation_type])

    for i, row in variants.iterrows():
        chrom = row[0]
        start = int(float(row[1]))
        end = int(float(row[2]))

        # Get the conservation score for the specified region
        score = bw.values(chrom, start, end)
        score_str = ''.join(str(s) for s in score)
        conservation.loc[i] = [score_str]

    bw.close()
    return conservation

def process_human_conservation(args, initial_dataset):
    """
    Process and annotate human conservation scores.
    :param args: Parsed command-line arguments.
    :param initial_dataset: DataFrame containing variant data.
    :return: DataFrame with annotated conservation scores.
    """
    conservation_scores = []

    if args.phastCons100way_path:
        conservation_scores.append(get_conservation_score(args.data_path, args.phastCons100way_path, "phastCons100way"))
    if args.phyloP100way_path:
        conservation_scores.append(get_conservation_score(args.data_path, args.phyloP100way_path, "phyloP100way"))
    if args.phastCons30way_path:
        conservation_scores.append(get_conservation_score(args.data_path, args.phastCons30way_path, "phastCons30way"))
    if args.phyloP30way_path:
        conservation_scores.append(get_conservation_score(args.data_path, args.phyloP30way_path, "phyloP30way"))

    annotation_result = pd.concat([initial_dataset] + conservation_scores, axis=1)
    return annotation_result

def process_cattle_conservation(args, initial_dataset):
    """
    Process and annotate cattle conservation scores.
    :param args: Parsed command-line arguments.
    :param initial_dataset: DataFrame containing variant data.
    :return: DataFrame with annotated conservation scores.
    """
    cattle_results = []

    for _, row in initial_dataset.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']

        phastCons_file = os.path.join(args.phastCons241way_path, f'phastCons_bosTau9_{chrom}.bw')
        phyloP_file = os.path.join(args.phyloP241way_path, f'phyloP_bosTau9_{chrom}.bw')

        phastCons_value = 0
        phyloP_value = 0

        if os.path.exists(phastCons_file):
            with pyBigWig.open(phastCons_file) as bw_phastCons:
                phastCons_value = bw_phastCons.values(chrom, start, end)[0]

        if os.path.exists(phyloP_file):
            with pyBigWig.open(phyloP_file) as bw_phyloP:
                phyloP_value = bw_phyloP.values(chrom, start, end)[0]

        cattle_results.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'variant_id': row['variant_id'],
            'phastCons241way': phastCons_value,
            'phyloP241way': phyloP_value
        })

    return pd.DataFrame(cattle_results)

def main():
    parser = argparse.ArgumentParser(description="Annotate variants with conservation scores.")
    parser.add_argument("species", help="Species used in the analysis (human or cattle)")
    parser.add_argument("data_path", help="Input variant position file path")
    parser.add_argument("output_path", help="Result output file path")
    parser.add_argument("--phastCons100way_path", help="Human phastCons100way bigwig file")
    parser.add_argument("--phastCons30way_path", help="Human phastCons30way bigwig file")
    parser.add_argument("--phyloP100way_path", help="Human phyloP100way bigwig file")
    parser.add_argument("--phyloP30way_path", help="Human phyloP30way bigwig file")
    parser.add_argument("--phastCons241way_path", help="Cattle phastCons241way bigwig files folder path")
    parser.add_argument("--phyloP241way_path", help="Cattle phyloP241way bigwig files folder path")
    args = parser.parse_args()

    # Load the variant data
    try:
        bed_file = pd.read_table(args.data_path, header=None, names=['chrom', 'start', 'end', 'variant_id'])
    except Exception as e:
        sys.exit(f"Error reading data file: {e}")

    # Annotate variants based on species
    if args.species.lower() == 'human':
        annotation_result = process_human_conservation(args, bed_file)
    elif args.species.lower() == 'cattle':
        annotation_result = process_cattle_conservation(args, bed_file)
    else:
        sys.exit(f"Unsupported species: {args.species}")

    # Save the results to the output file
    try:
        annotation_result.to_csv(args.output_path, header=True, index=False, sep='\t')
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")

if __name__ == "__main__":
    main()