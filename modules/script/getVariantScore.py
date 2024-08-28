#!/usr/bin/env python

import tensorflow as tf
import tensorflow_hub as hub
import joblib
import gzip
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from enformerModels import Enformer
from variantProcessing import FastaStringExtractor, one_hot_encode
import sys
import os

# Constants
SEQUENCE_LENGTH = 393216
TARGETS_TXT_URL = 'https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt'
MODEL_PATH = 'https://tfhub.dev/deepmind/enformer/1'
FASTA_FILE_PATH = '/exports/eddie/scratch/s1616612/Annotation_pipeline_230618/data/genome.fa'

# Load targets and model
df_targets = pd.read_csv(TARGETS_TXT_URL, sep='\t')
model = Enformer(MODEL_PATH)
fasta_extractor = FastaStringExtractor(FASTA_FILE_PATH)

def get_diff_vec(chrom, pos, ref, alt, variant_id):
    """
    Calculate the difference vector between reference and alternate sequences using the Enformer model.
    
    :param chrom: Chromosome
    :param pos: Position
    :param ref: Reference allele
    :param alt: Alternate allele
    :param variant_id: Variant identifier
    :return: DataFrames with differences, absolute differences, and proportional differences.
    """
    variant = kipoiseq.Variant(chrom, pos, ref, alt, id=variant_id)
    interval = Interval(variant.chrom, variant.start, variant.start).resize(SEQUENCE_LENGTH)
    
    seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
    center = interval.center() - interval.start
    
    reference_seq = seq_extractor.extract(interval, [], anchor=center)
    alternate_seq = seq_extractor.extract(interval, [variant], anchor=center)
    
    reference_prediction = model.predict_on_batch(one_hot_encode(reference_seq)[np.newaxis])['human'][0]
    alternate_prediction = model.predict_on_batch(one_hot_encode(alternate_seq)[np.newaxis])['human'][0]
    
    diff = alternate_prediction.mean(axis=0) - reference_prediction.mean(axis=0)
    diff_abs = np.abs(alternate_prediction - reference_prediction).mean(axis=0)
    max_df = np.maximum(alternate_prediction, reference_prediction)
    diff_proportional = (np.abs(alternate_prediction - reference_prediction) / max_df).max(axis=0)
    
    variant_scores = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(df_targets.description, diff))}
    variant_scores_abs = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(df_targets.description, diff_abs))}
    variant_scores_proportional = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(df_targets.description, diff_proportional))}
    
    diff_df = pd.DataFrame([variant_scores])
    diff_df.insert(0, 'id', variant_id)
    
    diff_abs_df = pd.DataFrame([variant_scores_abs])
    diff_abs_df.insert(0, 'id', variant_id)
    
    diff_proportional_df = pd.DataFrame([variant_scores_proportional])
    diff_proportional_df.insert(0, 'id', variant_id)
    
    return diff_df, diff_abs_df, diff_proportional_df

def main():
    if len(sys.argv) != 5:
        print("Usage: script.py <input_path> <output_path> <output_path_abs> <output_path_prop>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    output_path_abs = sys.argv[3]
    output_path_prop = sys.argv[4]

    if not os.path.exists(input_path):
        print(f"Error: Input file '{input_path}' does not exist.")
        sys.exit(1)
    
    data = pd.read_table(input_path, sep='\t')
    
    result = pd.DataFrame()
    result_abs = pd.DataFrame()
    result_proportional = pd.DataFrame()
    
    for index, row in data.iterrows():
        diff_vec, diff_vec_abs, diff_vec_proportional = get_diff_vec(
            row['chrom'], row['pos'], row['ref'], row['alt'], row['id']
        )
        print(f"Processing variant {index + 1}/{len(data)}: {row['id']}")
        
        result = pd.concat([result, diff_vec], ignore_index=True, axis=0)
        result_abs = pd.concat([result_abs, diff_vec_abs], ignore_index=True, axis=0)
        result_proportional = pd.concat([result_proportional, diff_vec_proportional], ignore_index=True, axis=0)
    
    result.to_csv(output_path, index=False, sep='\t')
    result_abs.to_csv(output_path_abs, index=False, sep='\t')
    result_proportional.to_csv(output_path_prop, index=False, sep='\t')
    
    print("Processing complete. Results saved to the specified output files.")

if __name__ == "__main__":
    main()