import pandas as pd
import sys

def process_input_file(data_path):
    """
    Processes the input file and returns the sorted DataFrame.
    Parameters:
        data_path (str): Path to the input data file.
    Returns:
        pd.DataFrame: A DataFrame containing the processed data.
    """
    try:
        # Load the file
        file = pd.read_table(data_path, header=None, names=['chr', 'start', 'end', 'variant_id'])
        # Sort the DataFrame by chromosome and end position
        file = file.sort_values(['chr', 'end'], ignore_index=True)
        return file
    except Exception as e:
        sys.stderr.write(f"Error processing input file: {e}\n")
        sys.exit(1)

def generate_vep_input(file):
    """
    Generates the VEP input DataFrame from the processed data.
    Parameters:
        file (pd.DataFrame): The processed input data.
    Returns:
        pd.DataFrame: A DataFrame formatted for VEP input.
    """
    vep_input = pd.DataFrame()
    vep_input['chr'] = file['chr'].apply(lambda x: int("".join(filter(str.isdigit, x))))
    vep_input['start'] = file['end']
    vep_input['end'] = file['end']

    new = file['variant_id'].str.split('_', n=4, expand=True)
    vep_input['allele'] = new[2] + '/' + new[3]
    vep_input['strand'] = '+'
    vep_input['identifier'] = file['variant_id']

    return vep_input

def generate_enformer_input(file):
    """
    Generates the Enformer input DataFrame from the processed data.
    Parameters:
        file (pd.DataFrame): The processed input data.
    Returns:
        pd.DataFrame: A DataFrame formatted for Enformer input.
    """
    enformer_input = pd.DataFrame()
    enformer_input['chrom'] = file['chr']
    enformer_input['pos'] = file['end']

    new = file['variant_id'].str.split('_', n=4, expand=True)
    enformer_input['ref'] = new[2]
    enformer_input['alt'] = new[3]
    enformer_input['id'] = file['variant_id']

    return enformer_input

def save_to_csv(df, output_path, header=None, sep='\t'):
    """
    Saves a DataFrame to a CSV file.
    Parameters:
        df (pd.DataFrame): The DataFrame to save.
        output_path (str): Path to the output file.
        header (bool or None): Whether to write column names.
        sep (str): Separator for the CSV file.
    """
    try:
        df.to_csv(output_path, index=False, header=header, sep=sep)
    except Exception as e:
        sys.stderr.write(f"Error saving file to {output_path}: {e}\n")
        sys.exit(1)

def main():
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: python script.py <data_path> <vep_path> <enformer_path>\n")
        sys.exit(1)

    data_path = sys.argv[1]
    vep_path = sys.argv[2]
    enformer_path = sys.argv[3]

    # Process the input data file
    file = process_input_file(data_path)

    # Generate VEP and Enformer inputs
    vep_input = generate_vep_input(file)
    enformer_input = generate_enformer_input(file)

    # Save outputs to CSV files
    save_to_csv(vep_input, vep_path, header=None, sep='\t')
    save_to_csv(enformer_input, enformer_path, header=True, sep='\t')

if __name__ == "__main__":
    main()