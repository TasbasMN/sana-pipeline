import pandas as pd


def reverse_complement_rna_to_dna(rna_sequence):
    complement = {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = rna_sequence[::-1]
    return ''.join(complement[base] for base in reverse_seq)


def get_nucleotides_in_interval(chrom, start, end):
    # sourcery skip: extract-method
    file_path = f"data/fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        start_offset = start - 1
        end_offset = end - 1
        num_start_new_lines = start_offset // line_length
        num_end_new_lines = end_offset // line_length
        start_byte_position = byte_position + start_offset + num_start_new_lines
        end_byte_position = byte_position + end_offset + num_end_new_lines
        file.seek(start_byte_position)

        # Read the nucleotides in the interval
        nucleotides = file.read(end_byte_position - start_byte_position + 1)

    # Remove newlines from the nucleotides
    nucleotides = nucleotides.replace('\n', '')

    return nucleotides

def get_nucleotide_at_position(chrom, position):
    file_path = f"data/fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        offset = position - 1
        num_new_lines = offset // line_length
        byte_position = byte_position + offset + num_new_lines
        file.seek(byte_position)

        # Read the nucleotide at the position
        nucleotide = file.read(1)
    return nucleotide


def generate_is_mirna_column(df, grch):
    coords = pd.read_csv(f"data//mirbase22_grch{grch}_coordinates.csv")
    df['is_mirna'] = 0
    df["mirna_accession"] = None
    
    # Iterate over each mutation in the mutations dataframe
    for index, row in df.iterrows():
        mutation_chr = row['chr']
        mutation_start = row['pos']

        # Check if the mutation falls into any of the miRNAs
        matching_rnas = coords[(coords['chr'] == mutation_chr) & (coords['start'] <= mutation_start) & (coords['end'] >= mutation_start)]

        if not matching_rnas.empty:
            # Update the 'is_mirna' column to 1 for the current mutation
            df.at[index, 'is_mirna'] = 1
            df.at[index, 'mirna_accession'] = matching_rnas['mirna_accession'].values[0]
    return df

def import_pyensembl(grch):
    if grch not in [37, 38]:
        raise ValueError("grch must be either 37 or 38")

    ens_release = 75 if grch == 37 else 109
    
    import os
    from pyensembl import EnsemblRelease
    
    os.environ['PYENSEMBL_CACHE_DIR'] = "../data"
    assembly = EnsemblRelease(ens_release)
    assembly.download()
    assembly.index()
    
    return assembly

def generate_transcript_id_and_gene_name_columns(df, grch):
    
    assembly = import_pyensembl(grch)
    
    df['transcript_id'] = df.apply(lambda x: assembly.transcript_ids_at_locus(x['chr'], x['pos']), axis=1)
    df["gene_name"] = df.apply(lambda x: assembly.gene_names_at_locus(x['chr'], x['pos']), axis=1)
    
    return df



