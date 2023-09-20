import pandas as pd

def calculate_au_content(sequence):
    au_count = sequence.count(
        'A') + sequence.count('T') + sequence.count('U')
    return None if len(sequence) == 0 else au_count / len(sequence)

def generate_positions_from_id(vcf_df):
    vcf_df['chr'] = vcf_df['id'].str.split('_').str[0]

    vcf_df['start_coordinate'] = vcf_df['id'].str.split('_').str[1].astype(int) - 30 + vcf_df["mrna_start"]
    vcf_df['end_coordinate'] = vcf_df['id'].str.split('_').str[1].astype(int) - 30 + vcf_df["mrna_end"]
    
    return vcf_df

def generate_alignment_string_from_dot_bracket(df):
    full_strings = []
    for _, row in df.iterrows():
        start_string = (row.mirna_start) * "0"
        mid_string = row["mirna_dot_bracket_5to3"].replace(".", "0").replace(")", "1")
        end_string = (len(row.mirna_sequence) - row.mirna_end -1) * "0"
        
        full_string = start_string + mid_string + end_string
        full_strings.append(full_string)

    df["alignment_string"] = full_strings

    return df

def generate_match_count_columns(df):

    def count_ones(str, seed=False):
        return str[1:7].count("1") if seed else str.count("1")

    df["pred_num_basepairs"] = df["alignment_string"].apply(count_ones)
    df["pred_seed_basepairs"] = df["alignment_string"].apply(
        count_ones, seed=True)

    return df

def generate_ta_sps_columns(df):
    # Generate temporary seed column
    df["seed"] = df["mirna_sequence"].str[1:8].str.replace("T", "U")
    # Read ta sps data
    ta_sps_df = pd.read_csv("data/ta_sps.csv", usecols=["seed_8mer", "ta_log10", "sps_mean"])
    ta_sps_df = ta_sps_df.rename(columns={"seed_8mer": "seed"})
    # Merge dataframes on seed column
    df = df.merge(ta_sps_df, on="seed", how="left")
    # Drop temporary column
    df.drop(columns=["seed"], inplace=True)

    return df

def generate_mre_sequence_for_vcf(vcf_df):

    def slice_column(row):
        return row["mrna_sequence"][row["mre_start"]:row["mre_end"]]
    
    # getting mirna length
    vcf_df["mirna_length"] = vcf_df["mirna_sequence"].str.len()

    # using mirna length to figure out mre coordinates
    vcf_df["mre_end"] = vcf_df["mrna_end"] + vcf_df["mirna_start"]
    vcf_df["mre_start"] = vcf_df["mre_end"] - vcf_df["mirna_length"]

    # some start values might be lower than zero, so we need to adjust
    vcf_df["mre_start"] = vcf_df["mre_start"].apply(lambda x: max(x, 0))

    # creating mre sequence column
    vcf_df["mre_region"] = vcf_df.apply(slice_column, axis=1)

    # dropping temp column
    vcf_df.drop(columns=["mirna_length"], inplace=True)
    
    return vcf_df

def generate_important_sites(df):
    df["anchor_a"] = (df["mre_region"].str[-1] == "A").astype(int)
    df["6mer_seed"] = (
        df["alignment_string"].str[1:7].str.count("0") == 0).astype(int)
    df["match_8"] = (df["alignment_string"].str[7] == "1").astype(int)
    df["6mer_seed_1_mismatch"] = (
        df["alignment_string"].str[1:7].str.count("0") == 1).astype(int)

    df["compensatory_site"] = (
        df["alignment_string"].str[12:17].str.count("0") == 0).astype(int)

    df["supplementary_site"] = (
        df["alignment_string"].str[12:16].str.count("0") == 0).astype(int)
    df["supplementary_site_2"] = (
        df["alignment_string"].str[16:21].str.count("0") == 0).astype(int)
    df["empty_seed"] = (
        df["alignment_string"].str[1:8].str.count("1") == 0).astype(int)

    df["9_consecutive_match_anywhere"] = (df["alignment_string"]
                                          .str
                                          .contains("1{" + str(9) + ",}")
                                          .astype(int))

    return df



def generate_mirna_conservation_column(df):
    mirna_df = (pd.read_csv("data/mirbase22.csv", 
                            usecols=["accession", "conservation"])
                            .rename(columns={"accession": "mirna_accession", "conservation": "mirna_conservation"})
                            [["mirna_accession", "mirna_conservation"]])

    df = df.merge(mirna_df, on="mirna_accession", how="left")
    return df

def generate_important_sites(df):
    df["anchor_a"] = (df["mre_region"].str[-1] == "A").astype(int)
    df["6mer_seed"] = (
        df["alignment_string"].str[1:7].str.count("0") == 0).astype(int)
    df["match_8"] = (df["alignment_string"].str[7] == "1").astype(int)
    df["6mer_seed_1_mismatch"] = (
        df["alignment_string"].str[1:7].str.count("0") == 1).astype(int)

    df["compensatory_site"] = (
        df["alignment_string"].str[12:17].str.count("0") == 0).astype(int)

    df["supplementary_site"] = (
        df["alignment_string"].str[12:16].str.count("0") == 0).astype(int)
    df["supplementary_site_2"] = (
        df["alignment_string"].str[16:21].str.count("0") == 0).astype(int)
    df["empty_seed"] = (
        df["alignment_string"].str[1:8].str.count("1") == 0).astype(int)

    df["9_consecutive_match_anywhere"] = (df["alignment_string"]
                                          .str
                                          .contains("1{" + str(9) + ",}")
                                          .astype(int))

    return df


def generate_seed_type_columns(df):
    df['seed_8mer'] = ((df['anchor_a'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_7mer_a1'] = ((df['anchor_a'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 0)).astype(int)
    df['seed_7mer_m8'] = ((df['anchor_a'] == 0) & (df['6mer_seed'] == 1) & (df['match_8'] == 1) & (
        df['supplementary_site'] == 0) & (df['supplementary_site_2'] == 0)).astype(int)
    df['seed_compensatory'] = ((df['compensatory_site'] == 1) & (
        df['6mer_seed_1_mismatch'] == 1) & (df['match_8'] == 1)).astype(int)

    df['seed_clash_2'] = ((df['supplementary_site'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_3'] = ((df['supplementary_site_2'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_4'] = ((df['empty_seed'] == 1) & (
        df['9_consecutive_match_anywhere'] == 1)).astype(int)
    df['seed_clash_5'] = ((df['pred_num_basepairs'] > 10)
                          & (df['6mer_seed'] == 0)).astype(int)

    return df


def generate_mre_au_content_column(df):
    df["mre_au_content"] = df['mre_region'].apply(calculate_au_content)

    return df

def generate_au_content_column_for_vcf(vcf_df):
    
    def calculate_au_content(sequence):
        au_count = sequence.count('A') + sequence.count('T') + sequence.count('U')
        return None if len(sequence) == 0 else au_count / len(sequence)

    vcf_df["local_au_content"] = vcf_df['mrna_sequence'].apply(calculate_au_content)
    
    return vcf_df
