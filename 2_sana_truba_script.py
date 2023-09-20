import argparse
import os
import csv
import time
from multiprocessing import cpu_count
import subprocess
from socket import gethostname
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from scripts.utils import *
from scripts.pipeline import *

# parse cli
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("-s", default=0, type=int, help="Start index")
parser.add_argument("-e", default=10000, type=int, help="End index")
parser.add_argument("--vcf", help="VCF file to analyze")
parser.add_argument("--output_dir", default="./results",
                    help="Output directory for runtime file, defaults to ./results")
args = parser.parse_args()
START = args.s
END = args.e
VCF_FILE = args.vcf
OUTPUT_DIR = args.output_dir
NUM_CORES = cpu_count()



# global variables
################################################################################
hostname = gethostname()
if hostname == "nazif-pc":
    RNADUPLEX_LOCATION = "/usr/bin/RNAduplex"
else:
    RNADUPLEX_LOCATION = "/truba/home/mtasbas/miniconda3/envs/venv/bin/RNAduplex"

ENERGY_CUTOFF = 5.0
CSV_WRITE_INTERVAL = 10000  # Set the number of steps to write to CSV
# Set the filename for the CSV file
RESULT_CSV_FILE = f"results/sana_results_{START}_{END}.csv"
MIRNA_CSV_FILE = "data/processed/mirbase/mirbase22.csv"


# functions
################################################################################
def run_jobs(job_list, binary_value):

    with ProcessPoolExecutor() as executor:
        future_jobs = [executor.submit(
            rnaduplex_worker_function, *job) for job in job_list]

        results = []

        for future in as_completed(future_jobs):
            result = future.result()
            result_with_binary = result + (binary_value,)
            results.append(result_with_binary)

        with open(RESULT_CSV_FILE, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(results)
    
    return results


def rnaduplex_worker_function(long_sequence, short_sequence, long_identifier, short_identifier):

    # handle the sequence input
    input_sequence = f"{long_sequence}\n{short_sequence}"

    result = subprocess.run(
        [RNADUPLEX_LOCATION, "-e", f"{ENERGY_CUTOFF}", "-s"],
        input=input_sequence,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,  # Capture stderr
        text=True)

    try:
        lines = result.stdout.split("\n")
        first_line = lines[0].split()
        brackets, long_start_end, _, short_start_end, energy = first_line
        long_bracket, short_bracket = brackets.split("&")
        start_long, end_long = map(int, long_start_end.split(","))
        start_short, end_short = map(int, short_start_end.split(","))
        energy_value = float(energy.strip("()"))

        return start_long - 1, end_long, long_bracket, start_short - 1, end_short, short_bracket, energy_value, long_identifier, short_identifier, long_sequence, short_sequence

    # predictions with positive energy prints out "( 5.0163)" with an extra whitespace in the beginning, so strip function adds another element "(" to the array.
    # we try to capture 6 elements with 5 assignments so the interpreter returns ValueError: too many values to unpack (expected 5). This part handles this issue.
    # predictions with negative energy prints out "(-5.0163)" so it works without problem.
    except ValueError as e:
        lines = result.stdout.split("\n")
        first_line = lines[0].split()

        if first_line == []:
            return 0, 1, ".", 0, 1, ".", 0, long_identifier, short_identifier, long_sequence, short_sequence
        
        brackets, long_start_end, _, short_start_end, _, energy = first_line
        long_bracket, short_bracket = brackets.split("&")
        start_long, end_long = map(int, long_start_end.split(","))
        start_short, end_short = map(int, short_start_end.split(","))
        energy_value = float(energy.strip("()"))

        return start_long - 1, end_long, long_bracket, start_short - 1, end_short, short_bracket, energy_value, long_identifier, short_identifier, long_sequence, short_sequence


def prepare_jobs_from_df(df):
    wt_sequences = df["sequence"].to_list()
    mutated_sequences = df["mutated_sequence"].to_list()
    identifiers = df["id"].to_list()
    mirna_identifiers = mirna_df["accession"].to_list()
    mirna_sequences = (mirna_df["sequence"]
                       .apply(reverse_complement_rna_to_dna)
                       .tolist())

    wt_pairs = [(i, j) for i in wt_sequences for j in mirna_sequences]
    mutated_pairs = [(k, l)
                     for k in mutated_sequences for l in mirna_sequences]
    id_pairs = [(m, n) for m in identifiers for n in mirna_identifiers]

    wt_jobs_set = {
        wt_pair + id_pair for wt_pair, id_pair in zip(wt_pairs, id_pairs)
    }
    mutated_jobs_set = {
        mutated_pair + id_pair
        for mutated_pair, id_pair in zip(mutated_pairs, id_pairs)
    }

    return wt_jobs_set, mutated_jobs_set


def handle_target_file():
    if os.path.isfile(RESULT_CSV_FILE):
        # Get the current timestamp
        timestamp = time.strftime("%Y%m%d%H%M%S")
        # Rename the existing CSV file with a ".bak" suffix and a timestamp
        bak_filename = f"{RESULT_CSV_FILE}_{timestamp}.bak"
        os.rename(RESULT_CSV_FILE, bak_filename)


def apply_pipeline(df):
    df = generate_positions_from_id(df)
    df = generate_alignment_string_from_dot_bracket(df)
    df = generate_match_count_columns(df)
    df = generate_ta_sps_columns(df)
    df = generate_mre_sequence_for_vcf(df)
    df = generate_important_sites(df)
    df = generate_mirna_conservation_column(df)
    df = generate_seed_type_columns(df)
    df = generate_mre_au_content_column(df)
    df = generate_au_content_column_for_vcf(df)
    return df
    # df_filtered = filter_columns_for_xgb_prediction(df)
    # return make_predictions_regressor(df, df_filtered)


def split_df_to_num_thread_chunks(df):
    chunk_size = (len(df) + NUM_CORES - 1) // NUM_CORES
    return [df.iloc[i:i+chunk_size] for i in range(0, len(df), chunk_size)]


def make_predictions_regressor(df, df_filtered, model_name="model_with_no_close_proximity_cv"):

    # Create an empty DMatrix model
    dmatrix_model = xgb.DMatrix(df_filtered)

    # Load the pre-trained model
    model = xgb.Booster()
    model.load_model(f"results/{model_name}.json")

    # Make predictions on the df_filtered DataFrame
    predictions = model.predict(dmatrix_model)

    # Append the predictions to the "prediction" column in the df DataFrame
    df["prediction"] = predictions
    df["binary_prediction"] = (df["prediction"] > 0.5).astype(int)
    df.sort_values(["id", "is_mutated"], inplace=True)

    # creating wt and mut dfs
    wt = df[df.is_mutated == 0].reset_index(drop=True)
    mut = df[df.is_mutated == 1].reset_index(drop=True)

    # Calculate the difference between wt and mut predictions
    wt['pred_difference'] = mut['prediction'] - wt['prediction']
    wt['pred_difference_binary'] = mut['binary_prediction'] - \
        wt['binary_prediction']

    # Merge the difference values back to the original df DataFrame
    df = df.merge(
        wt[['id', 'pred_difference', 'pred_difference_binary']], on='id', how='left')

    return df

# importing dfs
################################################################################
df = pd.read_csv(VCF_FILE)
mirna_df = pd.read_csv(MIRNA_CSV_FILE)
df = df[START:END]


# preprocesses
################################################################################
wt_jobs, mutated_jobs = prepare_jobs_from_df(df)
handle_target_file()


# running
################################################################################
wt_result_array = run_jobs(wt_jobs, 0)
mut_result_array = run_jobs(mutated_jobs, 1)


# creating dfs
################################################################################
wt_result_df = pd.DataFrame(wt_result_array)
mut_result_df = pd.DataFrame(mut_result_array)
rnaduplex_results_df = pd.concat([wt_result_df, mut_result_df])

colnames = ["mrna_start", "mrna_end", "mrna_dot_bracket_5to3", "mirna_start", "mirna_end", "mirna_dot_bracket_5to3", "pred_energy", "mutation_id", "mirna_accession", "mrna_sequence", "mirna_sequence", "is_mutated" ]
rnaduplex_results_df.columns = colnames

rnaduplex_results_df["id"] = rnaduplex_results_df["mutation_id"].astype(str) + "_" + rnaduplex_results_df["mirna_accession"].astype(str)
rnaduplex_results_df.drop(columns=["mutation_id"], inplace=True)

df_chunks = split_df_to_num_thread_chunks(rnaduplex_results_df)

# running prediction pipeline
with ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
    results = list(executor.map(apply_pipeline, df_chunks))
    
pipeline_results_df = pd.concat(results)

filtered_df = filter_columns_for_xgb_prediction(pipeline_results_df)
df = make_predictions_regressor(pipeline_results_df, filtered_df)

df.to_csv(f"results/sana_results_{START}_{END}_with_prediction.csv", index=False)

df2 = df[df.pred_difference_binary != 0]
df2.to_csv(f"results/sana_results_{START}_{END}_with_prediction_only_meaningful_results.csv", index=False)