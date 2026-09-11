# Imports
from scripts.anarcii_claude import process_sequence, cols_result
from scripts.helpers_claude import read_json
import pandas as pd
import os

spacer_string = "---------------------------------------------------------------"


# Importing configuration
print(spacer_string+"\n"+"> Importing paths from `config_claude.json` file.")

config_path = os.path.join(os.path.dirname(__file__), "config_claude.json")
config = read_json(config_path)["settings"]
dir_input, dir_output, file_input, file_output = config["dir_input"], config["dir_output"], config["file_input"], config["file_output"]
seq_col, limit_rows = config["seq_col"], config["limit_rows"]

# Defining row limit if in config
try:
    limit_rows = int(limit_rows)
except (TypeError, ValueError):
    limit_rows = None


# Configuring paths of imports
print(spacer_string+"\n"+f"> Loading input file ({file_input}) to python env.")

path_input = os.path.join(dir_input, file_input)
split_fname = file_output.split(".")
path_output = os.path.join(dir_output, f"{split_fname[0]}{limit_rows}.{split_fname[1]}")

# Loading input into memory
df_input = pd.read_csv(path_input, index_col=0)
if isinstance(limit_rows, int) and limit_rows <= df_input.shape[0]:
    df_input = df_input[seq_col].head(limit_rows)
else:
    df_input = df_input[seq_col]

# Applying the function on input file
print(spacer_string+"\n"+f"> Generating output results file.")

anarcii = df_input.apply(process_sequence, anarcii=True, output_type="table")
anarcii_df = pd.DataFrame(anarcii.tolist(), columns=cols_result, index=anarcii.index)
anarcii_df.to_csv(path_output)

print(spacer_string+"\n"+f"> output file saved as {path_output}.")
