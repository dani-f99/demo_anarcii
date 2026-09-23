if __name__ == "__main__":

    from scripts.helpers import read_json
    from scripts.anarcii import AssignCDR3
    import pandas as pd
    import os


    spacer = "-------------------------------------------------------------------------"

    # loading run information from `config.json`
    print(spacer,"\n> Loading configuration from `config.json`")
    config_settings = read_json("config.json")["settings"]
    config_cdr3 = read_json("config.json")["AssignCDR3"]
    input_dir, input_file, run_name = config_settings["dir_input"], config_settings["file_input"], config_settings["run_name"]
    germline_col, sequence_col, cdr3_col = config_cdr3["germline_column"], config_cdr3["sequence_column"], config_cdr3["cdr3_column"]

    # Defining paths and creating folders
    input_file = os.path.join(input_dir, input_file)
    output_dir = os.path.join("results", run_name)
    os.makedirs(output_dir, exist_ok=True)

    # Initiating the AnarcII class
    print(spacer,"\n> Initiating AnarcII class")
    aclass = AssignCDR3(pd.read_csv(input_file, index_col=0),
                            germline_column=germline_col, 
                            seq_column=sequence_col,
                            results_path=output_dir)

    # Performing the alligment 
    print(spacer,"\n> Performing alligment")
    anarci_allign = aclass.anarchii_allign()

    # Unifing data to our original table
    print(spacer,"\n> Unifing Data")
    anarci_cdr3 = aclass.implant_cdr3()

    # Printing final messege + showing saved files
    output_path = os.path.join("results", run_name)
    print(spacer,f"\n> Processing done. Output files saved at `{output_path}`.")
    print(f"> output files: {os.listdir(output_path)}")