if __name__ == "__main__":

    from scripts.helpers import read_json
    from scripts.anarcii import AssignCDR3
    from scripts.report import write_report
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
    write_report(output_dir, "Run parameter loaded, pipeline starting.", overwrite=True)

    # Initiating the AnarcII class
    write_report(output_dir, "Loading data into the AnarcII algorithm.", overwrite=False)
    print(spacer,"\n> Initiating AnarcII class")
    aclass = AssignCDR3(pd.read_csv(input_file, index_col=0),
                        germline_column=germline_col, 
                        seq_column=sequence_col,
                        results_path=output_dir)

    # Performing the alligment 
    write_report(output_dir, "Performing AnarcII sequence alligment.", overwrite=False)
    print(spacer,"\n> Performing alligment")
    anarci_allign = aclass.anarchii_allign()

    # Unifing data to our original table
    print(spacer,"\n> Unifing Data")
    write_report(output_dir, "Saving output table.", overwrite=False)
    anarci_cdr3 = aclass.implant_cdr3()

    # Printing final messege + showing saved files
    output_path = os.path.join("results", run_name)
    print(spacer,f"\n> Processing done. Output files saved at `{output_path}`.")
    print(f"> output files: {os.listdir(output_path)}")
    write_report(output_dir, "Run completed.", overwrite=False)