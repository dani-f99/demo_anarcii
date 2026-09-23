# Imports
from anarcii.output_data_processing import imgt_order, numbered_sequence_dict, required_residue_numbers
from anarcii import Anarcii # https://github.com/oxpig/ANARCII (! pip install anarcii)
from scripts.helpers import translateNT, get_region, get_digit
import pandas as pd
import os


# Bulding the anarcii model for the wanted type of sequence (antibody, tcr, shark or unknown) and instantiate the model. 
anarchii_model = Anarcii(seq_type="antibody", mode="accuracy")


################################
def use_anarchii(nt_seqs: str,
                 model = anarchii_model):
    """
    Helper function that use to anarcii model to number amino acids positions according to the IMGT numbering scheme

    list_seq: array-like -> list of AA sequence to number.
    model: str -> Which model of the anarcii algorithm to use (see doc for other options).
    """
    

    # Call the number method on a list of sequences, path to a fasta or PDB file.
    results = model.number(nt_seqs)

    return results

##########################################################
def aligned_frame(numbered_results: dict) -> pd.DataFrame:
    """
    Reindex ANARCII numbering results onto a shared IMGT position set so every
    sequence gets the same-length alignment, with '-' for unused CDR3 insertions.
    numbered_results: dict -> ANARCII model output from `Anarcii.number()`.
    """
     
    residue_numbers = set(required_residue_numbers)
    rows = {}
    for name, result in numbered_results.items():
        numbering = result["numbering"] or []
        residue_numbers.update(num for num, _ in numbering)
        rows[name] = numbered_sequence_dict(numbering)

    columns = [str(num) + ins.strip() for num, ins in imgt_order(residue_numbers)]
    return pd.DataFrame.from_dict(rows, orient="index", columns=columns).fillna("-")


###################
class AssignCDR3():
    def __init__(self, 
                 input_df:pd.DataFrame | str, 
                 germline_column:str,
                 seq_column:str,
                 results_path:str,
                 append_cdr3: tuple = (True, "cdr3_aa"),
                 model_seqtype:str = "antibody", 
                 model_mode: str = "accuracy"):
        
        """
        Custom class that takes a NT sequence column from dataframe and assign an AnarcII allighen CDR3 to it.

        input_df: pd.DataFrame | str -> Dataset in pd.DataFrame format or exact string path to the input dataset in CSV format.
        germline_column: str -> string name of the column which contains the germline NT sequence.
        seq_column: str -> string name of the column which contains the sequencing NT sequence.
        append_cdr3: tuple with boolean value and string, if True will try to append the cdr3 aa sequence originated from the specified column into the translated sequence.
        model_seqtype:str -> Which model the anarcii algorithm wil load.
        model_mode:str -> Processing mode of the anarcii algorithm.
        """

        # Defining column names
        self.germline_column = germline_column
        self.seq_column = seq_column
        self.append_cdr3 = append_cdr3
        self.results_path = results_path

        # Getting CDR3 information if required
        if append_cdr3[0]:
            self.cdr3_aa_col = append_cdr3[1] 


        # Importing the raw data into python
        if isinstance(input_df, pd.DataFrame):
            self.input_df = input_df

        elif isinstance(input_df, str):
            try:
                self.input_df = pd.read_csv(input_df, index_col=0)

            except:
                raise Exception(f"> Invalid input path (under `input_df` argument): `{input_df}`")

        self.anarchii_model = Anarcii(seq_type=model_seqtype, mode=model_mode)


    def anarchii_allign(self):
        """
        Initializing the class, getting the raw data, and alligining the sequence according to  
        the AnacrII algorithm.
        """
        # Getting the germline NT column -> translating to AA
        germline_nt = self.input_df[self.germline_column]
        germline_aa = germline_nt.apply(translateNT,aa_end=104).str.replace("-","").str.replace("X","")

        if self.append_cdr3[0]:
            germline_aa = germline_aa + self.input_df[self.append_cdr3[1]]

        # Creating anarcii results dataframe
        sequence_anarcii = self.anarchii_model.number(germline_aa)
        self.results_df = pd.DataFrame(sequence_anarcii).T


        # inserting alligned sequence to the results dataframe
        results_anarcii = aligned_frame(sequence_anarcii)
        results_anarcii.columns = pd.MultiIndex.from_tuples([(i, get_region(get_digit(i))) for i in results_anarcii.columns],
                                                            names=["position", "region"])

        results_anarcii.to_csv(os.path.join(self.results_path, "anarcii_sequences.csv"))

        self.results_df.insert(loc=1, 
                               column="germline_anarcii", 
                               value= results_anarcii.astype(str).agg("".join, axis=1))

        self.results_df.insert(loc=2,
                               column="cdr3_anacrii",
                               value=results_anarcii[[i for i in results_anarcii.columns if i[1] == "cdr3"]].astype(str).agg("".join, axis=1))
        
        
        self.results_df.to_csv(os.path.join(self.results_path, "anarcii_output.csv"))

        return self.results_df

        
    def implant_cdr3(self):
        """
        Inserting the AnarcII-alligned CDR3 sequence into the ImmuneDB seq (from aa position 104 onward).
        """
        # seq_id, sample_id, subject_id, clone_id, functional, germline_nt, sequance_nt, germline_aa, cdr3_aa, cdr3_anarcii
        
        final_output = pd.concat([self.input_df[["seq_id", "sample_id", "subject_id", "clone_id", "functional", "copy_number","germline", "sequence", "cdr3_aa"]],
                                      self.results_df[["germline_anarcii", "cdr3_anacrii"]]],
                                      axis=1)
        final_output = final_output[["seq_id", "sample_id", "subject_id", "clone_id", "functional", "copy_number","sequence","germline", "germline_anarcii", "cdr3_aa", "cdr3_anacrii"]]

        final_output.to_csv(os.path.join(self.results_path, "sequences_results.csv"))

        return final_output