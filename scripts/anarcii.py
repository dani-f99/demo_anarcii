# Imports
from scripts.helpers import translateNT
from anarcii import Anarcii # https://github.com/oxpig/ANARCII (! pip install anarcii)
import pandas as pd
import numpy as np


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
    
    # Select the type of sequence (antibody, tcr, shark or unknown) and instantiate the model. 
    #model = Anarcii(seq_type="antibody")

    # Call the number method on a list of sequences, path to a fasta or PDB file.
    results = model.number(nt_seqs)

    return results


######################################
# Known bug > anarcii sometimes fail to provide sequence (error, for example index 1)
cols_result = ['nt_og', 'seq_og', 'seq_anarcii', 'numbering', 'chain_type', 'score','query_start', 'query_end', 'error', 'scheme']

def process_sequence(seq: str,
                     anarcii: bool = True,
                     output_type: str = "str",
                     remove_x:bool = True) -> str:
    """
    Custom function that process our ImmuneDB sequence for the processing of the Anarcii algorighm.
    (Should work on any NT DNA sequence, it's just won't utilize all of the steps.)
    Steps:
   
    seq: str -> DNA NT sequence in string format.
    anarcii: bool -> if to use the anarcii algorithm so assign IMGT numbering.
    output_type: str -> must be ["str", "table"]. str for string output of sequence, table for information table about processed sequence.
    remove_x: bool -> will remove "X" unknown AAs from the translated sequence (for the anarcii processing).
    """


    # Regualr translation of the NT sequence
    if anarcii is False:
        return translateNT(seq)

    # Verifies that the sequence dosent have spacers that dosent divide by 3
    else:
        div3_spacers = seq.count("-") % 3
        if div3_spacers != 0:
             return np.nan

        # Removing spacers - will be indroduced after the anarcii run
        seq2translate = seq.replace("-","")              

        # Translation & anarcii utilization
        seq_aa = translateNT(seq2translate)  # translated sequence 
        dict_remX = {True:"", False:"X"}
        seq_aa_anarchii = use_anarchii(seq_aa.replace("X",dict_remX[remove_x])) # anarcii numbered sequence (IMGT numbering)

        # Incase anarcii fails to privide a sequence
        try:
            seq_anarcii = "".join([i[1] for i in seq_aa_anarchii["Sequence"]["numbering"]])

        except:
            seq_anarcii = np.nan

        if output_type in (["str", "table"]):
            l_output = output_type.lower()
            if l_output == "str":
                return seq_anarcii

            elif l_output  == "table":
                nt_og = seq
                seq_og = translateNT(nt_og)

                df_results = pd.concat([pd.DataFrame({"nt_og":nt_og, 
                                                      "seq_og":seq_og, 
                                                      "seq_anarcii":seq_anarcii}, index=[0]),
                                        pd.DataFrame(seq_aa_anarchii).T.reset_index(drop=True)], axis=1)

                
                return df_results.values[0]
                #return {i:j for i,j in zip(cols_result, df_results.values)}
