# Imports
from scripts.helpers_claude import translateNT
from anarcii import Anarcii # https://github.com/oxpig/ANARCII (! pip install anarcii)
import pandas as pd
import numpy as np


# Bulding the anarcii model for the wanted type of sequence (antibody, tcr, shark or unknown) and instantiate the model.
anarcii_model = Anarcii(seq_type="antibody", mode="accuracy")


######################################
cols_result = ['nt_og', 'seq_og', 'seq_anarcii', 'numbering', 'chain_type', 'score', 'query_start', 'query_end', 'error', 'scheme']


################################
def use_anarcii(nt_seqs: str,
                model = anarcii_model):
    """
    Helper function that uses the anarcii model to number amino acid positions according to the IMGT numbering scheme.

    nt_seqs: str -> AA sequence to number.
    model: Anarcii -> instantiated anarcii model to run the numbering with.
    """

    # Call the number method on a list of sequences, path to a fasta or PDB file.
    results = model.number(nt_seqs)

    return results


################################
def _failed_row(seq: str, reason: str) -> np.ndarray:
    """
    Builds a `cols_result`-shaped row for a sequence that couldn't be processed
    (e.g. spacer count not divisible by 3), so "table" output always has a
    consistent, fixed-length shape regardless of success/failure.
    """

    return np.array(
        [seq, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, reason, np.nan],
        dtype=object,
    )


######################################
def process_sequence(seq: str,
                     anarcii: bool = True,
                     output_type: str = "str",
                     remove_x: bool = True):
    """
    Custom function that processes our ImmuneDB sequence for the Anarcii algorithm.
    (Should work on any NT DNA sequence, it just won't utilize all of the steps.)

    seq: str -> DNA NT sequence in string format.
    anarcii: bool -> if to use the anarcii algorithm so assign IMGT numbering.
    output_type: str -> must be ["str", "table"]. str for string output of sequence, table for information table about processed sequence.
    remove_x: bool -> will remove "X" unknown AAs from the translated sequence (for the anarcii processing).
    """

    # Regular translation of the NT sequence
    if anarcii is False:
        return translateNT(seq)

    # Verifies that the sequence doesn't have spacers that don't divide by 3
    div3_spacers = seq.count("-") % 3
    if div3_spacers != 0:
        if output_type == "table":
            return _failed_row(seq, "Spacer count not divisible by 3")
        return np.nan

    # Removing spacers - will be introduced after the anarcii run
    seq2translate = seq.replace("-", "")

    # Translation & anarcii utilization
    seq_aa = translateNT(seq2translate)  # translated sequence
    dict_remX = {True: "", False: "X"}
    seq_aa_anarcii = use_anarcii(seq_aa.replace("X", dict_remX[remove_x]))  # anarcii numbered sequence (IMGT numbering)

    # In case anarcii fails to provide a sequence
    try:
        seq_anarcii = "".join([i[1] for i in seq_aa_anarcii["Sequence"]["numbering"]])

    except Exception:
        seq_anarcii = np.nan

    if output_type == "str":
        return seq_anarcii

    elif output_type == "table":
        nt_og = seq
        seq_og = translateNT(nt_og)

        df_results = pd.concat([pd.DataFrame({"nt_og": nt_og,
                                              "seq_og": seq_og,
                                              "seq_anarcii": seq_anarcii}, index=[0]),
                                pd.DataFrame(seq_aa_anarcii).T.reset_index(drop=True)], axis=1)

        return df_results.values[0]

    raise ValueError("> Invalid output_type argument. must be `str` or `table`")
