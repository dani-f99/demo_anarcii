# Imports
import numpy as np
import regex as re
import json



################################################################
codon_dict = {   'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                 'TCT': "S", 'TCC': "S", 'TCA': "S", 'TCG': "S",
                 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # * for STOP
                 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

                 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

                 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',

                 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                }


################################
# IMGT numbring positions scheme
# source -> https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
imgt_dict = {range(1,27):"fr1",
             range(27,39):"cdr1",
             range(39,56):"fr2",
             range(56,66):"cdr2",
             range(66,105):"fr3",
             range(105,118):"cdr3",
             range(118,129):"fr4"}

######################################################
# Custom function for region extraction per int number
def get_region(position_aa: int, 
               regions_dict: dict = imgt_dict) -> str:
    """
    Function that returns IMGT heavy chain position according to amino acid position input in int format.
    position_aa: int -> amino acid position.
    regions_dict: dictionary -> IMGT numbring scheme in range dict format.
    """
    for key_range, value in regions_dict.items():
        if position_aa in key_range:
            return value
    return np.nan


#####################################
def get_digit(pos_string:str) -> int:
    """
    Custom function that extract amino acid position in int format from string (for example, extract `111` from `111A`)
    pos_string:str -> amino acid position in string format"
    """
    re_found = re.search(pattern = r"[\d]+", 
                         string=pos_string)

    return int(re_found.group(0))


####################################################
def read_json(path:str = "config.json") -> dict:
    """
    Reading information from json file. Used to extract the parameters from the `config.json`.
    path : str -> path of the json file
    """

    with open(path) as config:
        config_f = json.load(config)

    return config_f


################################
def translateNT(seq:str,
                aa_start:int = 1,
                aa_end:int = -1) -> str:
    """
    Helper function that translate NT DNA sequence to AA sequence.
    seq:str -> string of AA sequence to be translated.
    """

    # Getting correct aa length for translation
    if aa_end == -1:
        aa_end = len(seq)

    nt_sequence = seq[aa_start*3 -3 : aa_end*3]
    len_seq = len(nt_sequence)
    len_aa = int((len_seq - len_seq % 3) / 3)
    translated = []

    # Translating sequence
    for i in range(1, len_aa + 1):
        codon = nt_sequence[i*3-3:i*3]

        if codon in list(codon_dict.keys()):
            aa = codon_dict[codon] 

        elif codon == "---":
            aa = "-"
            
        else:
            aa = "X"
         
        translated.append(aa)

    return "".join(translated)