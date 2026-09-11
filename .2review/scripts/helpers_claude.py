# Imports
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
def translateNT(seq:str) -> str:
    """
    Helper function that translate NT DNA sequence to AA sequence.
    seq:str -> string of AA sequence to be translated.
    """

    # Getting correct aa length for translation
    len_seq = len(seq)
    len_aa = int((len_seq - len_seq % 3) / 3)
    translated = []

    # Translating sequence
    for i in range(1, len_aa + 1):
        codon = seq[i*3-3:i*3]

        if codon in list(codon_dict.keys()):
            aa = codon_dict[codon]

        elif codon == "---":
            aa = "-"

        else:
            aa = "X"

        translated.append(aa)

    return "".join(translated)
