#!/home/miniconda3/envs/public_py37/bin/python
# -*- coding: utf-8 -*-
from Bio import SeqIO
import sys
import getopt


def exclude_stopcodon(filepath, outpath):
    """
    Main process that can delete stopcodon in the end of sequences or 
    replace the premature termination stop codon with "NNN".
    """
    recs = SeqIO.parse(filepath, "fasta")
    recs_list = []
    for rec in recs:
        if len(rec.seq) % 3 != 0:
            raise(TypeError(
                "The length of your file is not a multiple of 3,you must input codon file!"))
        else:
            for base in set(rec.seq):
                if base.upper() not in ["A", "T", "G", "C", "N"]:
                    raise (
                        TypeError("Your file includes illegal characters, please check it."))
        i = 0
        stopcodon = ["TAA", "TAG", "TGA"]
        while i <= len(rec.seq):
            codon = rec.seq[i:i+3]
            if codon.upper() in stopcodon:
                if i+3 != len(rec.seq):
                    print(
                        f"An premature termination codon was found in {rec.id}and replaced with NNN!")
                    rec.seq = rec.seq[0:i] + \
                        "NNN"+rec.seq[i+3:]
                else:
                    print(
                        f"Found a stop codon in the end of {rec.id} and deleted it.")
                    rec.seq = rec.seq[0:-3]
            i += 3
        recs_list.append(rec)
    SeqIO.write(recs_list, outpath, "fasta")


def get_input_and_output():
    """
    Get the paths of input and output from the command line.
    """
    opts, args = getopt.getopt(sys.argv[1:], "-i:-o:-h")
    for opt_name, opt_value in opts:
        if opt_name == "-i":
            inputfile = opt_value
        if opt_name == "-o":
            outputfile = opt_value
    return inputfile, outputfile


if __name__ == '__main__':
    inputfile, outputfile = get_input_and_output()
    exclude_stopcodon(inputfile, outputfile)
