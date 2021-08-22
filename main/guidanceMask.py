#!/usr/bin/python3
from Bio import SeqIO
import getopt
import sys


def usage():
    print("Usage:\t\n\tpython /PathToScript/rmGapsInCodonAlign.py -s <scorefile> -a <alignfile> -o <outfile>")


def mask(scfile, alignfile, outfile):
    aligns = SeqIO.parse(alignfile, "fasta")
    alignlist = [align for align in aligns]
    with open(scfile, "rt") as sc:
        for line in sc:
            if line.startswith("#"):
                continue
            recs = line.strip().split("\t")
            row = recs[1]
            col = recs[0]
            score = recs[2]
            if score != "-nan" and float(score) < 0.9:
                alignlist[int(row)-1].seq = alignlist[int(row)-1].seq[:int(col)-1]+"N"+ alignlist[int(row)-1].seq[int(col):]
    SeqIO.write(alignlist, outfile, "fasta")


def getoptions():
    opts, args = getopt.getopt(sys.argv[1:], "-s:-a:-o:-h")
    for opt_name, opt_value in opts:
        if opt_name == "-s":
            scfile = opt_value
        if opt_name == "-a":
            alignfile = opt_value
        if opt_name == "-o":
            outfile = opt_value
        if opt_name == "-h":
            usage()
            sys.exit(0)
    return scfile, alignfile, outfile


if __name__ == "__main__":
    scfile, alignfile, outfile = getoptions()
    mask(scfile, alignfile, outfile)
