#!/usr/bin/python3
from Bio import AlignIO
import getopt
import sys


def usage():
    print("Usage:\t\n\tpython <PathToScript> -i <alignfile> -o <outputfile>")


def rmgaps(infile, outfile):
    a = AlignIO.read(infile, "fasta")
    i = 0
    print(f"The length of originnal alignment is {a.get_alignment_length()}")
    while i//3*3+3 <= a.get_alignment_length():
        if "-" in a[:,i]:
            if i // 3 != 0:
                a = a[:,0:i//3*3]+a[:, i//3*3+3:]
            else:
                a = a[:,3:]
            i = 0
        else:
            i += 1
    print(f"The length of nogaps alignment is {a.get_alignment_length()}")
    AlignIO.write(a, outfile, "fasta")


def getinout():
    opts, args = getopt.getopt(sys.argv[1:], "-i:-o:-h")
    for opt_name, opt_value in opts:
        if opt_name == "-i":
            infile = opt_value
        if opt_name == "-o":
            outfile = opt_value
        if opt_name == "-h":
            usage()
            sys.exit(0)
    return infile, outfile

if __name__ == '__main__':
    infile, outfile = getinout()
    rmgaps(infile, outfile)
    
