#!/usr/bin/python
# -*- coding: utf-8 -*-
# This script needs you input file that include IDs of ncbi genes and outputs the fasta file to stdout.
import sys
from Bio import Entrez
from Bio import SeqIO


def getfasta():
    idlistpath = sys.argv[1]
    Entrez.email = "A.N.Other@example.com"
    with open(idlistpath, "rt") as fp:
        for geneid in fp:
            geneid = geneid.strip()
            handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=geneid)
            seq_record = SeqIO.read(handle, "gb")
            handle.close()
            for feature in seq_record.features:
                if feature.type == "CDS":
                    print(">"+geneid+"\n"+feature.location.extract(seq_record.seq))


if __name__ == "__main__":
    getfasta()
