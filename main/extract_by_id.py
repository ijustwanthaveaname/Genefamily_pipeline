#!/usr/bin/python3
# -*- coding: utf-8 -*-
### ----------------------------
# Extract fasta records by id list.
### ----------------------------
import getopt
import sys
id_list = []
flag = 0


def usage():
    print("Usage:")
    print("\textract_by_id.py -f <fasta_file> -i <id_list>")


def getoptions():
    opts, args = getopt.getopt(sys.argv[1:], "-f:-i:-h")
    for opt_name, opt_value in opts:
        if opt_name == "-f":
            fasta = opt_value
        if opt_name == "-i":
            name =opt_value
        if opt_name == "-h":
            usage()
            sys.exit(0)
    return name, fasta


def extract_fasta(name, fasta):
    with open(name, "r") as f:
        for name in f:
            id_list.append(">"+name.strip())
    with open(fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                if line.strip().split(" ")[0] in id_list:
                    sys.stdout.write(line)
                    flag = 1
                else:
                    flag = 0
            elif flag == 1:
                sys.stdout.write(line)


if __name__ == "__main__":
    try:
        name, fasta = getoptions()
        extract_fasta(name, fasta)
    except Exception as e:
        print(e)
    else:
        print("All done!")
