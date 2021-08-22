#!/usr/bin/python3
### ------------------------------------
# Exclude overlap records in sorted mergedgff file
### ------------------------------------
from collections import OrderedDict
import sys
import getopt
import os


def gffrecord(gffpath):
    """
    Convert gff file to a OrderedDict that keys are the rows that the 3rd column are "gene"
    """
    with open(gffpath, "rt") as f:
        recs_dict = OrderedDict()
        recs = ""
        first_rec = f.readline()
        recs += first_rec
        rec_name = (first_rec.split("\t")[0], int(first_rec.split("\t")[3]), int(first_rec.split("\t")[4]), first_rec.split("\t")[6])
        for line in f:
            recs_list = line.split("\t")
            if len(recs_list) != 9:
                continue
            elif recs_list[2] != "gene":
                recs += line
            else:
                recs_dict[rec_name] = recs
                recs = line
                rec_name = (recs_list[0], int(recs_list[3]), int(recs_list[4]), recs_list[6])
        recs_dict[rec_name] = recs
    return recs_dict


def filtergff(recs_dict, outputfile):
    outname_list = []
    recskeys = list(recs_dict.keys())
    recskeys = sorted(recskeys, key=lambda x:(x[0], x[3], x[1], x[2]))
    firstname = recskeys[0]
    loc = firstname[0]
    start = firstname[1]
    end = firstname[2]
    strand = firstname[3]
    outname_list.append(firstname)
    for name in recskeys[1:]:
        if name[0] == loc and name[3] == strand:
            if (start <= name[1] <= end or start <= name[2] <= end):
                if (name[2] - name[1]) > (end - start):
                    outname_list[-1] = name
                else:
                    continue
            else:
                loc = name[0]
                start = name[1]
                end = name[2]
                strand = name[3]
                outname_list.append(name)
        else:
            loc = name[0]
            start = name[1]
            end = name[2]
            strand = name[3]
            outname_list.append(name)
    with open(outputfile, "w") as fp:
        for key in outname_list:
            fp.write(recs_dict[key])
    return outname_list


def usage():
    print("Usage\n\tpython /pathToScript -g <gfffile> -o <output>")


def getoptfromsys():
    opts, args = getopt.getopt(sys.argv[1:], "-g:-o:-h")
    for opt_name, opt_value in opts:
        if opt_name == "-g":
            gffpath = opt_value
            if not os.path.isfile(gffpath):
                raise IOError("Can't find gff!Please specify the correct path!")
        if opt_name == "-o":
            outpath = opt_value
        if opt_name == "-h":
            usage()
            sys.exit(0)
    return gffpath, outpath


if __name__ == "__main__":
    gffpath, outpath = getoptfromsys()
    recs_dict = gffrecord(gffpath)
    outname_list = filtergff(recs_dict, outpath)
