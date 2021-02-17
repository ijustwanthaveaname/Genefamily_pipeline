#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import getopt
import numpy as np
import pandas as pd
import platform
import sys
try:
    import pysam
except (NameError, ModuleNotFoundError) as e:
    if platform.system() == "Linux":
        print("You'd better download the pysam module, which will speed up the processing.")
    else:
        from Bio import SeqIO


class homogff:
    def __init__(self, genomepath, homopath, gffpath):
        self.genomepath = genomepath
        self.homopath = homopath
        self.gffpath = gffpath
        self.homoID = []
        self.gffDf = None

    def readdata(self):
        """
        read fasta and gff from homopath and gffpath
        """
        # read homo protein fasta
        if "Bio" in sys.modules:
            homofa = SeqIO.parse(self.homopath, "fasta")
            for rec in homofa:
                self.homoID.append(rec.id)
        else:
            fa = pysam.FastaFile(self.homopath)
            self.homoID = fa.references
        try:
            # read gff of the reference genome
            self.gffDf = pd.read_table(self.gffpath,
                                       sep="\t",
                                       comment="#",
                                       header=None)
            self.gffDf.iloc[:, 3:5] = self.gffDf.iloc[:,
                                                      3:5].values.astype(int)
        except pd.errors.ParserError:
            print("Please specfy the correct number of rows to be skipped in gff")

    def proccessdata(self, gffDf, homoID):
        """
        proccess 
        """
        # Check input format
        if not isinstance(gffDf, pd.DataFrame):
            raise TypeError(
                "Please input the gff in DataFrame format read using pandas.")
        if not isinstance(homoID, list):
            raise TypeError("Please input the idlist in the format of list.")
        # get index of cds in ref gff
        psearchPattern = "|".join(homoID)
        features = gffDf.iloc[:, -1]
        protein_id_search = features.str.contains(
            f"protein_id={psearchPattern}")
        index_list = []
        for index, boolvalue in enumerate(protein_id_search):
            if boolvalue == True:
                index_list.append(index)
        homogff = gffDf.iloc[index_list, ]

        # get mRNA id of refgff
        rna_Parent_search = homogff.iloc[:, -
                                         1].str.findall(r"Parent=(.*?)(?:;|$)")
        rna_id_list = []
        for id in rna_Parent_search:
            if id[0] not in rna_id_list:
                rna_id_list.append(id[0])
        rsearchPattern = "|".join(rna_id_list)
        mRNA_Parent_search = features.str.contains(f"Parent={rsearchPattern}")
        mRNA_id_search = features.str.contains(f"ID={rsearchPattern}")
        final_search = mRNA_id_search | protein_id_search | mRNA_Parent_search
        final_gff = gffDf.loc[final_search, :]
        regionID = final_gff.iloc[:, 0].unique()

        # filter region that not contain homo
        filter_list = []
        if "Bio" in sys.modules:
            genomefa = SeqIO.parse(self.genomepath, "fasta")
            for rec in genomefa:
                if rec.id in regionID:
                    filter_list.append(rec)
        else:
            with pysam.FastxFile(self.genomepath) as genomefa:
                for rec in genomefa:
                    if rec.name in regionID:
                        filter_list.append(str(rec))
        return final_gff, filter_list

    @staticmethod
    def writegff(final_gff, writehomogff):
        if not isinstance(final_gff, pd.DataFrame):
            raise TypeError(
                "Please input the gff in DataFrame format read using pandas!")
        final_gff.to_csv(writehomogff,
                         header=None, sep="\t", index=False)

    @staticmethod
    def writefilterref(filter_list, filterregionpath):
        if "Bio" in sys.modules:
            SeqIO.write(filter_list,
                        filterregionpath, "fasta")
        else:
            with open(filterregionpath, "w") as fout:
                fout.write("\n".join(filter_list))


def getoptfromsys():
    opts, args = getopt.getopt(sys.argv[1:], "-g:-h:-a:-b:-f:")
    for opt_name, opt_value in opts:
        if opt_name == "-g":
            genomepath = opt_value
            print(f"Genome is {genomepath}")
        elif opt_name == "-h":
            homologous = opt_value
            print(f"Homologous is {homologous}")
        elif opt_name == "-a":
            annotation = opt_value
            print(f"Annotations is {annotation}")
        elif opt_name == "-b":
            outannotation = opt_value
            print(f"Output homologous gff is {outannotation}")
        elif opt_name == "-f":
            filterregion = opt_value
            print(f"The output file of the Genome regions that include homologous is {filterregion}")
        else:
            raise TypeError(
                "Your must specify options in -g, -h, -a, -b, -f!")
    return genomepath, homologous, annotation, outannotation, filterregion


def mainworkflow():
    genomepath, homologous, annotation, outannotation, filterregion = getoptfromsys()
    # main workflow
    main = homogff(genomepath, homologous, annotation)
    main.readdata()
    final_gff, filter_list = main.proccessdata(main.gffDf, main.homoID)
    main.writegff(final_gff, outannotation)
    main.writefilterref(filter_list, filterregion)


if __name__ == "__main__":
    mainworkflow()
