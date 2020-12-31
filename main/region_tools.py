#!/usr/bin/python3
from dataclasses import dataclass, replace
from operator import attrgetter
from Bio.SeqIO.FastaIO import SimpleFastaParser
import getopt
import sys
"""
用法:1. region_tools.py -u <int> -d <int> -i <blast_tabular> -o extend.txt -g <genome_file>
     2. region_tools.py -f [-s] -i extend.txt -t <region_length> -o filter.txt
"""


@dataclass
class Rec:
    __slots__ = ("q_id", "s_id", "identity", "alignment_length", "mismatches",
                 "gap_openings", "q_start", "q_end", "s_start", "s_end",
                 "e_value", "bit_score")
    q_id: str
    s_id: str
    identity: str
    alignment_length: str
    mismatches: str
    gap_openings: str
    q_start: str
    q_end: str
    s_start: int
    s_end: int
    e_value: str
    bit_score: float


def filter_overlap(sort_data, score=False):
    i = 0
    if sort_data:
        if sort_data[0].s_start <= sort_data[0].s_end:
            start = "s_start"
            end = "s_end"
        else:
            start = "s_end"
            end = "s_start"
        while i <= len(sort_data) - 2:
            if getattr(sort_data[i+1], start) < getattr(sort_data[i], end) and \
                    sort_data[i+1].s_id == sort_data[i].s_id:
                if not score:
                    if abs(sort_data[i+1].s_start - sort_data[i+1].s_end) > \
                            abs(sort_data[i].s_start - sort_data[i].s_end):
                        del sort_data[i]
                    else:
                        del sort_data[i+1]
                else:
                    if sort_data[i+1].bit_score > sort_data[i].bit_score:
                        del sort_data[i]
                    else:
                        del sort_data[i+1]
            else:
                i += 1


class RegionIO:
    def __init__(self, file_path):
        self.file_path = file_path
        self.sort_rf = []
        self.sort_rv = []
        self.recs_forward = []
        self.recs_reverse = []
        with open(self.file_path, "r") as f:
            for rec in f:
                rec = rec.strip().split("\t")
                bit_score = float(rec[11])
                rec_dict = Rec(
                    *rec[0:8],
                    int(rec[8]),
                    int(rec[9]),
                    rec[10],
                    int(bit_score) if bit_score == int(bit_score) else bit_score)
                if rec_dict.s_start <= rec_dict.s_end:
                    self.recs_forward.append(rec_dict)
                else:
                    self.recs_reverse.append(rec_dict)

    def extend_region(self, up=0, down=0, genome_path=""):
        if genome_path == "":
            for n in range(len(self.recs_forward)):
                start = self.recs_forward[n].s_start - up
                end = self.recs_forward[n].s_end + down
                self.recs_forward[n] = \
                    replace(self.recs_forward[n], s_start=start) \
                    if start >= 1 \
                        else replace(self.recs_forward[n], s_start=1)
                self.recs_forward[n] = \
                    replace(self.recs_forward[n], s_end=end) \
                    if start >= 1 \
                        else replace(self.recs_forward[n],
                                     s_end=abs(start)+end+1)
            for n in range(len(self.recs_reverse)):
                start = self.recs_reverse[n].s_end - up
                end = self.recs_reverse[n].s_start + down
                self.recs_reverse[n] = \
                    replace(self.recs_reverse[n], s_end=start) \
                        if start >= 1 \
                        else replace(self.recs_reverse[n], s_end=1)
                self.recs_reverse[n] = \
                    replace(self.recs_reverse[n], s_start=end) \
                        if start >= 1 \
                        else replace(self.recs_reverse[n],
                                     s_start=abs(start)+end+1)
        else:
            id_len = {}
            with open(genome_path,"r") as in_handle:
                for title, seq in SimpleFastaParser(in_handle):
                    id_len[title.split(" ")[0]] = len(seq)
            for n in range(len(self.recs_forward)):
                start = self.recs_forward[n].s_start - up
                end = self.recs_forward[n].s_end + down
                self.recs_forward[n] = \
                    replace(self.recs_forward[n], s_start=start) \
                        if start >= 1 \
                        else replace(self.recs_forward[n], s_start=1)
                if end <= id_len[self.recs_forward[n].s_id]:
                    self.recs_forward[n] = \
                        replace(self.recs_forward[n], s_end=end) \
                            if start >= 1 \
                            else replace(self.recs_forward[n],
                                         s_end=abs(start) + end + 1 if
                                         abs(start) + end + 1 <=
                                         id_len[self.recs_forward[n].s_id] else
                                         id_len[self.recs_forward[n].s_id])
                else:
                    self.recs_forward[n] = \
                        replace(self.recs_forward[n],
                                s_end=id_len[self.recs_forward[n].s_id])
                    self.recs_forward[n] = \
                        replace(self.recs_forward[n],
                                s_start=self.recs_forward[n].s_start
                                             -(end-id_len[self.recs_forward[n].s_id])
                            if self.recs_forward[n].s_start
                                -(end-id_len[self.recs_forward[n].s_id]) >= 1 else 1)
            for n in range(len(self.recs_reverse)):
                start = self.recs_reverse[n].s_end - up
                end = self.recs_reverse[n].s_start + down
                self.recs_reverse[n] = \
                    replace(self.recs_reverse[n], s_end=start) \
                        if start >= 1 \
                        else replace(self.recs_reverse[n], s_end=1)
                if end <= id_len[self.recs_reverse[n].s_id]:
                    self.recs_reverse[n] = \
                        replace(self.recs_reverse[n], s_start=end) \
                            if start >= 1 \
                            else replace(self.recs_reverse[n],
                                         s_start=abs(start) + end + 1 if
                                         abs(start) + end + 1 <=
                                         id_len[self.recs_forward[n].s_id] else
                                         id_len[self.recs_forward[n].s_id])
                else:
                    self.recs_reverse[n] = \
                        replace(self.recs_reverse[n],
                                s_start=id_len[self.recs_reverse[n].s_id])
                    self.recs_reverse[n] = \
                        replace(self.recs_reverse[n],
                                s_end=self.recs_reverse[n].s_end
                                           -(end - id_len[self.recs_reverse[n].s_id])
                            if (self.recs_reverse[n].s_end
                                -(end - id_len[self.recs_reverse[n].s_id])) >= 1 else 1)

    def region_sort(self):
        self.sort_rf = sorted(self.recs_forward,
                              key=attrgetter("s_id", "s_start", "s_end"))
        self.sort_rv = sorted(self.recs_reverse,
                              key=attrgetter("s_id", "s_end", "s_start"))

    def filter_threshold(self, threshold):
        if int(threshold) > 0:
            self.sort_rf = [record for record in self.sort_rf
                         if abs(record.s_start-record.s_end) > threshold]
            self.sort_rv = [record for record in self.sort_rv
                         if abs(record.s_start-record.s_end) > threshold]

    def write(self, output_path):
        with open(output_path, "w") as o:
            for rec in self.sort_rf+self.sort_rv:
                o.write(f"{rec.q_id}\t{rec.s_id}"
                        f"\t{rec.identity}\t{rec.alignment_length}"
                        f"\t{rec.mismatches}\t{rec.gap_openings}"
                        f"\t{rec.q_start}\t{rec.q_end}"
                        f"\t{rec.s_start}\t{rec.s_end}"
                        f"\t{rec.e_value}\t{rec.bit_score}\n")


class StepError(Exception):
    def __init__(self, error):
        self.error = error


if __name__ == "__main__":
    path, up, down, filt, genome, output, threshold, score = "", 0, 0, False, "", "", 0, False
    try:
        opts, args = getopt.getopt(sys.argv[1:], "-i:-g:-o:-u:-d:-f-t:-s")
        for opt_name, opt_value in opts:
            if opt_name == "-i":
                path = opt_value
            if opt_name == "-u":
                up = int(opt_value)
            if opt_name == "-d":
                down = int(opt_value)
            if opt_name == "-f":
                filt = True
            if opt_name == "-t":
                threshold = int(opt_value)
            if opt_name == "-g":
                genome = opt_value
            if opt_name == "-o":
                output = opt_value
            if opt_name == "-s":
                score = True
    except getopt.GetoptError as e:
        for error in e.args:
            print("".join(error))
    results = RegionIO(path)
    if up or down:
        if genome:
            results.extend_region(up=up, down=down, genome_path=genome)
        else:
            results.extend_region(up=up, down=down)
    results.region_sort()
    if filt:
        filter_overlap(results.sort_rf, score)
        filter_overlap(results.sort_rv, score)
    if threshold > 0:
        results.filter_threshold(threshold)
    if output:
        results.write(output)
    else:
        results.write(path+"_region_tools")
