#!/usr/bin/python3
from Bio import SeqIO
import sys
import getopt


def usage():
    print("Usage:\npython <script_file> -g <fasta_file> -n <nums>")


def get_genome_and_nums():
    opts, args = getopt.getopt(sys.argv[1:], "-g:-n:-h")
    for opt_name, opt_value in opts:
        try:
            if opt_name == "-g":
                my_genome_file = opt_value
                my_genome_fasta = SeqIO.parse(my_genome_file, "fasta")
            if opt_name == "-n":
                my_nums = opt_value
            if opt_name == "-h":
                usage()
                sys.exit(0)
        except Exception as e:
            print("Error details are :", e)
            print("The error type is :", e.__class__.__name__)
        else:
            if opt_name == "-g":
                print(f"Successfully get {opt_name}!\nYour genome file is {opt_value}!")
    return my_genome_fasta, int(my_nums), my_genome_file


def split_genome(genome_fasta, fold_nums, genome_file):
    record_list = [record for record in genome_fasta]  # 保存序列结果为列表，方便切片操作
    record_num = len(record_list)  # 得到序列总数量
    i = 0
    split_nums = record_num // fold_nums  # 每次输入的序列数
    if record_num <= fold_nums:  # 如果记录数小于分割倍数，直接按记录数一一分割
        fold_nums = record_num
        split_nums = 1
    while i < fold_nums-1:  # 考率到会有剩余序列，分fold-1段和last段输出
        SeqIO.write(record_list[i * split_nums:(i + 1) * split_nums],
                    f"{genome_file}_{i + 1}", "fasta")
        i += 1
    SeqIO.write(record_list[i * split_nums:],  # 输入剩余序列
                f"{genome_file}_{i + 1}", "fasta")


if __name__ == "__main__":
    sequence, fold, file = get_genome_and_nums()
    split_genome(sequence, fold, file)
