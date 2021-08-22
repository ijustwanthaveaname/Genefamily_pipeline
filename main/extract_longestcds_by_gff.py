#------------------------------------------------------------------
#------------------------------------------------------------------
#作者：香波地海
#时间：2019/09/12
#脚本名称：TAIR10.1_sub_longest_cds.py
#脚本功能：根据拟南芥注释文件的features，提取最长编码蛋白序列
#运行模式：python Script.py genome.fasta genome.gff
#返回结果格式：
#   >geneid1
#   ATGTTTGGGAAACCCTGCGATGCTACGCT
#   >geneid2
#   ATGCCCGTAGCTAGCGATCGTAGCTAGCTAGCT
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#举例说明提取过程
#---------------------------------------------------------------------------------------------------

# cds = {}
#如果不是第一个gene信息的话，则处理上个基因里面的多个编码蛋白序列，计算长度提取最长的编码序列；

# cds['gene1'] = {}
# cds['gene1']['pro1'] = [('chr1','+',1,100)]
# cds['gene1']['pro1'] = [('chr1','+',1,100),('chr1','+',150,200),('chr1','+',250,300),('chr1','+',350,400),('chr1','+',450,500),]

# cds['gene1']['pro2'] = [('chr1','+',1,100)]
# cds['gene1']['pro2'] = [('chr1','+',1,100),('chr1','+',120,200),('chr1','+',220,300),('chr1','+',320,400),('chr1','+',420,500),]

# cds['gene2'] = {}
# cds['gene2']['pro3'] = [('chr1','+',1,100)]
# cds['gene2']['pro3'] = [('chr1','+',1,100),('chr1','+',110,200),('chr1','+',210,300),('chr1','+',310,400),('chr1','+',410,500),]

# cds['gene2']['pro4'] = [('chr1','+',1,100)]
# cds['gene2']['pro4'] = [('chr1','+',1,100),('chr1','+',100,200),('chr1','+',200,300),('chr1','+',300,400),('chr1','+',400,500),]

#处理最后一个基因的多个编码蛋白序列，提取最长的蛋白编码序列

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import sys
import re
from Bio import SeqIO


#---------------------------------------------------------------------------------------------------
#写一个函数readgff(),读取基因组注释文件，返回gene、mrna、exon、cds、longest_cds的字典
#---------------------------------------------------------------------------------------------------
def readgff(gff_file):
    '''

    ''' 
    g_id = ''
    gene = {}
    mrna = {}
    exon = {}
    cds = {}
    longest_cds = {}

    with open(gff_file) as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            else:
                line = line.strip()
                (chrs,ref,features,start,end,dot,strand,start_point,note)= line.split('\t')

                if features == 'gene' :
                #-----------------------------------------------------------------------------------
                #在读取下一个基因之前，筛选上一个基因的转录本，构建每个基因对应的最长转录本字典
                #-----------------------------------------------------------------------------------
                    if g_id != '':
                        if cds.get(g_id,None):
                            temp = {}
                            for pro, pos in cds[g_id].items():
                                sum = 0
                                for i in pos:
                                    sum += int(i[-1]) - int(i[-2]) + 1
                                temp[pro] = sum
                            longest_pro = max(temp,key=temp.get)
                            longest_cds[g_id] = cds[g_id][longest_pro]
                    #-------------------------------------------------------------------------------
                    #如果g_id为空值，说明这是注释的第一个基因，那么读取基因ig，构建新的基因字典；如果g_id不是空值，
                    #那么找到上一个基因的最长转录本之后，记录和存储下一个基因的各种信息
                    #-------------------------------------------------------------------------------
                    g_id = note.split(';')[0].replace('ID=','')
                    gene[g_id] = (chrs,strand,start,end)
                    cds[g_id] = {}
                    exon[g_id] = {}
                    mrna[g_id] = {}
    
                if features == 'mRNA':
                    r_id = note.split(';')[0].replace('ID=','')
                    mrna[g_id][r_id] = (chrs,strand,start,end)

                if features == 'exon':
                    r_id = note.split(';')[1].replace('Parent=','')
                    if exon[g_id].get(r_id,None):
                        exon[g_id][r_id].append((chrs,strand,start,end))
                    else:
                        exon[g_id][r_id] = [(chrs,strand,start,end)]

                if features == 'CDS':   
                    c_id = note.split(';')[0].replace('ID=cds-','')
                    if cds[g_id].get(c_id,None):
                        cds[g_id][c_id].append((chrs,strand,start,end))
                    else:
                        cds[g_id][c_id] = [(chrs,strand,start,end)]
                else:
                    continue
    #-----------------------------------------------------------------------------------------------
    #分析最后一个基因的最长转录本
    #-----------------------------------------------------------------------------------------------
    if cds.get(g_id,None):
        temp = {}
        for pro, pos in cds[g_id].items():
            sum = 0
            for i in pos:
                sum += int(i[-1]) - int(i[-2]) + 1
            temp[pro] = sum
        longest_pro = max(temp,key=temp.get)
        longest_cds[g_id] = cds[g_id][longest_pro]
    return gene,mrna,exon,cds,longest_cds

#---------------------------------------------------------------------------------------------------
#定义一个函数用于根据基因位置列表，提取并拼接得到orf的函数sequence(postion)
# longest_cds['gene2'] = [('chr1','+',1,100),('chr1','+',100,200),('chr1','+',200,300),
#                                                       ('chr1','+',300,400),('chr1','+',400,500),]
#---------------------------------------------------------------------------------------------------

def sequence(position):
    '''
    根据函数参数元组position提供的位置参数，提取基因的碱基序列，"+"正链的按照起始位置切片提取，"-"负链的转换成反向互补序列提取，
    最后输出dna字符串
    5'----ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCAT----3'
    3'----TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATCGTA----5'

    '''
    hash = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    scaffold,strand,start,end = position
    if strand == '+':
        dna = Fa[str(scaffold)][int(start)-1:int(end)].upper()
    if strand == '-':
        dna1 = Fa[str(scaffold)][int(start)-1:int(end)][::-1].upper()
        dna = ''.join([hash[i] for i in dna1])
    return dna
#---------------------------------------------------------------------------------------------------
#定义一个函数dna_to_aa()，将核苷酸序列转变成氨基酸序列。输入为dna序列，输出为aa序列。
#---------------------------------------------------------------------------------------------------

def dna_to_aa(dna_seq):
    '''
    定义密码子与氨基酸对应的字典，其中终止密码子用*表示

    '''
    aa_dict = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I','ATA':'I', 
    'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'CCT':'P', 'CCC':'P',
    'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A','TAT':'Y',
    'TAC':'Y',  'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 
    'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C', 'TGG':'W', 'CGT':'R', 'CGC':'R','CGA':'R','CGG':'R', 'AGT':'S', 'AGC':'S', 
    'AGA':'R', 'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'TGA':'*', 'TAA':'*', 'TAG':'*',
    }

    
    #定一个空列表，用来存储翻译完成的氨基酸序列
    aa_seq = []
    for i in range(0,len(dna_seq)-2,3):
        codons = dna_seq[i:i+3]
        if aa_dict.get(codons):
            aa_seq.append(aa_dict[codons])
        else:
            aa_seq.append('@')
    aa = "".join(aa_seq)
    return aa
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------


try:
    genome = sys.argv[1]
    gff = sys.argv[2]
    gene,mrna,exon,cds,longest_cds = readgff(gff)
    Fa = {rec.id:rec.seq for rec in SeqIO.parse(genome,"fasta")}
    print('There are {} genes in the genome.'.format(len(gene)))
    print('There are {} genes transcripted mRNA.'.format(len(mrna)))
    print('There are {} genes transcripted cds.'.format(len(exon)))
    print('There are {} longest proteins produced.'.format(len(longest_cds)))
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#输出最长转录本
# longest_cds['gene2'] = [('chr1','+',1,100),('chr1','+',100,200),('chr1','+',200,300),('chr1','+',300,400),('chr1','+',400,500),]
except IndexError:
    print('please input an gff and fasta file!')
finally:
    #---------------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------------
    with open("longest_cds.faa",'w') as W:
        for name,positions in longest_cds.items():
            seq = ''
            for p in positions:
                seq += sequence(p)
            aa = dna_to_aa(seq)
            W.write(">" + name + "\n")
            W.write(aa + "\n")

    print('The programs of {} have completed!'.format(sys.argv[0]))

#------------------------------------------------------------------#------------------------------------------------------------------
#结束了!
#------------------------------------------------------------------#------------------------------------------------------------------
