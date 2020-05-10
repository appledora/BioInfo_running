import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


data = pd.read_csv('revr2.csv')

my_read = data["Final Read"]
my_read2 = data["Final Read after GA"]

a_read2 = []
purine = 0
pyrimidine = 0
other = 0


def match_process(s1, s2):
    # print(len(s1),"   ", len(s2))
    temp = ""
    global purine
    global pyrimidine
    global other
    minlen = len(s1)
    if(minlen > len(s2)):
        minlen = len(s2)
    for i in range(minlen):
        if (s1[i] != s2[i]):
            temp = temp + s2[i].lower()
            if (s1[i] == 'A' or s1[i] == 'G'):
                if (s2[i] == 'A' or s2[i] == 'G'):
                    purine += 1
                else:
                    other += 1
            elif (s1[i] == 'T' or s1[i] == 'C'):
                if(s2[i] == 'T' or s2[i] == 'C'):
                    pyrimidine += 1
                else:
                    other += 1

        else:
            temp += s2[i]
    a_read2.append(temp)
    # print(temp)


def process(lines=None):
    ks = ['name', 'sequence','optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


def count_ga(s1, s2):
    # print(s1,"\n",s2)
    indelCount = 0
    subCount = 0
    for i in range(len(s2)):
        if(s1[i] != s2[i]):
            if (s1[i] == '-' or s2[1] == '-'):
                indelCount += 1
            else:
                subCount += 1

    return indelCount, subCount


a_read1 = []


with open("pandaEnv/merged_pandaseq.fastq", "rU") as f1:
    lines1 = []
    for line1 in f1:
        lines1.append(line1.rstrip())
        if len(lines1) == 4:
            record = process(lines1)
            a_read1.append(record["sequence"])
            lines1 = []

print(len(my_read))
print(len(a_read1))
count = 0
isMatch = []
GA_score = []
in_presto = []
sub_presto = []
accuracy_presto = []

for i in range(len(my_read)):
    alignments = pairwise2.align.globalms(my_read[i], a_read1[i], 1, -1, -1, 0)
    inC, subC = count_ga(alignments[0][0], alignments[0][1])
    # print(inC,"----",subC)
    # in_presto.append(inC)
    # sub_presto.append(subC)
    accuracy = ((len(alignments[0][0]) - (inC+subC))/len(alignments[0][0]))*100
    accuracy_presto.append(accuracy)
    GA_score.append(alignments[0][2])
    # print(alignments[0][2])
    # print()
    match_process(my_read[i], a_read1[i])
    if (my_read[i] == a_read1[i] or my_read2[i] == a_read1[i]):
        count += 1
        isMatch.append(True)
    else:
        isMatch.append(False)

# print(GA_score)

# gcount = 0
# for i in range (len(GA_score)):
#     if(GA_score[i] <= 168):
#         gcount += 1

data['panda string'] = a_read2
data['Matches'] = isMatch;
data['panda vs Code GA_score'] = GA_score
data['Accuracy after aligning with panda'] = accuracy_presto
data.to_csv("revr1-panda.csv")
# print(isMatch)
# print(len(a_read2))
# print(data.columns)
# match_process("ABCDEFGH","ABEFPQGS")
print(count)
print(purine+pyrimidine)
print(purine, "------", pyrimidine, "------", other)
