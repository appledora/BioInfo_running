import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


data = pd.read_csv('paired_reads.csv')

my_read = data["Final Read"]
my_read2 = data["Final Read after GA"]

my_read_name = data["name"]

purine = 0
pyrimidine = 0
other = 0
match_index = []
my_name2 = []

a_read2 = []


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
    mixcr_length.append(len(temp))


def process(lines=None):
    ks = ['sequence', 'quality', 'name']
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
mixcr_read_name = []

with open("mixcrEnv/paired_mixcr.txt", "rU") as f1:
    lines1 = []
    for line1 in f1:
        lines1.append(line1.rstrip())
        if len(lines1) == 3:
            record = process(lines1)
            a_read1.append(record["sequence"])
            mixcr_read_name.append(record["name"])
            #print("#", len(a_read1))
            for itt in range (len(my_read)):
                if(record["name"] == my_read_name[itt]):
                    match_index.append(itt)
                    break;
            lines1 = []
            #print("##", len(a_read1))


print(len(my_read))
print(len(a_read1))
print(len(match_index))

count = 0
isMatch = []
GA_score = []
in_mixcr = []
sub_mixcr = []
accuracy_mixcr = []
mixcr_length = []

for i in range(len(my_read)):
    if i not in match_index:
        print(i," NOT IN MATCH INDEX")
        accuracy_mixcr.append("NULL")
        GA_score.append("NULL")
        isMatch.append("NULL")
        mixcr_length.append("NULL")
        a_read2.append("NULL")
    else:
        for j in range(len(a_read1)):
            if(my_read_name[i] == mixcr_read_name[j]):
                alignments = pairwise2.align.globalms(
                    my_read[i], a_read1[j], 1, -1, -1, 0)
                inC, subC = count_ga(alignments[0][0], alignments[0][1])

                accuracy = ((len(alignments[0][0]) -
                             (inC+subC))/len(alignments[0][0]))*100
                accuracy_mixcr.append(accuracy)
                GA_score.append(alignments[0][2])

                match_process(my_read[i], a_read1[j])
                if (my_read[i] == a_read1[j] or my_read2[i] == a_read1[j]):
                    count += 1
                    isMatch.append(True)
                else:
                    isMatch.append(False)


print(len(a_read2))
print(len(GA_score))
data['presto string'] = a_read2
data['presto length'] = mixcr_length
data['Matches'] = isMatch
data['pRESTO vs Code GA_score'] = GA_score
data['Accuracy after aligning with pResto'] = accuracy_mixcr
data.to_csv("paired_reads-mixcr.csv")


print(count)
print(purine+pyrimidine)
print(purine, "------", pyrimidine, "------", other)
