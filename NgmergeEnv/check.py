import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# this script is used to compare my algorithm with NgMerge's output

data = pd.read_csv('paired_reads.csv')

my_read = data["Final Read"]
my_read2 = data["Final Read after GA"]
my_read_name = data["name"]

a_read2 = []
purine = 0
pyrimidine = 0
other = 0


def match_process(s1, s2):
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
    ng_length.append(len(temp))
    # print(temp)


def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
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
ng_read_name = []


with open("paired_NG.fastq", "rU") as f1:
    lines1 = []
    for line1 in f1:
        lines1.append(line1.rstrip())
        if len(lines1) == 4:
            record = process(lines1)
            a_read1.append(record["sequence"])
            ng_read_name.append(record["name"])
            lines1 = []

print(len(my_read))
print(len(a_read1))

count = 0
countNameMismatch = 0

isMatch = []
GA_score = []
in_ng = []
sub_ng = []
accuracy_ng = []
ng_length = []
missing_id = []

data_names = ''.join(ng_read_name)

for i in my_read_name:
    if data_names.find(str(i)) < 0:
        missing_id.append(i)

print(len(missing_id))
j = 0
for i in range(len(my_read)):
    if(my_read_name[i] not in missing_id):

        alignments = pairwise2.align.globalms(
            my_read[i], a_read1[j], 1, -1, -1, 0)
        inC, subC = count_ga(alignments[0][0], alignments[0][1])
        accuracy = ((len(alignments[0][0]) -
                     (inC+subC))/len(alignments[0][0]))*100
        accuracy_ng.append(accuracy)
        GA_score.append(alignments[0][2])


        match_process(my_read[i], a_read1[j])
        if (my_read[i] == a_read1[j] or my_read2[i] == a_read1[j]):
            count += 1
            isMatch.append(True)
        else:
            isMatch.append(False)
        j+=1;
    else :
        print("MISMATCH")
        ng_length.append("NULL")
        a_read2.append("NULL")
        accuracy_ng.append("NULL")
        GA_score.append("NULL")
        isMatch.append("NULL")
        continue


data['ng string'] = a_read2
data['ng length'] = ng_length
data['Matches'] = isMatch
data['ng vs Code GA_score'] = GA_score
data['Accuracy after aligning with ng'] = accuracy_ng
data.to_csv("paired_reads-ng.csv")

print(count)
# print("===> ", countNameMismatch)
print(purine+pyrimidine)
print(purine, "------", pyrimidine, "------", other)