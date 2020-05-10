import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
#global head_r1 , head_r2 , endpos1 ,endpos2 , tail_r1 ,tail_r2 ,tailpos_R1 ,tailpos_R2 ,headpos_R1 ,headpos_R2


class SuffixAutomaton:

    def __init__(self):
        self.states = []
        self.link = []
        self.length = []

    def build(self, string):
        # Creates the automata's initial state
        self.states.append({})
        self.link.append(-1)
        self.length.append(0)
        last = 0
        # Looping over the input string
        for i in range(len(string)):
            # For each character in the string creates a new state
            c = string[i]
            self.states.append({})
            self.length.append(i+1)
            self.link.append(0)
            r = len(self.states) - 1

            p = last
            # Iterate from the last state to the next until there is no more state or no more c-value transitions from p

            while p >= 0 and not self.states[p].get(c):
                # Add a transition from the iterated state to the current state (state created right after the string iteration) with value c

                self.states[p][c] = r
                p = self.link[p]
            if p != -1:
                # q is the last state iterated from p with transition c
                q = self.states[p][c]
                if self.length[p] + 1 == self.length[q]:
                    # If the size is the size of p + 1, the created state r is added as a child of q
                    self.link[r] = q
                else:
                    # Creates a clone named qq of state q, which will be the parent of q. The size will be p.length + 1
                    self.states.append(self.states[q].copy())
                    self.length.append(self.length[p] + 1)
                    self.link.append(self.link[q])
                    qq = len(self.states) - 1

                    self.link[q] = qq
                    self.link[r] = qq

                    # Iterate starting from p until the initial state, as long as p has transition to q with c value
                    while p >= 0 and self.states[p].get(c) == q:
                        # Changes the transitions of p with value c pointing to the state q that was cloned from q
                        self.states[p][c] = qq
                        p = self.link[p]
            last = r
        p = last


def lcs(s, t):
    a = SuffixAutomaton()
    a.build(s)
    v = 0
    l = 0
    best = 0
    bestpos = 0
    for i in range(len(t)):
        while v > 0 and not a.states[v].get(t[i]):
            v = a.link[v]
            l = a.length[v]
        if t[i] in a.states[v]:
            v = a.states[v][t[i]]
            l += 1
        if l > best:
            best = l
            bestpos = i
    start = bestpos - best + 1
    return t[start: start+best], start


def collapse(read1, read2, reg1, reg2, qpos1, qpos2, qq1, qq2, substring, qep1, qep2):
    if (qpos1 > qpos2):
        mergedRead = read1[0:headpos_R1]
    else:
        mergedRead = read2[0: headpos_R2]  # prefix region
    prefix.append(mergedRead)

    # print("merge1 =>", mergedRead)
    qScore_region_R1 = qq1[headpos_R1:headpos_R1 +
                           len(head_r1)+len(substring) + len(tail_r1)]
    qScore_region_R2 = qq2[headpos_R2: headpos_R2 +
                           len(head_r1) + len(substring) + len(tail_r2)]

    for i in range(len(qScore_region_R1)):  # overlap region
        if reg1[i] != reg2[i]:
            if ord(qScore_region_R1[i]) > ord(qScore_region_R2[i]):
                mergedRead += reg1[i]
            else:
                mergedRead += reg2[i]

        else:
            mergedRead += reg1[i]
    # print("merge2 =>", mergedRead)
    if (qpos1 > qpos2):  # suffix region
        remaining = read2[len(reg2):len(r)]
    else:
        remaining = read1[len(reg1): len(firstDNA)]

    suffix.append(remaining)
    mergedRead += remaining
    return mergedRead


def collapse_ga(read1, read2, reg1, reg2, qa1, qa2, a_1, a_2, qpos1, qpos2):
    m = 0
    n = 0
    if (qpos1 > qpos2):
        mergedGA = read1[0:headpos_R1]
    else:
        mergedGA = read2[0: headpos_R2]
    for i in range(len(a_1)):
        if(a_1[i] != a_2[i]):
            if(a_1[i] == '-'):
                n += 1
                mergedGA += a_2[i]
            elif (a_2[i] == '-'):
                m += 1
                mergedGA += a_1[i]
            else:
                if(qa1[m] > qa2[n]):
                    mergedGA += a_1[i]
                else:
                    mergedGA += a_2[i]

        else:
            mergedGA += a_1[i]

    if (qpos1 > qpos2):
        remaining = read2[len(reg2):len(r)]
    else:
        remaining = read1[len(reg1): len(firstDNA)]
    mergedGA += remaining
    return mergedGA


def count_ga(s1, s2):
    # print(s1,"\n",s2)
    indelCount = 0
    subCount = 0
    for i in range(len(s1)):
        if(s1[i] != s2[i]):
            if (s1[i] == '-' or s2[1] == '-'):
                indelCount += 1
            else:
                subCount += 1

    return indelCount, subCount


def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


read1 = []
read2 = []
q1 = []
q2 = []
name = []
FinalRead = []
FinalReadLen = []
LCS = []
LCS_len = []
LCS_Score = []
GA_Score = []
index_R1 = []
index_R2 = []
GA_1 = []
GA_2 = []
accuracy_mm = []
indel_count = []
sub_count = []
FinalReadGA = []
FinalReadGALen = []
accuracy_GA = []
prefix = []
suffix = []

with open("paired_reads1.fq", "rU") as f1, open("paired_reads2.fq", "rU") as f2:
    lines1 = []
    lines2 = []
    for line1, line2 in zip(f1, f2):
        lines1.append(line1.rstrip())
        lines2.append(line2.rstrip())
        if len(lines1) == 4:
            record = process(lines1)
            read1.append(record["sequence"])
            q1.append(record["quality"])  # reverse quality of read 1
            name.append(record["name"])
            lines1 = []
        if len(lines2) == 4:
            record = process(lines2)
            read2.append(record["sequence"])
            q2.append(record["quality"][::-1])
            lines2 = []

head_r1 = None
head_r2 = None
endpos1 = int()
endpos2 = int()
tail_r1 = None
tail_r2 = None
tailpos_R1 = int()
tailpos_R2 = int()
headpos_R1 = 0
headpos_R2 = 0


for i in range(len(read1)):
    print("#", i)
    firstDNA = read1[i]
    secondDNA = read2[i]
    qual1 = q1[i]
    qual2 = q2[i]

    # firstDNA = "abcdefghijklmux"
    # r = "mnjkpqrstu"
    # qual1 = "~!@#$%^&*()_+:<>?"
    # qual2 = "|}{+_&^%$#"

    ReverseSecond = Seq(secondDNA)
    r = ReverseSecond.reverse_complement()

    # ReverseFirst = Seq (rr)
    # firstDNA = ReverseFirst.reverse_complement()

    substring, pos1 = lcs(r, firstDNA)  # finds first occurance in R1
    substring2, pos2 = lcs(firstDNA, r)  # finds first occurance in R2

    index_R1.append(pos1)
    index_R2.append(pos2)

    # print("LCS : ", substring, "\nlength of LCS : ", len(substring))
    # print("Found in Read1 at index => ", pos1)
    # print("Found in Read2 at index => ", pos2)
    LCS.append(substring)
    LCS_len.append(len(substring))
    # if (pos1 != 0 and pos2 != 0):
    if(pos1 > pos2):
        headpos_R1 = pos1 - pos2  # start point of head substring in R1
        headpos_R2 = 0
        head_r1 = firstDNA[headpos_R1: pos1]
        head_r2 = r[0: pos2]

    elif (pos2 >= pos1):
        headpos_R2 = pos2 - pos1  # start point of head substring in R1
        headpos_R1 = 0
        head_r2 = r[headpos_R2: headpos_R2 + pos1]
        head_r1 = firstDNA[0: pos1]

    # to-do : else if (pos 1 == 0 ) / else if (pos2 == 0)
    endpos2 = pos2 + len(substring)  # where the match ends in string t
    endpos1 = pos1 + len(substring)  # do in string s
    # print("ep1: ",endpos1," -- ep2 : ",endpos2, " lenR : ",len(r)-1)
    if (endpos1 <= len(firstDNA)-1 and endpos2 <= len(r)-1):
        if(endpos1 > endpos2):
            # endposition of tail substring in R2
            tailpos_R2 = len(firstDNA) - endpos1  # length of tail in r2
            tail_r1 = firstDNA[endpos1: endpos1+tailpos_R2]  # tail of R1
            tail_r2 = r[endpos2: endpos2+tailpos_R2]  # tail of R2
        elif(endpos2 >= endpos1):
            tailpos_R1 = len(r) - endpos2
            tail_r2 = r[endpos2: endpos2+tailpos_R1]
            tail_r1 = firstDNA[endpos1: endpos1+tailpos_R1]
    else:
        if endpos2 > len(r)-1:
            tailpos_R2 = 0
            tail_r2 = ""
            tail_r1 = ""
        elif endpos1 >= len(firstDNA)-1:
            tailpos_R1 = 0
            tail_r2 = ""
            tail_r1 = ""

    count = 0
    mcount = 0
    # print(tail_r1,"\n",tail_r2)
    for j in range(len(head_r1)):
        if(head_r1[j] != head_r2[j]):
            count = count + 1
        else:
            mcount += 1
    for k in range(len(tail_r1)):
        if(tail_r1[k] != tail_r2[k]):
            count = count + 1
        else:
            mcount += 1
    # print("Mismatch Count in region => ", count)
    # print("Match Count in region => ", mcount)
    # print("final Score => ", len(substring) + mcount - count)
    # print("\n ------------------------------------------ \n")
    finalScore = len(substring) + mcount - count
    LCS_Score.append(finalScore)
    # print(headpos_R1, "----", headpos_R1 +
    #       len(head_r1)+len(substring) + len(tail_r1))
    region_r1 = firstDNA[headpos_R1:headpos_R1 +
                         len(head_r1)+len(substring) + len(tail_r1)]
    region_r2 = r[headpos_R2:headpos_R2 +
                  len(head_r2) + len(substring) + len(tail_r2)]

    finalRead = collapse(firstDNA, r, region_r1, region_r2, pos1, pos2, qual1, qual2,
                         substring, endpos1, endpos2)

    FinalRead.append(finalRead)
    FinalReadLen.append(len(finalRead))
    accuracy = ((len(region_r1) - count) * 100) / len(region_r1)
    # print("Accuracy : ", accuracy, "%")
    accuracy_mm.append(accuracy)
    alignments = pairwise2.align.globalms(
        region_r1, region_r2, 1, -1, 0, 0)
    GA_Score.append(alignments[0][2])
    indelCount, subCount = count_ga(alignments[0][0], alignments[0][1])
    indel_count.append(indelCount)
    sub_count.append(subCount)
    accuracy_ga = (
        ((len(alignments[0][0]) - (indelCount+subCount)) / len(alignments[0][0]))*100)
    accuracy_GA.append(accuracy_ga)
    # print(finalRead)
    finalGARead = collapse_ga(firstDNA, r, region_r1, region_r2,
                              qual1, qual2, alignments[0][0], alignments[0][1], pos1, pos2)
    # print(finalGARead)
    FinalReadGA.append(finalGARead)
    FinalReadGALen.append(len(finalGARead))
    # print(indelCount,".......", subCount)
    # print(format_alignment(*alignments[0]))
    # print(len(name),"---",len(index_R1),"---",len(LCS),"---",len(accuracy_mm))
    # print(alignments[0][0])
    # print(alignments[0][1])
    # print(alignments[0][2])
df = pd.DataFrame({'name': name,
                   'LCS': LCS,
                   'LCS length': LCS_len,
                   'Index_R1': index_R1,
                   'Index_R2': index_R2,
                   'Prefix': prefix,
                   'Suffix': suffix,
                   'Overlap Score': LCS_Score,
                   'Global Align Score in Overlapping area': GA_Score,
                   'InDels in GA': indel_count,
                   'SUbstitutions in GA': sub_count,
                   'Final Read': FinalRead,
                   'Final Read Length': FinalReadLen,
                   'Final Read after GA':  FinalReadGA,
                   'Final Read GA Length': FinalReadGALen,
                   'Total accuracy in region': accuracy_mm,
                   'Total accuracy in region after GA': accuracy_GA})
df.to_csv('paired_reads.csv')
