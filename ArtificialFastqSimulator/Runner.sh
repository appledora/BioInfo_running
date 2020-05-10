#!/usr/bin/env bash
AssemblePairs.py align -1 ERR346600_2.fastq -2 ERR346600_1.fastq \
    --coord sra --rc tail --outname M1 --log AP.log 
FilterSeq.py quality -s M1_assemble-pass.fastq -q 20 --outname M1 --log FS.log
MaskPrimers.py score -s M1_quality-pass.fastq -p Greiff2014_VPrimers.fasta \
    --start 4 --mode mask --pf VPRIMER --outname M1-FWD --log MPV.log
MaskPrimers.py score -s M1-FWD_primers-pass.fastq -p Greiff2014_CPrimers.fasta \
    --start 4 --mode cut --revpr --pf CPRIMER --outname M1-REV --log MPC.log
CollapseSeq.py -s M1-REV_primers-pass.fastq -n 20 --inner --uf CPRIMER \
    --cf VPRIMER --act set --outname M1
SplitSeq.py group -s M1_collapse-unique.fastq -f DUPCOUNT --num 2 --outname M1
ParseHeaders.py table -s M1_atleast-2.fastq -f ID DUPCOUNT CPRIMER VPRIMER
ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE
ParseLog.py -l FS.log -f ID QUALITY
ParseLog.py -l MPV.log MPC.log -f ID PRIMER ERROR
