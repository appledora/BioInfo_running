#!/usr/bin/env bash
AssemblePairs.py align -1 R1.fastq -2 R2.fastq --coord illumina --rc head --outname HeadRC --log HeadRC.log
AssemblePairs.py align -1 R1.fastq -2 R2.fastq --coord illumina --rc tail --outname TailRC --log TailRC.log
AssemblePairs.py align -1 R1.fastq -2 R2.fastq --coord illumina --rc both --outname BothRC --log BothRC.log
AssemblePairs.py align -1 R1.fastq -2 R2.fastq --coord illumina --rc none --outname NoneRC --log NoneRC.log