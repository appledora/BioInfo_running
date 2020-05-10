# BioInfo_running
<h1> Summary of Work done so far </h1>
The whole explaination of overall process is here. The work has been divided in **two** modules. In the *first* one,information were extracted from the sequence files. In the *second* part, They were compared with pRESTO output TailRC_assemble-pass.fastq. Datasets for the whole process is added here too. 
     
<h3>MODULE 1</h3>

- Take read1 and read2 from the sequence files. Reverse complement read2.
- Find LCS(Longest Common Subsequence) of the two sequences using Suffix Automata.
- Calculate the region to merge around the LCS. Call it areaX from now on.

![img](https://user-images.githubusercontent.com/36191408/78833865-6b35af80-7a0f-11ea-8dea-6a2e9055b484.jpeg)

- Produce final string in two ways:
    1. Merge areaX normally. In case of mismatches choose base depending on quality scores. -> FinalRead column in CSV
    1. Perform Global Alignment on areaX. Merge based on quality score. -> **FinalReadGA** column in CSV
- Calculate the merge accuracy in areaX. 
    1. In normal case , ![img](http://latex.codecogs.com/svg.latex?Accuracy%3D%5Cfrac%7Blen%28areaX%29-number.of.mismatches%7D%7Blen%28areaX%29%7D%2A100%5C%25) -> Total accuracy in region column in CSV 
    2. in Global Align case, ![img](http://latex.codecogs.com/svg.latex?accuracy%3D%5Cfrac%7Blen%28aligned.region%29-%28indels%2Bsubstitutions%29%7D%7Blen%28aligned.region%29%7D%2A100%5C%25) -> Total accuracy in region after GlobalAlign column in CSV
    In most cases, these two columns have the same value.
    
<h2>MODULE 2</h2> 

- Read strings from pRESTO **TailRC_assembled-pair.fastq**. 
- Compare them with **FinalRead** and **FinalReadGA** column.  -> pRESTO vs Code GA_score column
- Calculate the number of strings that match ( 420 strings in this case). -> **isMatch** column
- Align pRESTO strings  with **FinalRead**. Calculate *indel* and *substitution numbers*.
- Calculate accuracy of alignment, as : 
  ![img](http://latex.codecogs.com/svg.latex?%5Cfrac%7Blen%28whole.aligned.string%29-%28indels%2Bsubstitutions%29%7D%7Blen%28whole.aligned.string%29%2A100%5C%25%7D)  -> Accuracy after aligning with pResto column
