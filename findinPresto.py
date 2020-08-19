import numpy as np
import pandas as pd

data = pd.read_csv('/paired_reads-presto.csv')

# r1_prefix = data["prefix"]
# r2_suffix = data ["suffix"]
# presto_string = data[]

no_match = data.loc[data['Matches'] == False]
# no_match.Prefix = no_match.Prefix.astype(str)
# no_match.Suffix = no_match.Suffix.astype(str)
# no_match["presto string"] = no_match["presto string"].astype(str)
r1_prefix = no_match["Prefix"].tolist()
r2_suffix = no_match["Suffix"].tolist()
presto_string = no_match["presto string"].tolist()

prefix_pos = []
suffix_pos = []

# print(presto_string)
for i in range(len(r1_prefix)):
    print("pRESTO string length => ", len(presto_string[i]) )
    print("prefix found at index => ", presto_string[i].index(r1_prefix[i]))
    print("presto overlap region starts at => ", presto_string[i].index(r1_prefix[i]) + len(r1_prefix))
    print("suffix found at index => ", presto_string[i].index(r2_suffix[i]))
