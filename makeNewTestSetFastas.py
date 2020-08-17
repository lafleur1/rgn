import pandas as pd
import numpy as np

OUT = "./hairpin/positives/"
def makefasta(row, outputDir = OUT):
    id = row['protein_id']
    seq = row['protein_seq']
    with open(outputDir + id + ".fasta", "w") as f:
        f.write(">" + str(id) + "\n" + seq + "\n")



if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    data_df = pd.read_csv("coiled_coil_test_set_len_72.csv", sep = "\t")
    data_df.apply(makefasta, axis=1)

