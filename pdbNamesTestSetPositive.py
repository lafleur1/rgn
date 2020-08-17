#getting actual pdb to split the chains of
import pandas as pd
import os

testSet = pd.read_csv("coiled_coil_test_set_len_72.csv", sep = "\t")
pdbs = testSet.protein_id.to_list()
pdbs = [name[:-2]+".pdb" for name in pdbs]
print (pdbs)
file = open("testSetPositives.txt", "w")
file.write("\n".join(pdbs))
file.close()

#check existance of all in orig
pdbs = testSet.protein_id.to_list()
found = 0
notfound = 0
badNames = []
for name in pdbs:
    #see if the chain exists now
    print ("ON: ", name)
    if os.path.isfile("./chainsFromDesignedStructures/" + name + ".pdb"):
        found += 1
        os.rename("./chainsFromDesignedStructures/" + name + ".pdb", "./chainsFromDesignedStructures/testSetChains/" + name  + ".pdb")
    else:
        notfound += 1
        badNames.append(name)
print ("found ", found, ' not found ', notfound)
for n in badNames:
    print (n)


seqs = testSet.protein_seq.to_list()
for s in seqs:
    print (len(s))