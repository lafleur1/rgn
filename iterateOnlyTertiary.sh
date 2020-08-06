#!/bin/bash 

#runs R script for all tertiary files, then runs PULCHAR to make it full atom representation to be cleaned up in PyRosetta 
#generates and cleans up after itself 
#run from teh main rgnModified folder 

TERTIARYFOLDER=$1 #folder with tertiary output from RGN
FASTAFOLDER=$2
PDBFOLDER=$3

#test if the pdb folder already exists, if it doesn't, create it
if test -d $PDBFOLDER; then
        echo "pdb folder already exists"
else
        echo "creating times txt file"
        mkdir $PDBFOLDER
fi

for f in $TERTIARYFOLDER*; do
        start_time=$SECONDS
        type=${f##*.}
        if [ $type == "tertiary" ]; then
        	echo "IS TERTIARY "$f
        	#run the R script -> parse tertiary coordinates into simple backbone structure
                name=$(basename $f .tertiary)
                totalFastaName=$FASTAFOLDER$name.fasta
                totalPDB=$PDBFOLDER$name.pdb
                echo "fasta: "$totalFastaName
                echo "pdb: "$totalPDB
                ./tertiary2pdb.R -t $f -f $totalFastaName -p $totalPDB
                #after generating backbone pdb, add full atom with PULCHAR
                ./pulchra304/bin/linux/pulchra $totalPDB
                #delete old backbone only pdb 
                rm $totalPDB
                #find and rename the full atom rep of the pdb 
                pulchraOutput=$(basename $totalPDB .pdb).rebuilt.pdb
                echo "renamed: "$pulchraOutput
                mv $PDBFOLDER$pulchraOutput $totalPDB
        fi
        elapsed_time=$(($SECONDS - $start_time))
        echo "$(($elapsed_time/60)) min $(($elapsed_time%60)) sec"    
       
done 



