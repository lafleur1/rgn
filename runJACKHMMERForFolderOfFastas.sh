#!/bin/bash 

#run HHBlits from command line, try to speed up the MSA generation by avoiding opening pipes and etc. q/ subprocess.Popen
#arg 1 is a folder name with the fastas to generate
#arg 2 is the destination folder 
#arg 3 is the fasta database being used
#arg 4 is the text file to output times to run a jackhmmer in 
NUM_ITERATIONS=5
Z_VALUE=1e8
E_VALUE=1e-10
NUM_THREADS=8
FASTAFILES=$1
OUTFOLDER=$2
DB=$3
OUTFILE=$4

#if the times file does not exist, create it 
if test -f $OUTFILE; then
        echo "times file already exists"
else
        echo "creating times txt file"
        touch $OUTFILE
fi


for f in ./$FASTAFILES/*; do
        start_time=$SECONDS
        echo "Processing $f file..."
        echo $OUTFOLDER
        STEM=./$OUTFOLDER/$(basename $f .fasta)
        number=$(basename $f .fasta)
        echo $number
        if test -f $STEM.sto; then
                echo "file already exists"
        else
                echo "starting jackhmmer"
                jackhmmer -N ${NUM_ITERATIONS} -Z ${Z_VALUE} --incE ${E_VALUE} --incdomE ${E_VALUE} -E ${E_VALUE} --domE ${E_VALUE} --cpu ${NUM_THREADS} -o /dev/null -A $STEM.sto --tblout $STEM.tblout $f $DB
                echo "done w/ jackhmmer"
                esl-reformat -o $STEM.a2m a2m $STEM.sto
                echo "done with esl-reformat"
                esl-weight -p --amino --informat a2m -o $STEM.weighted.sto $STEM.a2m
                echo "done with esl-weight"
                esl-alistat --weight --amino --icinfo $STEM.icinfo --cinfo $STEM.cinfo $STEM.weighted.sto > /dev/null
                echo "done with esl-alistat"
                elapsed_time=$(($SECONDS - $start_time))
                echo "$(($elapsed_time/60)) min $(($elapsed_time%60)) sec"    
                echo $number","$elapsed_time >> $OUTFILE
        fi
        
done


