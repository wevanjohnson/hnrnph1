
##first you must make a vim script file to do so...

vim create_htseq_qsub.sh

##type in "i" and then paste in your script below 

for i in {1..8}
do
touch CB-${i}_L005.qsub
echo "#!/bin/bash -l

#$ -cwd
#$ -N HTseq_CB-${i}_L005
#$ -j y
#$ -o HTseq_CB-${i}_L005.qlog
#$ -m be
#$ -M neema@bu.edu
#$ -l h_rt=120:00:00


echo \"==========================================================\"
echo \"Starting on       : \$(date)\"
echo \"Running on node   : \$(hostname)\"
echo \"Current directory : \$(pwd)\"
echo \"Current job ID    : \$JOB_ID\"
echo \"Current job name  : \$JOB_NAME\"
echo \"Task index number : \$SGE_TASK_ID\"
echo \"==========================================================\"


module load python2.7/Python-2.7.3_gnu446
module load htseq/0.6.1p1

htseq-count -m union -s no -t exon -i gene_id -q \\
/restricted/projectnb/rufy1/alignment/tophat_andi/thout_CB-${i}_L005/accepted_hits.sam \\
/restricted/projectnb/rufy1/reference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf \\
> /restricted/projectnb/rufy1/alignment/HTseqprac_output/htseq_CB-${i}_L005.txt \\  

echo \"==========================================================\"
echo \"Finishing on       : \$(date)\"
echo \"==========================================================\"

" > CB-${i}_L005.qsub
done

##push ESC and then type ":wq!" to exit
##Next, the script has to be executable to the user, so type the following...

chmod u+x create_htseq_qsub.sh

##Now submit this file to create ALL the .qsub files of your samples at once...

./create_htseq_qsub.sh

##To submit files without running it one by one manually, NOTE: that you are now going to define ${NSLOTS} here. 

for i in {1..8} ##This will change with samples CB-#
do
qsub -P rufy1 -pe single_node 8 CB-${i}_L005.qsub 
done
##To submit files individually, NOTE: that you are now going to define ${NSLOTS} here. 

qsub -P rufy1 -pe single_node 8 CB-#_L005.qsub 




