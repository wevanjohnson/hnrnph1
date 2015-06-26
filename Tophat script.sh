##first you must make a vim script file to do so...

vim create_tophat_qsub.sh

##type in "i" and then paste in your script below 


for i in {1..8}
do
touch CB-${i}_L005.qsub
echo "#!/bin/bash -l

#$ -cwd
#$ -N CB-${i}_L005
#$ -j y
#$ -o CB-${i}_L005.qlog
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

module load python/2.7.5
module load boost/boost_1_51_0_gnu446
module load samtools/samtools-0.1.19_gnu446
module load bowtie2/2.2.2

tophat -p \${NSLOTS} \\
        --tmp-dir=\${TMPDIR} \\
        -G /restricted/projectnb/rufy1/reference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf \\
        --no-coverage-search \\
        -o /restricted/projectnb/rufy1/alignment/tophat_andi/thout_CB-${i}_L005 \\
        /restricted/projectnb/rufy1/reference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \\
        /restricted/project/rufy1/data/zero/fastq_trimmed/CB-${i}_L005*fastq.gz \\

echo \"==========================================================\"
echo \"Finishing on       : \$(date)\"
echo \"==========================================================\"

" > CB-${i}_L005.qsub
done

##push ESC and then type ":wq!" to exit
##Next, the script has to be executable to the user, so type the following...

chmod u+x create_tophat_qsub.sh

##Now submit this file to create ALL the .qsub files of your samples at once...

./create_tophat_qsub.sh


##To submit files without running it one by one manually, NOTE: that you are now going to define ${NSLOTS} here. 

for i in {9..16} ##This will change with samples CB-#
do
qsub -P rufy1 -pe single_node 8 CB-${i}_L008.qsub 
done

##To submit files individually, NOTE: that you are now going to define ${NSLOTS} here. 
qsub -P rufy1 -pe single_node 8 CB-#_L008.qsub 


