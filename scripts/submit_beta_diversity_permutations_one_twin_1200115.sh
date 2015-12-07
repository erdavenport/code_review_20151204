#!/bin/bash

#$ -S /bin/bash
#$ -q regular.q 
#$ -j n
#$ -N beta_div_perms
#$ -M ed379@cornell.edu
#$ -m be
#$ -pe bscb 4
#$ -l h_rt=350:00:00

# script inputs:
# $1 abo.file <- "../results/pABO/abo.table.pABO.txt"
# $2 cov.file <- "../results/pABO/covs.table.pABO.txt"
# $3 taxa.file <- "../results/pABO/taxa.table.pABO.txt"
# $4 bc.file <- "../results/pABO/beta.div.bc.pABO.txt"
# $5 wu.file <- "../results/pABO/beta.div.wu.pABO.txt"
# $6 uu.file <- "../results/pABO/beta.div.uu.pABO.txt"
# $7 n <- 100
# $8 nproc <- 4
# $9 outpath <- "../results/pABO/8_beta_diversity/"

# date
echo $1
d1=$(date +%s)

# Move to scripts folder:
cd Projects/ABO/scripts/

# Make a directory on the local drive to run computations:
mkdir -p /SSD/$USER/$JOB_ID/scripts/OUT

# Copy results folder so files are available:
rsync -a ../results /SSD/$USER/$JOB_ID/

# ls /SSD/$USER/$JOB_ID/

# Switch to that directory to run scripts:
cd /SSD/$USER/$JOB_ID/scripts/

# Run script:
Rscript /home/ed379/Projects/ABO/scripts/beta_diversity_permutations_one_twin_120115ERD.R $1 $2 $3 $4 $5 $6 $7 $8 $9

# Move out files to correct location:
rsync -a /SSD/$USER/$JOB_ID/results/ ~/Projects/ABO/results/ 

rm -r /SSD/$USER/$JOB_ID/

# date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
