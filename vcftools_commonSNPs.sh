#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p priority                          # Partition to run in
#SBATCH --mem=1000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kianhongkock@g.harvard.edu   # Email to which notifications will be sent

module load gcc/6.2.0
module load vcftools/0.1.15
vcftools --vcf 00-common_all.vcf  --remove-indels --recode --recode-INFO VC --recode-INFO PM --recode-INFO NSF --recode-INFO NSM --recode-INFO NSN --recode-INFO REF --recode-INFO SYN --recode-INFO U3 --recode-INFO U5 --recode-INFO ASS --recode-INFO DSS --recode-INFO INT --recode-INFO R3 --recode-INFO R5 --recode-INFO CAF --out SNPs_caf