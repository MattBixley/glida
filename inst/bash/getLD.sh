#! /bin/bash

#### getLD.sh ####
#
# Using 1000 Genomes data, calculate LD    
# for a given subset of people (population-based),
# a given subset of SNPs (variant filtering).
#
# Calculates the pairwise LD for all SNPs
#
# Nick Burns
# April 2016

# Arguments to this script used to define the chromosome and region
VCFFILE=$1
SAMPLES=$2

# Calculate LD
/home/nickb/GenomeTools/PLINK/plink --vcf ${VCFFILE} --keep-fam ${SAMPLES} --snps-only no-DI --r2 --out ${VCFFILE%.vcf}       


