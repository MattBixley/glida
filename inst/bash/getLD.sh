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

# POPFILE=../extdata/1KGPopulations.panel
POPFILE=$1

# Arguments to this script used to define the chromosome and region
CHR=$2
START=$3
END=$4

# input / output files
VCFFILE=Genotype_${CHR}_${START}-${END}.vcf
SAMPLES=$5

# Download 1000 Genomes data
tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${CHR}:${START}-${END} > ${VCFFILE}

# Calculate LD
/home/nickb/GenomeTools/PLINK/plink --vcf ${VCFFILE} --keep-fam ${SAMPLES} --snps-only no-DI --r2 --out ${VCFFILE%.vcf}       


