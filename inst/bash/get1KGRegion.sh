#! /bin/bash

#### get1KGRegion.sh ####
#
# Download a defined region from 1000 Genomes
#
#
# Nick Burns
# April 2016

# Arguments to this script used to define the chromosome and region
CHR=$1
START=$2
END=$3
VCFFILE=$4

echo ${VCFFILE}
# Download 1000 Genomes data
tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${CHR}:${START}-${END} > ${VCFFILE}
