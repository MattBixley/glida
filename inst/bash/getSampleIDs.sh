#! /bin/bash

# GetSampleIDs.sh
#
# Given a population ($1) extract all relevant sample IDs 
# and write to samples_list.txt

POPFILE=$1
POP=$2
PIPE=$3
OUTFILE=$4

grep ${POP} ${POPFILE} | cut -f1 ${PIPE} ${OUTFILE}

