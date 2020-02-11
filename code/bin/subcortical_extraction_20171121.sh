#!/bin/bash

# extracts subcortical ts for input files

subid=${1}
input_dtseries=${2}
templates_dir=${3}
outputdir=${4}


## extract the left thalamus
mkdir ${outputdir}/${subid}
seed_files=`ls ${templates_dir}/*.nii.gz`
parallel -j 4 "ciftify_meants --outputcsv ${outputdir}/${subid}/${subid}_s0_{/.}_meants.csv ${input_dtseries} {}" ::: ${seed_files}
