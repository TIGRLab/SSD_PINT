#!/bin/bash

# This is an example of a script that will three post-processing ciftify functions using singularity exec

# The three functions are:
# ciftify_clean_img - to denoise and smooth the scans
# ciftify_PINT_vertices - to run PINT
# cifti_vis_PINT - to build PINT qc pages

# the inputs are:
#   subject -> the surbect id example: "sub-02"
#   func_base -> the bids bold base string example: "task-rhymejudgment_bold"
#   outdir -> the base directory for all the derivaties example: $SCRATCH/bids_outputs/${dataset}/fmriprep_p05
#   sing_home -> a ramdom empty folder to bind singularity's /home to example: sing_home=$SCRATCH/sing_home/ciftify/$dataset
#   ciftify_container -> full path to the singularty container with the ciftify env inside

subject=$1
session=$2
func_base=$3
outdir=$4
sing_home=$5
ciftify_container=$6

# subject="sub-02"
# func_base="task-rhymejudgment_bold"
# outdir=$SCRATCH/bids_outputs/${dataset}/fmriprep_p05
# ciftify_container
# sing_home=$SCRATCH/sing_home/ciftify/$dataset

##
# module load singularity
# module load gnu-parallel/20180322

export ciftify_container


# Step 2. Extract the spherical timeseries from the cleaned files


# Note that the subcortical templates are generated from the main template using the commands in the templates/README.md from this repo..
# example: cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz ${outdir}/ciftify_meants/templates
# example: cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order.dlabel.nii ${outdir}/ciftify_meants/templates
if [[ ${session} = *"ses"* ]]; then
  bids_func_base=${subject}/${session}/${subject}_${session}_${func_base}
else
  bids_func_base=${subject}/${subject}_${func_base}
fi



singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    ${ciftify_container} ciftify_meants \
      --outputcsv /output/ciftify_meants/${bids_func_base}_desc-volcleansm8_atlas-Shaefer7N100P_timeseries.csv \
      /output/ciftify_clean_img/${bids_func_base}_desc-volcleansm8_bold.nii.gz \
      /output/ciftify_meants/templates/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz

singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    ${ciftify_container} ciftify_meants \
      --outputcsv /output/ciftify_meants/${bids_func_base}_desc-clean_atlas-Shaefer7N100P_timeseries.csv \
      /output/ciftify_clean_img/${bids_func_base}_desc-clean_bold.dtseries.nii \
      /output/ciftify_meants/templates/Schaefer2018_100Parcels_7Networks_order.dlabel.nii

