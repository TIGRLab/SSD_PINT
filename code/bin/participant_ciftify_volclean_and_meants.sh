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
archive_pipedir=$4
outdir=$5
sing_home=$6
ciftify_container=$7

# subject="sub-02"
# func_base="task-rhymejudgment_bold"
# outdir=$SCRATCH/bids_outputs/${dataset}/fmriprep_p05
# ciftify_container
# sing_home=$SCRATCH/sing_home/ciftify/$dataset

##
# module load singularity
# module load gnu-parallel/20180322

export ciftify_container

# mkdir -p ${outdir}/ciftify_clean_img/${subject}

# Step 1. Run the cleaning and smoothing script
# note: before calling this script, you need to place a clean_config file into ${outdir}/ciftify_clean_img
#       sample clean_config.json files can be found in https://github.com/edickie/ciftify/tree/master/ciftify/data/cleaning_configs
#    example: cp ~/code/ciftify/ciftify/data/cleaning_configs/24MP_8Phys_4GSR.json ${outdir}/ciftify_clean_img
# In this case I manually edited to be:
# {
#   "--detrend": true,
#   "--standardize": true,
#   "--cf-cols": "X,Y,Z,RotX,RotY,RotZ,CSF,WhiteMatter,GlobalSignal",
#   "--cf-sq-cols": "X,Y,Z,RotX,RotY,RotZ,CSF,WhiteMatter,GlobalSignal",
#   "--cf-td-cols": "X,Y,Z,RotX,RotY,RotZ,CSF,WhiteMatter,GlobalSignal",
#   "--cf-sqtd-cols": "X,Y,Z,RotX,RotY,RotZ,CSF,WhiteMatter,GlobalSignal",
#   "--low-pass": 0.1,
#   "--high-pass": 0.01,
#   "--drop-dummy-TRs": 3,
#   "--smooth-fwhm": 8

if [[ ${session} = *"ses"* ]]; then
    mkdir -p ${outdir}/ciftify_clean_img/${subject}/${session}
    singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} ciftify_clean_img \
        --output-file=/output/ciftify_clean_img/${subject}/${session}/${subject}_${session}_${func_base}_desc-volcleansm8_bold.nii.gz \
        --clean-config=/output/ciftify_clean_img/24MP_8Phys_4GSR.json \
        --confounds-tsv=/archiveout/fmriprep/${subject}/${session}/func/${subject}_${session}_${func_base}_confounds.tsv \
        /archiveout/ciftify/${subject}/MNINonLinear/Results/${session}_${func_base}/${session}_${func_base}.nii.gz
else
  singularity exec \
      mkdir -p ${outdir}/ciftify_clean_img/${subject}
      singularity exec \
      -H ${sing_home} \
      -B ${outdir}:/output \
      -B ${archive_pipedir}:/archiveout \
      ${ciftify_container} ciftify_clean_img \
          --output-file=/output/ciftify_clean_img/${subject}/${subject}_${func_base}_desc-volcleansm8_bold.nii.gz \
          --clean-config=/output/ciftify_clean_img/24MP_8Phys_4GSR.json \
          --confounds-tsv=/archiveout/fmriprep/${subject}/func/${subject}_${func_base}_confounds.tsv \
          /archiveout/ciftify/${subject}/MNINonLinear/Results/${func_base}/${func_base}.nii.gz
fi


# Step 2. Extract the spherical timeseries from the cleaned files


# Note that the subcortical templates are generated from the main template using the commands in the templates/README.md from this repo..
# example: cp ~/code/SZ_PINT/templates/tpl-MNI152NLin6Asym_res-02_desc-6mmYeo780_dseg.nii.gz ${outdir}/ciftify_meants/templates
if [[ ${session} = *"ses"* ]]; then
  bids_func_base=${subject}/${session}/${subject}_${session}_${func_base}
else
  bids_func_base=${subject}/${subject}_${func_base}
fi



singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    ${ciftify_container} ciftify_meants \
      --outputcsv /output/ciftify_meants/${bids_func_base}_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv \
      /output/ciftify_clean_img/${bids_func_base}_desc-volcleansm8_bold.nii.gz \
      /output/ciftify_meants/templates/tpl-MNI152NLin6Asym_res-02_desc-6mmYeo780_dseg.nii.gz
