#!/bin/bash

working_dir=/external/SchizConnect/COBRE

cd ${working_dir}

mkdir -p epitome/Resting_State
subjects=`cd data/nii; ls -1d A00*`

for subj in ${subjects}; do
    fmri_data=$(find data/nii/${subj} -name Resting_State*.nii.gz)
    parallel "mkdir -p ${working_dir}/epitome/Resting_State/${subj}/RST/SESS01/RUN0{#}" ::: $fmri_data
    parallel "cd ${working_dir}/epitome/Resting_State/${subj}/RST/SESS01/RUN0{#}; ln -s ../../../../../../{} {/}" ::: $fmri_data

    t1=$(find data/nii/${subj} -name Resting_State*.nii.gz)
    parallel "mkdir -p ${working_dir}/epitome/Resting_State/${subj}/RST/SESS01/RUN0{#}" ::: $fmri_data
    parallel "cd ${working_dir}/epitome/Resting_State/${subj}/RST/SESS01/RUN0{#}; ln -s ../../../../../../{} {/}" ::: $fmri_data
    mkdir -p epitome/EXP/${subjname}/REST/SESS01/RUN01/
    mkdir -p epitome/EXP/${subjname}/T1/SESS01/RUN01/
    ln -s ../../../../../../${subj} epitome/Resting_State/${subjname}/REST/SESS01/RUN01/func.nii.gz
    t1=$(find raw/*/struct -name "${subjname}*")
    ln -s ../../../../../../${t1} epitome/Resting_State/${subjname}/T1/SESS01/RUN01/t1.nii.gz
done
