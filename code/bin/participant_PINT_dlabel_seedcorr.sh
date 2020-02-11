#!/bin/bash

# This is an example of a script that will generate seedcorrelation maps for the last PINT analysis bit..

# The three functions are:
# want to use ciftify_surface_rois to also output a roi dlabel file
# want to used ciftify_surface_rois to output a network cortical rois dlabel file
# want to also the container to calculate the seedcorr for each network label

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

mkdir -p ${outdir}/ciftify_clean_img/${subject}

# Step 1. Run the cleaning and smoothing script
if [[ ${session} = *"ses"* ]]; then
  bids_func_base=${subject}/${session}/${subject}_${session}_${func_base}
else
  bids_func_base=${subject}/${subject}_${func_base}
fi


singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} ciftify_surface_rois \
        --vertex-col "pvertex" \
        --labels-col "NETWORK" \
        --overlap-logic "EXCLUDE" \
        /output/ciftify_PINT/${bids_func_base}_desc-clean_bold_summary.csv \
        6 \
        /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.L.midthickness.32k_fs_LR.surf.gii \
        /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.R.midthickness.32k_fs_LR.surf.gii \
        /output/ciftify_PINT/${bids_func_base}_desc-pvertexNET_atlas.dscalar.nii

singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} ciftify_surface_rois \
        --vertex-col "pvertex" \
        --labels-col "roiidx" \
        --overlap-logic "EXCLUDE" \
        /output/ciftify_PINT/${bids_func_base}_desc-clean_bold_summary.csv \
        6 \
        /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.L.midthickness.32k_fs_LR.surf.gii \
        /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.R.midthickness.32k_fs_LR.surf.gii \
        /output/ciftify_PINT/${bids_func_base}_desc-pvertexROI_atlas.dscalar.nii

for vertextype in tvertex pvertex; do
singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} ciftify_surface_rois \
        --vertex-col ${vertextype} \
        --labels-col "roiidx" \
        --overlap-logic "EXCLUDE" \
        /output/ciftify_PINT/${bids_func_base}_desc-clean_bold_summary.csv \
        6 \
        /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.L.midthickness.32k_fs_LR.surf.gii \
        /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.R.midthickness.32k_fs_LR.surf.gii \
        /output/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dscalar.nii

## note we for this to run properly, the "Yeo7_2011_80verts_roiidx_LUT.txt" needs to be copied into the folder from the ciftify repo
##
singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    ${ciftify_container} wb_command -cifti-label-import \
    /output/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dscalar.nii \
    /output/ciftify_PINT/Yeo7_2011_80verts_roiidx_LUT.txt \
    /output/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dlabel.nii

singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} ciftify_meants \
      --outputcsv /output/ciftify_PINT/${bids_func_base}_desc-sulc_atlas-${vertextype}ROI_morph.csv \
      /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.sulc.32k_fs_LR.dscalar.nii \
      /output/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dlabel.nii

rm ${outdir}/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dscalar.nii
done

## now building the seedcorrelation maps..
for netlabel in 2 3 4 5 6 7; do
  singularity exec \
      -H ${sing_home} \
      -B ${outdir}:/output \
      ${ciftify_container} ciftify_seed_corr --fisher-z \
       --roi-label ${netlabel} \
       --outputname /output/ciftify_PINT/${bids_func_base}_atlas-pvertexNET_roi-${netlabel}_fcmap.dscalar.nii \
        /output/ciftify_clean_img/${bids_func_base}_desc-clean_bold.dtseries.nii \
        /output/ciftify_PINT/${bids_func_base}_desc-pvertexNET_atlas.dscalar.nii
done

## while I'm at is I might as well output sulcal depth..
