#!/bin/bash

# This is an example of a script that do all the post-processing ciftify functions using singularity exec

# The three functions are:
# ciftify_clean_img - to denoise and smooth the scans
# ciftify_PINT_vertices - to run PINT
# cifti_vis_PINT - to build PINT qc pages
# extracting the subcortical timeseries
# output dlabel files for the PINT roi locations
# calculate seed correlation maps from the pvertices by network
# extract sulcal depth from the tvertex and ptervex locations

# the inputs are:
#   subject -> the surbect id example: "sub-02"
#   session -> the session (or "none" if no session label is present)
#   func_base -> the bids bold base string example: "task-rest_bold"
#   archive_pipedir -> the top directory above the "ciftify" and "fmriprep" derivaties folder
#   outdir -> the base directory for all the derivaties example: $SCRATCH/bids_outputs/${dataset}/ciftify_PINT
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

# sets three bits for the filenames according to wether a session is present
if [[ ${session} = *"ses"* ]]; then
  bids_func_base=${subject}/${session}/${subject}_${session}_${func_base}
  fmriprep_base=${subject}/${session}/func/${subject}_${session}_${func_base}
  name_of_fmri=${session}_${func_base}
else
  bids_func_base=${subject}/${subject}_${func_base}
  fmriprep_base=${subject}/func/${subject}_${func_base}
  name_of_fmri=${func_base}
fi


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

mkdir -p ${outdir}/ciftify_clean_img/$(dirname ${bids_func_base})

singularity exec \
-H ${sing_home} \
-B ${outdir}:/output \
-B ${archive_pipedir}:/archiveout \
${ciftify_container} ciftify_clean_img \
    --output-file=/output/ciftify_clean_img/${bids_func_base}_desc-clean_bold.dtseries.nii \
    --clean-config=/output/ciftify_clean_img/24MP_8Phys_4GSR.json \
    --confounds-tsv=/archiveout/fmriprep/${fmriprep_base}_confounds.tsv \
    --left-surface=/archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.L.midthickness.32k_fs_LR.surf.gii \
    --right-surface=/archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.R.midthickness.32k_fs_LR.surf.gii \
    /archiveout/ciftify/${subject}/MNINonLinear/Results/${name_of_fmri}/${name_of_fmri}_Atlas_s0.dtseries.nii

# Step 2. Run PINT, this will run with the default radii or 6 6 12
mkdir -p ${outdir}/ciftify_PINT/$(dirname ${bids_func_base})

# for simplicity I have also moved the PINT <input-vertices.csv> file into this folder
# example: cp ~/code/ciftify/ciftify/data/PINT/Yeo7_2011_80verts.csv ${outdir}/ciftify_PINT/


mkdir -p ${outdir}/ciftify_PINT/${subject_outdir}
singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} ciftify_PINT_vertices --pcorr \
      /output/ciftify_clean_img/${bids_func_base}_desc-clean_bold.dtseries.nii \
      /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.L.midthickness.32k_fs_LR.surf.gii \
      /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.R.midthickness.32k_fs_LR.surf.gii \
      /output/ciftify_PINT/Yeo7_2011_80verts.csv \
      /output/ciftify_PINT/${bids_func_base}_desc-clean_bold

# Step 3. Generate PINT QC
singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} cifti_vis_PINT subject \
      --ciftify-work-dir /output/ciftify/ \
      --qcdir /output/ciftify_PINT/qc \
      /output/ciftify_clean_img/${bids_func_base}_desc-clean_bold.dtseries.nii \
      ${subject} \
      /output/ciftify_PINT/${bids_func_base}_desc-clean_bold_summary.csv


# Step 4. Run the cleaning and smoothing script w out smoothing before subcortical extraction
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

singularity exec \
-H ${sing_home} \
-B ${outdir}:/output \
-B ${archive_pipedir}:/archiveout \
${ciftify_container} ciftify_clean_img \
    --output-file=/output/ciftify_clean_img/${bids_func_base}_desc-cleansm0_bold.dtseries.nii \
    --clean-config=/output/ciftify_clean_img/24MP_8Phys_4GSR_sm0.json \
    --confounds-tsv=/archiveout/fmriprep/${fmriprep_base}_confounds.tsv \
    /archiveout/ciftify/${subject}/MNINonLinear/Results/${name_of_fmri}/${name_of_fmri}_Atlas_s0.dtseries.nii


# Step 5. Extract the subcortical timeseries from the cleaned files

# Note that the subcortical templates are generated from the main template using the commands in the templates/README.md from this repo..
# example: cp ~/code/ciftify/ciftify/data/PINT/Yeo7_2011_80verts.csv ${outdir}/ciftify_PINT/

mkdir -p ${outdir}/ciftify_meants/$(dirname ${bids_func_base})

template_dlabels=`cd ${outdir}/ciftify_meants/templates ; ls *.dlabel.nii | sed 's/_atlas.dlabel.nii//g' `

for template_dlabel in ${template_dlabels}; do
  singularity exec \
      -H ${sing_home} \
      -B ${outdir}:/output \
      ${ciftify_container} ciftify_meants \
        --outputcsv /output/ciftify_meants/${bids_func_base}_desc-cleansm0_atlas-${template_dlabel}_timeseries.csv \
        /output/ciftify_clean_img/${bids_func_base}_desc-cleansm0_bold.dtseries.nii \
        /output/ciftify_meants/templates/${template_dlabel}_atlas.dlabel.nii
done


## Step 6: build an dscalar file of the pvertices (starts as a dscalar) - output with network labels
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

## Step 6: build an dlabel file of the pvertices (starts as a dscalar) - output with roi labels

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

## writing the dscalar to a dlabel
## note we for this to run properly, the "Yeo7_2011_80verts_roiidx_LUT.txt" needs to be copied into the folder from the ciftify repo
##
singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    ${ciftify_container} wb_command -cifti-label-import \
    /output/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dscalar.nii \
    /output/ciftify_PINT/Yeo7_2011_80verts_roiidx_LUT.txt \
    /output/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dlabel.nii

## then we need to use the ROI versions to extract the sulcal depth
singularity exec \
    -H ${sing_home} \
    -B ${outdir}:/output \
    -B ${archive_pipedir}:/archiveout \
    ${ciftify_container} ciftify_meants \
      --outputcsv /output/ciftify_PINT/${bids_func_base}_desc-sulc_atlas-${vertextype}ROI_morph.csv \
      /archiveout/ciftify/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.sulc.32k_fs_LR.dscalar.nii \
      /output/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dlabel.nii

## remove the dscalar files becasue they are redundant
rm ${outdir}/ciftify_PINT/${bids_func_base}_desc-${vertextype}ROI_atlas.dscalar.nii
done

## Step 7: build seed correlation maps
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
