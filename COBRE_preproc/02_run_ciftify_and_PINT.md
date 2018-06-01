## Running ciftify_subject_fmri and PINT

So we tried to run the ciftify_subject_fmri step on our queue but it failed in many ways.

So now we will rerun it all on the SCC

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/external/SchizConnect/COBRE/hcp
fmri_basedir=/KIMEL/tigrlab/external/SchizConnect/COBRE/epitome/COBRE

subject_list=`cd ${fmri_basedir}; ls -1d A0*`
run_list=""
for subject in ${subject_list}; do
  lowpassfile=${fmri_basedir}/${subject}/REST/SESS01/func_lowpass.run1.01.nii.gz
  if [ -f ${lowpassfile} ] ; then
    run_list="${run_list} ${subject}"
  fi
done

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/REST/SESS01/func_lowpass.run1.01.nii.gz \
  COBRE_{}_SESS01 \
  REST_01" ::: ${run_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 11 -N  cfyfmri_COBRE1 -

run_list2=""
for subject in ${subject_list}; do
    lowpassfile=${fmri_basedir}/${subject}/REST/SESS01/func_lowpass.run1.02.nii.gz
    if [ -f ${lowpassfile} ] ; then
      run_list2="${run_list2} ${subject}"
    fi
done

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/REST/SESS01/func_lowpass.run1.02.nii.gz \
  COBRE_{}_SESS01 \
  REST_02" ::: ${run_list2} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 11 -N  cfyfmri_COBRE2 -
```
### COBRE_A00021598_SESS01 - everything needs to be rerun..

## Running PINT

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/external/SchizConnect/COBRE/hcp
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12_CMH_20171116

basesubs=`ls -1d ${HCP_DATA}/COBRE_A*`
subject_list=""
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
  bsubject=$(basename $subject)
  pint_summary=${pint_outputdir}/${bsubject}/${bsubject}_summary.csv
  if [ -f ${dtseriesfile} ] ; then
    if [ ! -f ${pint_summary} ] ; then
      subject_list="${subject_list} ${subject}"
    fi
  fi
done

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  {}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
  {}/MNINonLinear/fsaverage_LR32k/{/}.L.midthickness.32k_fs_LR.surf.gii \
  {}/MNINonLinear/fsaverage_LR32k/{/}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{/}/{/}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_COBRE -
```

## Extracting the subcortical covariates

```sh
ssh downie.camhres.ca
module load /projects/edickie/privatemodules/ciftify/201707

## specific files for this bit
HCP_DATA=/external/SchizConnect/COBRE/hcp
SCRIPTS_DIR=/projects/edickie/code/SZ_PINT
ts_outputs=/scratch/edickie/saba_PINT/subcortical_ts/

## go through the file system searching for scans where both ciftify_subject_fmri and PINT has run
basesubs=`ls -1d ${HCP_DATA}/COBRE_A*`
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s0.dtseries.nii
  bsubject=$(basename $subject)
  if [ -f ${dtseriesfile} ] ; then
    echo $subject
    ${SCRIPTS_DIR}/bin/subcortical_extraction_20171121.sh ${bsubject} ${dtseriesfile} ${SCRIPTS_DIR}/templates ${ts_outputs}
  fi
done
```

## build qc pages for all the data

```sh
ssh ogawa.camhres.ca
module load /projects/edickie/privatemodules/ciftify/201707
HCP_DATA=/external/SchizConnect/COBRE/hcp
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data

pint_outputdir=/scratch/edickie/saba_PINT/PINT_pcorr6-6-12_CMH_20171116
qc_dir_base=/scratch/edickie/saba_PINT/qc

## go through the file system searching for scans where both ciftify_subject_fmri and PINT has run
basesubs=`ls -1d ${HCP_DATA}/COBRE_A*`
subject_list=""
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
  bsubject=$(basename $subject)
  if [ -f ${dtseriesfile} ] ; then
      subject_list="${subject_list} ${subject}"
  fi
done

## use parallel to run many fmri qc_pages
parallel -j 16 "echo cifti_vis_fmri snaps --hcp-data-dir ${HCP_DATA} --SmoothingFWHM 8 REST_01 {/}" ::: $subject_list
cifti_vis_fmri index --hcp-data-dir ${HCP_DATA}


subject_list=""
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
  bsubject=$(basename $subject)
  pint_summary=${pint_outputdir}/${bsubject}/${bsubject}_summary.csv
  if [ -f ${dtseriesfile} ] ; then
    if [ -f ${pint_summary} ] ; then
      subject_list="${subject_list} ${subject}"
    fi
  fi
done

## use parallel to run the PINT QC
parallel -j 10 " echo cifti_vis_PINT snaps --qcdir ${qc_dir_base}/qc_REST_01_PINT-pcorr6-6-12 \
--hcp-data-dir {//} \
{}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
{/} ${pint_outputdir}/{/}/{/}_summary.csv" ::: ${subject_list}


```

# 2018-02-27

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/external/SchizConnect/COBRE/hcp
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr10-6-12_CMH_20170227

basesubs=`ls -1d ${HCP_DATA}/COBRE_A*`
subject_list=""
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
  bsubject=$(basename $subject)
  pint_summary=${pint_outputdir}/${bsubject}/${bsubject}_summary.csv
  if [ -f ${dtseriesfile} ] ; then
    if [ ! -f ${pint_summary} ] ; then
      subject_list="${subject_list} ${subject}"
    fi
  fi
done

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 10 --search-radius 6 --padding-radius 12 \
  {}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
  {}/MNINonLinear/fsaverage_LR32k/{/}.L.midthickness.32k_fs_LR.surf.gii \
  {}/MNINonLinear/fsaverage_LR32k/{/}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{/}/{/}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_COBRE -
```
