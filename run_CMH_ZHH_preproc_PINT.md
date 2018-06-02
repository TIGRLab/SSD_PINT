# 2017-10-3

Note: the phenotypic data is in:
/projects/saba/HCP_fMRI_SCZ/new_analyses_20170412/CAMH_ZHH_COBRE/NEWallSubjects_completeData3_DM_not_sexmatched.csv
Dayton has run most of what we need - except for ciftify_subject_fmri and PINT in ZHH

Rerunning ciftify_recon_all because Dayton's version has bad .mat files..

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

export HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/ZHH/hcp
export SUBJECTS_DIR=/KIMEL/tigrlab/external/miklos/freesurfer

subject_list=`cd $SUBJECTS_DIR; ls -1d EXP_2????_SESS01`

cd $HCP_DATA
parallel "echo ciftify_recon_all {}" ::: ${subject_list} | \
  qbatch --walltime 2:00:00 -c 1 -j 1 --ppj 6 -N  cfyfs_ZHH -
```

# Starting with ciftify_subject_fmri

Raw data is in the archive:
Arguments:
	InputfMRI: /external/miklos/epitome/EXP/21368/REST/SESS01//func_lowpass.run1.01.nii.gz
	HCP_DATA: /external/miklos/hcp
	hcpSubject: EXP_21368_SESS01
	NameOffMRI: REST_01
	SmoothingFWHM: 8.0

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/ZHH/hcp
fmri_basedir=/KIMEL/tigrlab/external/miklos/epitome/EXP

subject_list=`cd ${fmri_basedir}; ls -1d 2????`

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}//REST/SESS01/func_lowpass.run1.01.nii.gz \
  EXP_{}_SESS01 \
  REST_01" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  cfyfmri_ZHH -
```

## 2017-11-15 Need to rerun ciftify_subject_fmri for COBRE and CAMH samples

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/SPINS/hcp/
fmri_basedir=/KIMEL/tigrlab/archive/data-2.0/SPINS/pipelines/fmri/rest

subject_list=`cd ${HCP_DATA}; ls -1d SPN01_CMH*`

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/*_lowpass.nii.gz \
  {} \
  REST_01" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  cfyfmri_SPINS -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/ASDD/hcp/
fmri_basedir=/KIMEL/tigrlab/archive/data-2.0/ASDD/pipelines/fmri/rest

subject_list=`cd ${HCP_DATA}; ls -1d ASDD_CMH*`

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/*_lowpass.nii.gz \
  {} \
  REST_01" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  cfyfmri_ASDD -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/RTMSWM/hcp/
fmri_basedir=/KIMEL/tigrlab/archive/data-2.0/RTMSWM/pipelines/fmri/rest

subject_list=`cd ${HCP_DATA}; ls -1d RTMSWM_CMH*`

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/*_lowpass.nii.gz \
  {} \
  REST_01" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  cfyfmri_RTMSWM -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/PNSC/hcp/
fmri_basedir=/KIMEL/tigrlab/archive/data-2.0/PNSC/pipelines/fmri/rest

subject_list=`cd ${HCP_DATA}; ls -1d PNS_CMH*`

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/*_lowpass.nii.gz \
  {} \
  REST_01" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  cfyfmri_PNS -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/DTI3T/hcp/
fmri_basedir=/KIMEL/tigrlab/archive/data-2.0/DTI3T/pipelines/fmri/rest

subject_list=`cd ${HCP_DATA}; ls -1d DTI_CMH*`

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/*_lowpass.nii.gz \
  {} \
  REST_01" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  cfyfmri_DTI -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/ZHH/hcp
fmri_basedir=/KIMEL/tigrlab/external/miklos/epitome/EXP
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d EXP_2????_SESS01`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/data/ZHH/PINT_pcorr6-6-12

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_ZHH -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

BASE_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk


basesubs=`ls -1d ${BASE_DATA}/*/hcp/*_CMH_*`
subject_list=""
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
  if [ -f ${dtseriesfile} ] ; then
    subject_list="${subject_list} ${subject}"
  fi
done

pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12_CMH_20171116

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  {}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
  {}/MNINonLinear/fsaverage_LR32k/{/}.L.midthickness.32k_fs_LR.surf.gii \
  {}/MNINonLinear/fsaverage_LR32k/{/}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{/}/{/}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_CMH -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/DTI3T/hcp/
fmri_basedir=/KIMEL/tigrlab/archive/data-2.0/DTI3T/pipelines/fmri/rest

subject_list=`cd ${HCP_DATA}; ls -1d DTI_CMH*`
for subject in ${subject_list}; do
  dtseriesfile=${HCP_DATA}/${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
  if [ ! -f ${dtseriesfile} ] ; then
    rerun_list="${rerun_list} ${subject}"
    rm -r ${HCP_DATA}/${subject}/MNINonLinear/Results/REST_01
  fi
done


cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/*_lowpass.nii.gz \
  {} \
  REST_01" ::: ${rerun_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  rcfyfmri_DTI -
```

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/SPINS/hcp/
fmri_basedir=/KIMEL/tigrlab/archive/data-2.0/SPINS/pipelines/fmri/rest

subject_list=`cd ${HCP_DATA}; ls -1d SPN01_CMH*`
rerun_list=""
for subject in ${subject_list}; do
  dtseriesfile=${HCP_DATA}/${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
  if [ ! -f ${dtseriesfile} ] ; then
    rerun_list="${rerun_list} ${subject}"
    rm -r ${HCP_DATA}/${subject}/MNINonLinear/Results/REST_01
  fi
done

cd $HCP_DATA

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}/*_lowpass.nii.gz \
  {} \
  REST_01" ::: ${rerun_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  rcfyfmri_SPINS -
```

## 2017-11-17 Rerunning PINT subjects who failed

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

BASE_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12_CMH_20171116

basesubs=`ls -1d ${BASE_DATA}/*/hcp/*_CMH_*`
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
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_CMH -
```

```sh
export HCP_DATA=/scratch/edickie/saba_PINT/data/PNSC/hcp
subject_list=`cd ${HCP_DATA}; ls -1d PNS_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_fmri snaps --qcdir ${PWD}/qc_PNS_fmri --SmoothingFWHM 8 REST_01 ${subject}
done
```
## WIP - Do all the QC

```sh
ssh deckard.camhres.ca
module load /projects/edickie/privatemodules/ciftify/201707
BASE_DATA=/scratch/edickie/saba_PINT/data/
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data

pint_outputdir=/scratch/edickie/saba_PINT/PINT_pcorr6-6-12_CMH_20171116
qc_dir_base=/scratch/edickie/saba_PINT/qc

## go through the file system searching for scans where both ciftify_subject_fmri and PINT has run
basesubs=`ls -1d ${BASE_DATA}/*/hcp/*_CMH_*`
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

## use parallel to run many fmri qc_pages
parallel -j 20 "echo cifti_vis_fmri snaps --qcdir ${qc_dir_base}/qc_REST_01_fmri --hcp-data-dir {//} --SmoothingFWHM 8 REST_01 {/}" ::: $subject_list
cifti_vis_fmri index --qcdir /scratch/edickie/saba_PINT/qc/qc_REST_01_fmri --hcp-data-dir /scratch/edickie/saba_PINT/data//SPINS/hcp

## use parallel to run the PINT QC
parallel -j 20 " echo cifti_vis_PINT snaps --qcdir ${qc_dir_base}/qc_REST_01_PINT-pcorr6-6-12 \
--hcp-data-dir {//} \
{}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
{/} ${pint_outputdir}/{/}/{/}_summary.csv" ::: ${subject_list}


```

# 2017-11-18 Extracting the seed timeseries from the striatum

subid=${1}
input_dtseries=${2}
templates_dir=${3}
outputdir=${4}

```sh
ssh deckard.camhres.ca
module load /projects/edickie/privatemodules/ciftify/201707

## specific files for this bit
BASE_DATA=/scratch/edickie/saba_PINT/data/
SCRIPTS_DIR=/projects/edickie/code/SZ_PINT
ts_outputs=/scratch/edickie/saba_PINT/subcortical_ts/

## go through the file system searching for scans where both ciftify_subject_fmri and PINT has run
basesubs=`ls -1d ${BASE_DATA}/*/hcp/*_???_*`
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s0.dtseries.nii
  bsubject=$(basename $subject)
  if [ -f ${dtseriesfile} ] ; then
    echo $subject
    ${SCRIPTS_DIR}/bin/subcortical_extraction_20171121.sh ${bsubject} ${dtseriesfile} ${SCRIPTS_DIR}/templates ${ts_outputs}
  fi
done
```

## 2017-12-15 Extracting from the ZHH data

```sh
ssh deckard.camhres.ca
module load /projects/edickie/privatemodules/ciftify/201707

## specific files for this bit
BASE_DATA=/scratch/edickie/saba_PINT/data/
SCRIPTS_DIR=/projects/edickie/code/SZ_PINT
ts_outputs=/scratch/edickie/saba_PINT/subcortical_ts/

## go through the file system searching for scans where both ciftify_subject_fmri and PINT has run
basesubs=`ls -1d ${BASE_DATA}/ZHH/hcp/EXP*`
for subject in ${basesubs}; do
  dtseriesfile=${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s0.dtseries.nii
  bsubject=$(basename $subject)
  if [ -f ${dtseriesfile} ] ; then
    echo $subject
    ${SCRIPTS_DIR}/bin/subcortical_extraction_20171121.sh ${bsubject} ${dtseriesfile} ${SCRIPTS_DIR}/templates ${ts_outputs}
  fi
done
```


## errors that were thrown during qc file generation

subs to double check.. (maybe final smoothing failed??)
/scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_H164_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
/scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_H169_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
/scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_H170_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii
/scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_H171_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii  /scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_H172_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii   /scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_S143_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii  /scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_S149_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii  /scratch/edickie/saba_PINT/data//DTI3T/hcp/DTI_CMH_S149_01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii



################################ this is the old stuff.. ##################
ls: cannot access 'EXP_21456_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21457_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21486_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21488_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21489_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21505_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21507_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21532_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21534_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21538_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21543_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21546_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21547_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory
ls: cannot access 'EXP_21550_SESS01/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii': No such file or directory



```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/ZHH/hcp
fmri_basedir=/KIMEL/tigrlab/external/miklos/epitome/EXP

subject_list="21456
21457
21486
21488
21489
21505
21507
21532
21534
21538
21543
21546
21547
21550
"

cd $HCP_DATA

parallel "rm -r EXP_{}_SESS01/MNINonLinear/Results/REST_01" ::: ${subject_list}

parallel "echo ciftify_subject_fmri \
  --SmoothingFWHM 8 \
  --hcp-data-dir ${HCP_DATA} \
  ${fmri_basedir}/{}//REST/SESS01/func_lowpass.run1.01.nii.gz \
  EXP_{}_SESS01 \
  REST_01" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 7 -N  cfyfmri2_ZHH -
```



cifti_vis_recon_all index
Freesurfer fails
COBRE_A00028806_SESS01 - bad brain masking..
EXP_21768_SESS01 - poor quality T1 - freesurfer fails
EXP_22050_SESS01 - poor quality T1 - freesurfer fails
EXP_22361_SESS01 - poor quality T1 - freesurfer fails
EXP_21434_SESS01 - poor quality T1 - freesurfer fails
EXP_22590_SESS01 - poor quality T1 - freesurfer fails

Freesurfer questionable
EXP_22300_SESS01 - questionalable - some temporal lobe missing
EXP_21605_SESS01 - questionalable - some temporal lobe missing
EXP_21849_SESS01 - questionalable - some temporal lobe missing

[cifti_vis_fmri] ERROR: Expected fmri file /mnt/tigrlab/projects/dmiranda/NEW_HCP/DTI3T/hcp_lgi/DTI_CMH_S156_01/MNINonLinear/Results/REST/REST_Atlas_s0.dtseries.nii not found.
[cifti_vis_fmri] ERROR: Expected fmri file /mnt/tigrlab/projects/dmiranda/NEW_HCP/DTI3T/hcp_lgi/DTI_CMH_S171_01/MNINonLinear/Results/REST/REST_Atlas_s0.dtseries.nii not found.
[cifti_vis_fmri] ERROR: Expected fmri file /mnt/tigrlab/projects/dmiranda/NEW_HCP/PNSC/hcp_lgi/PNS_CMH_0023_01/MNINonLinear/Results/REST/REST_Atlas_s0.dtseries.nii not found.

# 2018-02-27

## Running PINT with other sampling radius options

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

BASE_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr10-6-12_CMH_20170227

basesubs=`ls -1d ${BASE_DATA}/*/hcp/*_CMH_*`
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
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_CMH -
```

## now the ZHH data

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/scratch/edickie/saba_PINT/data/ZHH/hcp
fmri_basedir=/KIMEL/tigrlab/external/miklos/epitome/EXP
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d EXP_2????_SESS01`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/data/ZHH/PINT_pcorr8-6-12_CMH_20170227

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 8 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_ZHH -
```

## copy all the outputs into two folders

```sh
ssh dev01

/imaging/scratch/kimel/edickie/saba_PINT
/KIMEL/tigrlab/scratch/edickie/saba_PINT/PINT_outputs_s8_8-6-12
/KIMEL/tigrlab/scratch/edickie/saba_PINT/PINT_outputs_s8_10-6-12

