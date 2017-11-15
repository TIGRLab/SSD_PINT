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

## The pint step

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
## so it turns out that Dayton ran the corr version instead of the pcorr version of PINT

Therefore..we are gonna rerun it for all networks..

# 2017-11-06 running COBRE pint

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/projects/dmiranda/NEW_HCP/COBRE/HCP_COBRE
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d COBRE_*_SESS01`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/data/ZHH/PINT_pcorr6-6-12

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_COBRE -
```



## 2017-11-07

Note: this running everything into a sub-directory of ZHH was stupid..

```sh
cd /imaging/scratch/kimel/edickie/saba_PINT/data/ZHH/
mv PINT_pcorr6-6-12 ../..
```

### first ASDD

```sh
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/projects/dmiranda/NEW_HCP/ASDD/HCP_ASDD
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d ASDD_CMH_*_01`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_ASDD -
```

### now DTI3T

```sh
ssh dev02
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/projects/dmiranda/NEW_HCP/DTI3T/hcp_lgi
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d DTI_CMH_*`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_DTI -
```
### now PNS

```sh
ssh dev01
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/projects/dmiranda/NEW_HCP/PNSC/hcp_lgi
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d PNS_CMH_*`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_PNS -
```

### now RTMSWM

```sh
ssh dev01
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/projects/dmiranda/NEW_HCP/RTMSWM/hcp_lgi
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d RTMSWM_CMH_*`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_RTMS -
```

### now SPINS

```sh
ssh dev01
module load /KIMEL/quarantine/modules/quarantine
module load Freesurfer/6.0.0
module load FSL/5.0.9-ewd
module load connectome-workbench/1.2.3
module load python/3.6_ciftify_01
module load GNU_PARALLEL/20170122

HCP_DATA=/KIMEL/tigrlab/projects/dmiranda/NEW_HCP/SPINS/hcp_lgi
export CIFTIFY_TEMPLATES=${HOME}/code/ciftify/ciftify/data
export TMPDIR=/export/ramdisk

subject_list=`cd ${HCP_DATA}; ls -1d SPN01_CMH_*`
pint_outputdir=/imaging/scratch/kimel/edickie/saba_PINT/PINT_pcorr6-6-12

mkdir -p  ${pint_outputdir}
cd ${pint_outputdir}

parallel "echo ciftify_PINT_vertices \
  --pcorr --sampling-radius 6 --search-radius 6 --padding-radius 12 \
  ${HCP_DATA}/{}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.L.midthickness.32k_fs_LR.surf.gii \
  ${HCP_DATA}/{}/MNINonLinear/fsaverage_LR32k/{}.R.midthickness.32k_fs_LR.surf.gii \
  ${CIFTIFY_TEMPLATES}/PINT/Yeo7_2011_80verts.csv \
  ${pint_outputdir}/{}/{}" ::: ${subject_list} |
  qbatch --walltime 1:00:00 -c 1 -j 1 --ppj 6 -N  pint_SPINS -
```

## running all the qc...

We did this local on nissl

for ZHH

```sh
module load /projects/edickie/privatemodules/ciftify/201707
export HCP_DATA=/scratch/edickie/saba_PINT/data/ZHH/hcp

cd ${HCP_DATA}

subject_list=`ls -1d EXP_2*`

for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_recon_all snaps ${subject}
done
cifti_vis_recon_all index
```
## running COBRE

```sh
module load /projects/edickie/privatemodules/ciftify/201707
export HCP_DATA=/projects/dmiranda/NEW_HCP/COBRE/HCP_COBRE

subject_list=`cd ${HCP_DATA}; ls -1d COBRE_*_SESS01`

cd /scratch/edickie/saba_PINT/qc

for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_recon_all snaps --qcdir ${PWD}/qc_COBRE_recon_all ${subject}
done
cifti_vis_recon_all index

export HCP_DATA=/projects/dmiranda/NEW_HCP/ASDD/HCP_ASDD
subject_list=`cd ${HCP_DATA}; ls -1d ASDD_CMH_*_01`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_recon_all snaps --qcdir ${PWD}/qc_ASDD_recon_all ${subject}
done
cifti_vis_recon_all index

export HCP_DATA=/projects/dmiranda/NEW_HCP/DTI3T/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d DTI_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_recon_all snaps --qcdir ${PWD}/qc_DTI3T_recon_all ${subject}
done
cifti_vis_recon_all index

export HCP_DATA=/projects/dmiranda/NEW_HCP/SPINS/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d SPN01_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_recon_all snaps --qcdir ${PWD}/qc_SPINS_recon_all ${subject}
done
cifti_vis_recon_all index

export HCP_DATA=/projects/dmiranda/NEW_HCP/PNSC/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d PNS_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_recon_all snaps --qcdir ${PWD}/qc_PNS_recon_all ${subject}
done

export HCP_DATA=/projects/dmiranda/NEW_HCP/RTMSWM/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d RTMSWM_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_recon_all snaps --qcdir ${PWD}/qc_RTMSWM_recon_all ${subject}
done

```
## fMRI qc

```sh
module load /projects/edickie/privatemodules/ciftify/201707
export HCP_DATA=/projects/dmiranda/NEW_HCP/COBRE/HCP_COBRE

subject_list=`cd ${HCP_DATA}; ls -1d COBRE_*_SESS01`

cd /scratch/edickie/saba_PINT/qc

for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_fmri snaps --qcdir ${PWD}/qc_COBRE_fmri --SmoothingFWHM 8 REST ${subject}
done

for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_fmri snaps --qcdir ${PWD}/qc_RTMSWM_fmri --SmoothingFWHM 8 REST ${subject}
done

export HCP_DATA=/projects/dmiranda/NEW_HCP/SPINS/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d SPN01_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_fmri snaps --qcdir ${PWD}/qc_SPINS_fmri --SmoothingFWHM 8 REST ${subject}
done

export HCP_DATA=/projects/dmiranda/NEW_HCP/PNSC/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d PNS_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_fmri snaps --qcdir ${PWD}/qc_PNS_fmri --SmoothingFWHM 8 REST ${subject}
done

export HCP_DATA=/projects/dmiranda/NEW_HCP/DTI3T/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d DTI_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_fmri snaps --qcdir ${PWD}/qc_DTI3T_fmri --SmoothingFWHM 8 REST ${subject}
done

export HCP_DATA=/projects/dmiranda/NEW_HCP/ASDD/HCP_ASDD
subject_list=`cd ${HCP_DATA}; ls -1d ASDD_CMH_*_01`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_fmri snaps --qcdir ${PWD}/qc_ASDD_fmri --SmoothingFWHM 8 REST ${subject}
done
```
## running PINTqc

```sh
pint_outputdir=/scratch/edickie/saba_PINT/PINT_pcorr6-6-12
export HCP_DATA=/projects/dmiranda/NEW_HCP/ASDD/HCP_ASDD
subject_list=`cd ${HCP_DATA}; ls -1d ASDD_CMH_*_01`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_PINT snaps --qcdir ${PWD}/qc_ASDD_PINT-pcorr6-6-12 \
  ${HCP_DATA}/${subject}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${subject} ${pint_outputdir}/${subject}/${subject}_summary.csv
done

pint_outputdir=/scratch/edickie/saba_PINT/PINT_pcorr6-6-12
export HCP_DATA=/projects/dmiranda/NEW_HCP/DTI3T/hcp_lgi
subject_list=`cd ${HCP_DATA}; ls -1d DTI_CMH_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_PINT snaps --qcdir ${PWD}/qc_DTI3T_PINT-pcorr6-6-12 \
  ${HCP_DATA}/${subject}/MNINonLinear/Results/REST/REST_Atlas_s8.dtseries.nii \
  ${subject} ${pint_outputdir}/${subject}/${subject}_summary.csv
done

pint_outputdir=/scratch/edickie/saba_PINT/PINT_pcorr6-6-12
export HCP_DATA=/scratch/edickie/saba_PINT/data/ZHH/hcp
subject_list=`cd ${HCP_DATA}; ls -1d EXP_*`
for subject in ${subject_list}; do
  echo ${subject}
  cifti_vis_PINT snaps --qcdir ${PWD}/qc_ZHH_PINT-pcorr6-6-12 \
  ${HCP_DATA}/${subject}/MNINonLinear/Results/REST_01/REST_01_Atlas_s8.dtseries.nii \
  ${subject} ${pint_outputdir}/${subject}/${subject}_summary.csv
done

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
