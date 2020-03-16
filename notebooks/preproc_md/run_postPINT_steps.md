# commands to run the postPINT bits

## some postPINT concat bits..

re-run 2020-02-24 on kandel - to get the `all_clinicalplusqa_group` results



```sh
outdir=/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/tmp/sing_home-lasdfjoasie
ciftify_container=/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
export ciftify_container

mkdir $sing_home 
singularity shell \
-H ${sing_home} \
-B ${outdir}:/output \
${ciftify_container}
```

### than inside the container

```sh
cd /output
mkdir all_clinicalplusqa_group/postPINT
ciftify_postPINT1_concat \
  /output/all_clinicalplusqa_group/postPINT/postPINT1_concat_all_qa_passes.csv \
  `cat /output/all_clinicalplusqa_group/pint_summary_filelist.csv`
```



### then we can send command to run the sub2sub calculation 

(this takes a while - was run overnight)


```sh
singularity exec \
-H ${sing_home} \
-B ${outdir}:/output \
${ciftify_container} \
ciftify_postPINT2_sub2sub \
/output/all_clinicalplusqa_group/postPINT/postPINT1_concat_all_qa_passes.csv \
/output/all_clinicalplusqa_group/postPINT/postPINT2_sub2sub_all_qa_passes.csv 

```

## checking seedcorr data sizes

We also wrote a little script to do this..so we will try to run that..
```sh
ssh kandel
outdir=/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/tmp/sing_home-lasdfjoasie
sz_pint_scripts=/projects/edickie/code/SZ_PINT
ciftify_container=/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
export ciftify_container

singularity exec \
-H ${sing_home} \
-B ${outdir}:/output \
-B ${sz_pint_scripts}:/scripts \
${ciftify_container} python /scripts/code/bin/seedcor_number_voxels_above_threshold.py  \
  /output/ \
  /scripts/data/processed/pheno/20200221_pheno_clinicalplusqa.csv \
  /output/all_clinicalplusqa_group/postPINT/seedcorr_numvxabove_all_qa_passes.csv
```

note: one participant failed this last step - `DTI3T/sub-CMHH170/ses-01` - if we do need the seedcorr info we will have to reload these things..
