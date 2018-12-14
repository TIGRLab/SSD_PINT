# commands to run the postPINT bits

## some postPINT concat bits..

run 2018-11-17 on kandel

```sh
outdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
export ciftify_container

singularity shell \
-H ${sing_home} \
-B ${outdir}:/output \
${ciftify_container}
```

#than inside the container

```sh
ciftify_postPINT1_concat /output/postPINT1_concat_all_qa_passes.csv `cat /output/qa_passes_summary_filelist.csv`
```

## then we can send command to run the sub2sub calculation

```sh
outdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
export ciftify_container

echo singularity exec \
-H ${sing_home} \
-B ${outdir}:/output \
${ciftify_container} \
ciftify_postPINT2_sub2sub \
/output/postPINT1_concat_all_qa_passes.csv \
/output/postPINT2_sub2sub_all_qa_passes.csv | \
qsub -V -l walltime=4:00:00,nodes=1:ppn=6 -N sub2sub -j oe -o ${outdir}
```
