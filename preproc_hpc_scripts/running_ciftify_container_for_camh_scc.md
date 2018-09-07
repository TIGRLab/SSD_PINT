# running the ZHH data on the SCC


```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/scratch/dmiranda/BIDS_ZHH
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
mkdir -p $sing_home


SIDlist=`cd ${bids_dir}; ls -1d sub* | sed 's/sub-//g'`
cd /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep
mkdir -p ZHH/out ZHH/work ZHH/logs
cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/input \
    -B ${outputdir}:/output \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img \
    /input /output/ZHH/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /output/ZHH/work \
    --fs-license /license_file.txt \
    --n_cpus 4 \
    --fmriprep-args="--use-aroma"  | \
    qsub -l walltime=23:00:00,nodes=1:ppn=4 -N ciftify_$SID -j oe -o ${outputdir}/ZHH/logs;
done
```

## but..23hrs may not have been enough for many subjects..

so we need a script to delete half completed outputs and rerun

There are three types of half completed outputs..

1. freesurfer is still running
    + in which case we need to delete `*Running` file as well as ciftify outputs
2. ciftify_recon_all didn't finish
    + we need to delete all ciftify outputs for that participant
    + if this is due to unifinished freesurfer we should remove the `*IsRunning`
3. ciftify_subject_fmri didn't finish
    + we need to delete the fmri/Results for this person

Then we can just rerun all participants to insure that all qc images are generated..

## 1. remove downstream from unifinished recon_all..

Note: if recon_all did finish..it is on the output folder..if not it is in the workdir..

```sh
bidsraw=/KIMEL/tigrlab/scratch/dmiranda/BIDS_ZHH
outdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/
cd ${bidsraw} ; ls -1d sub* | sort > ${outdir}/subsIn.txt

outdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/
cd ${outdir}/freesurfer ; ls -1d sub* | sort > ${outdir}/doneReconAll.txt

cd $outdir; diff subsIn.txt doneReconAll.txt
```

This prints the list of subjects who failed recon_all for some reason..

I got:

```
106,108d105
< sub-7740
< sub-7749
< sub-7793
159,161d155
< sub-8799
< sub-8831
< sub-8838
```

Now.. let us look through everyone who may not have finished ciftify_recon_all

```sh
for subject in `cat ${outdir}/doneReconAll.txt`; do
 cralog_lastline=`tail -3 ${outdir}/ciftify/$subject/cifti_recon_all.log | head -1`
 echo $subject,"$cralog_lastline"
done | grep -v Done > ${outdir}/ciftify_recon_all_errors.csv
```
Note: this file I opened in libreoffice and looked over. I decided that most of these files I should just delete the whole ciftify folder for:

```sh
todelete="sub-10391
sub-10563
sub-108
sub-6891
sub-7252
sub-8695
sub-8780
sub-9284
sub-9566
sub-10505
sub-10807
sub-11382
sub-6999
sub-7576
sub-7735
sub-7870
sub-8500
sub-9045
sub-9066
sub-9104
sub-9141
sub-9332
sub-9373
sub-9404
sub-9405
sub-9606
sub-9714
sub-9718
sub-10311
sub-11316
sub-9186
sub-9612
sub-11277
sub-10245
sub-9333
sub-8887
sub-8403
sub-8017
sub-9275
sub-10870
sub-11410
sub-6904
sub-6942
sub-7051
sub-7807
sub-8349
sub-8395
sub-84
sub-8418
sub-8484
sub-8557
sub-8974
sub-8999
sub-9027
sub-9229
sub-9436
sub-9591
sub-9724
sub-9799
sub-7991
sub-9577"

## note: I always run this first to check that everything looks good!
for sub in $todelete; do
  echo rm -r ${outdir}/ciftify/${sub}
done

## then this deletes
for sub in $todelete; do
  rm -r ${outdir}/ciftify/${sub}
done
```



```sh
for subject in `cat ${outdir}/doneReconAll.txt`; do
 cralog_lastline=`tail -3 ${outdir}/ciftify/$subject/cifti_recon_all.log | head -1`
 echo $subject $cralog_lastline
done | grep Done | sort | cut -d " " -f 1 > ${outdir}/done_ciftify_anat.txt


for subject in `cat ${outdir}/done_ciftify_anat.txt`; do
  ResultsLogs=`cd ${outdir}/ciftify/${subject}/MNINonLinear/Results/*/*.log`
  for ResultLog in ${ResultLogs}; do
  cralog_lastline=`tail -3 ${ResultLog} | head -1`
  echo ${subject} ${ResultLog} ${cralog_lastline}
done
done
```

```sh
fmrilogs=`ls ${outdir}/ciftify/sub-*/MNINonLinear/Results/*/ciftify_subject_fmri.log`
for log in $fmrilogs; do
  cralog_lastline=`tail -3 ${log} | head -1`
  echo ${log} ${cralog_lastline}
done | grep -v Done
```

# The above line gave me a small number to delete:

```
/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-10245/MNINonLinear/Results/ses-03_task-rest_bold/ciftify_subject_fmri.log Path: /output/ZHH/out/ciftify/sub-10245/MNINonLinear/fsaverage_LR32k
/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-11277/MNINonLinear/Results/ses-01_task-rest_bold/ciftify_subject_fmri.log Path: /output/ZHH/out/ciftify/sub-11277/MNINonLinear/fsaverage_LR32k
/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-7852/MNINonLinear/Results/ses-01_task-rest_bold/ciftify_subject_fmri.log -------------------------------------------------------------
/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-9725/MNINonLinear/Results/ses-01_task-rest_bold/ciftify_subject_fmri.log Running: wb_command -volume-to-surface-mapping /output/ZHH/out/ciftify/sub-9725/MNINonLinear/Results/ses-01_task-rest_bold/ses-01_task-rest_bold.nii.gz /output/ZHH/out/ciftify/sub-9725/MNINonLinear/Native/sub-9725.R.midthickness.native.surf.gii /tmp/tmpyyx2rat2/MNINonLinear/native/sub-9725.R.ses-01_task-rest_bold.native.func.gii -ribbon-constrained /output/ZHH/out/ciftify/sub-9725/MNINonLinear/Native/sub-9725.R.white.native.surf.gii /output/ZHH/out/ciftify/sub-9725/MNINonLinear/Native/sub-9725.R.pial.native.surf.gii -volume-roi /tmp/tmpyyx2rat2/goodvoxels.nii.gz
[edickie@scclogin02 out]$ rm -r /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-10245/MNINonLinear/Results/ses-03_task-rest_bold
[edickie@scclogin02 out]$ rm -r /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-11277/MNINonLinear/Results/ses-01_task-rest_bold
[edickie@scclogin02 out]$ rm -r /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-7852/MNINonLinear/Results/ses-01_task-rest_bold
[edickie@scclogin02 out]$ rm -r /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out//ciftify/sub-9725/MNINonLinear/Results/ses-01_task-rest_bold
```

## Now that certain files have been deleted..we are set to rerun

```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/scratch/dmiranda/BIDS_ZHH
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
mkdir -p $sing_home


SIDlist=`cd ${bids_dir}; ls -1d sub* | sed 's/sub-//g'`
cd /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep
mkdir -p ZHH/out ZHH/work ZHH/logs
cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/input \
    -B ${outputdir}:/output \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img \
    /input /output/ZHH/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /output/ZHH/work \
    --fs-license /license_file.txt \
    --n_cpus 4 \
    --fmriprep-args="--use-aroma"  | \
    qsub -l walltime=23:00:00,nodes=1:ppn=4 -N ciftify2_$SID -j oe -o ${outputdir}/ZHH/logs;
done
```

## now let's try running the cleaning and PINT bit..

```sh
ssh dev02
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/bids-on-scinet/examples/participant_ciftify_clean_and_PINT.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

for preprocfile in `ls ${outputdir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_task-rest_bold_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
echo ${cleaning_script} ${subject} ${session} task-rest_bold ${outputdir} ${sing_home} ${ciftify_container} | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N pint_${subject}_${session} -j oe -o ${outputdir}/../../ZHH/logs;
done
```

## add it turns out that we need to run the COBRE stuff on the SCC too...

```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
mkdir -p $sing_home


SIDlist=`cd ${bids_dir}/COBRE/COBRE; ls -1d sub* | sed 's/sub-//g'`
cd /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep
mkdir -p COBRE/out COBRE/work COBRE/logs
cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img \
    /bids/COBRE/COBRE /bids/COBRE/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /bids/COBRE/work \
    --fs-license /license_file.txt \
    --n_cpus 4 \
    --fmriprep-args="--use-aroma"  | \
    qsub -V -l walltime=23:00:00,nodes=1:ppn=4 -N ciftify_$SID -j oe -o ${bids_dir}/COBRE/logs;
done
```

## and now running subcortical extraction

```sh
ssh dev02
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_ciftify_clean_and_subcortical.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

for preprocfile in `ls ${outputdir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_task-rest_bold_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
echo ${cleaning_script} ${subject} ${session} task-rest_bold ${outputdir} ${sing_home} ${ciftify_container} | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N subts_${subject}_${session} -j oe -o ${outputdir}/../../ZHH/logs;
done
```

## running ghetto (no container) PINT qc for ZHH
```sh
ssh kandel
module load /projects/edickie/privatemodules/ciftify/201803
outputdir=/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out
func_base="task-rest_bold"
for preprocfile in `ls ${outputdir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_task-rest_bold_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
 cifti_vis_PINT subject \
      --ciftify-work-dir ${outputdir}/ciftify/ \
      --qcdir ${outputdir}/ciftify_PINT/qc \
      ${outputdir}/ciftify_clean_img/${subject}/${session}/${subject}_${session}_${func_base}_desc-clean_bold.dtseries.nii \
      ${subject} \
      ${outputdir}/ciftify_PINT/${subject}/${session}/${subject}_${session}_${func_base}_desc-clean_bold_summary.csv
done
```

## so using my newer metrics for determining if stuff finished - there are some more ZHH peeps to rerun

people to rerun only for ciftify_fmri

```sh
rm /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/freesurfer/sub*/scripts/IsRunning*
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-9061
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8873
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-11316
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-9612
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-7927
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8350
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-7960
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-9390
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8619
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8586
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-6891
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-9141
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-9125
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8946
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-7937
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8889
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-9194
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8949
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-11242/MNINonLinear/Results/ses-02_task-rest_bold
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-9718/MNINonLinear/Results/ses-01_task-rest_bold
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-8603/MNINonLinear/Results/ses-01_task-rest_bold
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify/sub-6729/MNINonLinear/Results/ses-01_task-rest_bold
```

```sh
ssh dev01
bids_dir=/KIMEL/tigrlab/scratch/dmiranda/BIDS_ZHH
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home

SIDlist="7927
7937
7960
8619
8889
8946
8949
9061
9125
9390
9612
11316
9194
6891
8350
8586
8873
9141
6729
8603
9718
11242"

mkdir -p $sing_home
cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/input \
    -B ${outputdir}:/output \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img \
    /input /output/ZHH/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /output/ZHH/work \
    --fs-license /license_file.txt \
    --n_cpus 4 \
    --fmriprep-args="--use-aroma"  | \
    qsub -l walltime=23:00:00,nodes=1:ppn=4 -N ciftify3_$SID -j oe -o ${outputdir}/ZHH/logs;
done

```

people to rerun

```sh
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/out/ciftify/sub-A00028052
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/out/ciftify/sub-A00000456
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/out/ciftify/sub-A00004087
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/out/ciftify/sub-A00000909/MNINonLinear/Results/ses-20110101_task-rest_run-01_bold
rm -r /scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/out/ciftify/sub-A00024301/MNINonLinear/Results/ses-20090101_task-rest_bold
```

```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home

SIDlist="A00000456
A00004087
A00028052
A00000909
A00024301"

cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img \
    /bids/COBRE/COBRE /bids/COBRE/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /bids/COBRE/work \
    --fs-license /license_file.txt \
    --n_cpus 4 \
    --fmriprep-args="--use-aroma"  | \
    qsub -V -l walltime=23:00:00,nodes=1:ppn=4 -N ciftify2_$SID -j oe -o ${bids_dir}/COBRE/logs;
done
```