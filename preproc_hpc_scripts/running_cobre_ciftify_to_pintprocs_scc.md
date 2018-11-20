# running cobre (all the way..)

So it turns out we were missing a lot of cobre participants in the first run
So we will rerun then all. We will see if upping the number of CPU's from 4 to 6 on the scc decreases the walltime


```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/external/SchizConnect/COBRE/bids/
bids_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
mkdir -p $sing_home


SIDlist=`cd ${bids_dir}/COBRE/; ls -1d sub* | sed 's/sub-//g'`
cd ${bids_out}
mkdir -p COBRE/out COBRE/work COBRE/logs

cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B ${bids_out}:/out \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    ${ciftify_container} \
    /bids/COBRE/ /out/COBRE/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /out/COBRE/work \
    --fs-license /license_file.txt \
    --n_cpus 6  | \
    qsub -V -l walltime=23:00:00,nodes=1:ppn=6 -N ciftify_$SID -j oe -o ${bids_dir}/COBRE/logs;
done
```

## 2018-11-11 Running the cobre mriqc on the scc too

```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/external/SchizConnect/COBRE/bids/
bids_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
mkdir -p $sing_home


SIDlist=`cd ${bids_dir}/COBRE/; ls -1d sub* | sed 's/sub-//g'`
cd /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep
mkdir -p COBRE/out COBRE/work COBRE/logs
cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B ${bids_out}:/out \
    /KIMEL/tigrlab/archive/code/containers/MRIQC/poldracklab_mriqc_0.11.0-2018-06-05-1e4ac9792325.img \
    /bids/COBRE/ /out/COBRE/out participant \
    --participant_label=$SID \
    -w /out/COBRE/work \
    --n_procs 2 \
    --no-sub | \
    qsub -V -l walltime=4:00:00,nodes=1:ppn=2 -N mriqc_$SID -j oe -o ${bids_out}/COBRE/logs;
done
```

## Note these are the things you need to set up in the output folder before this scripts will work

1. cleaning configs need to be in the ciftify_clean_img output dir
2. PINT template csvs and the color look up table should be in the ciftify_PINT output folder
3. the subcortical ROI templates need to be in the output folders

```sh
ssh dev02
dataset="COBRE"
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/logs
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_ciftify_clean_PINT_subcortts_dlabels.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

mkdir -p ${outputdir}/ciftify_meants/templates ${outputdir}/ciftify_PINT ${outputdir}/ciftify_clean_img
cp ${outputdir}/../../ZHH/out/ciftify_clean_img/*.json ${outputdir}/ciftify_clean_img/
cp ${outputdir}/../../ZHH/out/ciftify_PINT/Yeo* ${outputdir}/ciftify_PINT/
cp ${outputdir}/../../ZHH/out/ciftify_meants/templates/*..nii ${outputdir}/ciftify_meants/templates/

for preprocfile in `ls ${outputdir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_task-rest_run-01_bold_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_bold_desc-cleansm0_atlas-7RSN_roi-Rthalamus_timeseries.csv ]; then
echo ${cleaning_script} ${subject} ${session} task-rest_run-01_bold ${outputdir} ${outputdir} ${sing_home} ${ciftify_container}  | qsub -V -l walltime=00:40:00,nodes=1:ppn=4 -N subts_${subject}_${session}_run-01 -j oe -o ${logsdir};
fi
done

for preprocfile in `ls ${outputdir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_task-rest_run-02_bold_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_bold_desc-cleansm0_atlas-7RSN_roi-Rthalamus_timeseries.csv ]; then
echo ${cleaning_script} ${subject} ${session} task-rest_run-02_bold ${outputdir} ${outputdir} ${sing_home} ${ciftify_container}  | qsub -V -l walltime=00:40:00,nodes=1:ppn=4 -N subts_${subject}_${session}_run-02 -j oe -o ${logsdir};
fi
done

for preprocfile in `ls ${outputdir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_task-rest_bold_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_bold_desc-cleansm0_atlas-7RSN_roi-Rthalamus_timeseries.csv ]; then
echo ${cleaning_script} ${subject} ${session} task-rest_bold ${outputdir} ${outputdir} ${sing_home} ${ciftify_container}  | qsub -V -l walltime=00:40:00,nodes=1:ppn=4 -N subts_${subject}_${session} -j oe -o ${logsdir};
fi
done
```


----------

# The ugly - stuff that needed to be rerun..

# 2018-11-03

Dunno why 5 subjects have failed so far.. and there's a typo in the logs so I'm not getting good info.
Will delete and rerun so that they are fixed

The good news (from the logs is that most are finishing the full pipeline in under 23 hrs..)
So it looks like increaseing the CPU's to 6 was a good idea..

215545.headnode            ...ify_A00000300 edickie         16:10:35 C medium
215546.headnode            ...ify_A00000368 edickie         25:04:39 C medium
215547.headnode            ...ify_A00000456 edickie         01:02:11 C medium
215551.headnode            ...ify_A00001181 edickie         01:04:24 C medium
215552.headnode            ...ify_A00001243 edickie         15:11:10 C medium
215553.headnode            ...ify_A00001251 edickie         17:34:18 C medium
215556.headnode            ...ify_A00002405 edickie         26:29:44 C medium
215557.headnode            ...ify_A00002480 edickie         14:08:46 C medium
215558.headnode            ...ify_A00003150 edickie         16:03:58 C medium
215560.headnode            ...ify_A00004507 edickie         23:22:40 C medium
215564.headnode            ...ify_A00009656 edickie         00:00:40 C medium
215565.headnode            ...ify_A00009946 edickie         24:43:25 C medium
215566.headnode            ...ify_A00010150 edickie         19:03:04 C medium
215567.headnode            ...ify_A00010684 edickie         18:14:47 C medium
215568.headnode            ...ify_A00011107 edickie         00:00:18 C medium
215569.headnode            ...ify_A00011265 edickie         17:40:43 C medium
215572.headnode            ...ify_A00012995 edickie         18:41:29 C medium
215573.headnode            ...ify_A00013140 edickie         17:33:39 C medium
215583.headnode            ...ify_A00014636 edickie         00:32:26 C medium

A00000456
A00001181
A00009656
A00011107
A00014636
```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/external/SchizConnect/COBRE/bids/
bids_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
mkdir -p $sing_home

SIDlist="A00000456
A00001181
A00009656
A00011107
A00014636"

cd ${bids_out}
for SID in  $SIDlist; do
  ls COBRE/out/freesurfer/sub-${SID}
  ls COBRE/out/fmriprep/sub-${SID}/
  ls COBRE/out/fmriprep/sub-${SID}.html
  rm -r COBRE/out/freesurfer/sub-${SID}
  rm -r COBRE/out/fmriprep/sub-${SID}/
  rm COBRE/out/fmriprep/sub-${SID}.html
done
mkdir -p COBRE/work2

cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B ${bids_out}:/out \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    ${ciftify_container} \
    /bids/COBRE/ /out/COBRE/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /out/COBRE/work2 \
    --fs-license /license_file.txt \
    --n_cpus 6  | \
    qsub -V -l walltime=23:00:00,nodes=1:ppn=6 -N ciftify_$SID -j oe -o ${bids_out}/COBRE/logs;
done
```

output from qstat | grep C on 2019-11-04

215548.headnode            ...ify_A00000541 edickie         21:17:39 C medium
215549.headnode            ...ify_A00000838 edickie         26:18:05 C medium
215550.headnode            ...ify_A00000909 edickie         25:10:11 C medium
215552.headnode            ...ify_A00001243 edickie         15:11:10 C medium
215553.headnode            ...ify_A00001251 edickie         17:34:18 C medium
215554.headnode            ...ify_A00001452 edickie         23:47:36 C medium
215555.headnode            ...ify_A00002198 edickie         23:20:38 C medium
215556.headnode            ...ify_A00002405 edickie         26:29:44 C medium
215558.headnode            ...ify_A00003150 edickie         16:03:58 C medium
215559.headnode            ...ify_A00004087 edickie         23:28:01 C medium
215560.headnode            ...ify_A00004507 edickie         23:22:40 C medium
215561.headnode            ...ify_A00006754 edickie         22:16:39 C medium
215562.headnode            ...ify_A00007409 edickie         23:51:36 C medium
215563.headnode            ...ify_A00009280 edickie         26:18:45 C medium
215565.headnode            ...ify_A00009946 edickie         24:43:25 C medium
215566.headnode            ...ify_A00010150 edickie         19:03:04 C medium
215567.headnode            ...ify_A00010684 edickie         18:14:47 C medium
215569.headnode            ...ify_A00011265 edickie         17:40:43 C medium
215570.headnode            ...ify_A00011725 edickie         24:34:39 C medium
215571.headnode            ...ify_A00012767 edickie         24:35:40 C medium
215572.headnode            ...ify_A00012995 edickie         18:41:29 C medium
215573.headnode            ...ify_A00013140 edickie         17:33:39 C medium
215574.headnode            ...ify_A00013216 edickie         25:01:43 C medium
215575.headnode            ...ify_A00013363 edickie         23:31:23 C medium
215576.headnode            ...ify_A00013816 edickie         24:52:27 C medium
215577.headnode            ...ify_A00014120 edickie         23:41:44 C medium
215578.headnode            ...ify_A00014175 edickie         20:10:08 C medium
215579.headnode            ...ify_A00014225 edickie         21:59:09 C medium
215580.headnode            ...ify_A00014522 edickie         16:58:31 C medium
215581.headnode            ...ify_A00014590 edickie         17:05:59 C medium
215582.headnode            ...ify_A00014607 edickie         18:52:29 C medium
215584.headnode            ...ify_A00014719 edickie         17:29:50 C medium
215585.headnode            ...ify_A00014804 edickie         24:11:43 C medium
215586.headnode            ...ify_A00014830 edickie         23:26:36 C medium
215587.headnode            ...ify_A00014839 edickie         15:57:02 C medium
215588.headnode            ...ify_A00014898 edickie         17:59:04 C medium
215589.headnode            ...ify_A00015201 edickie         19:25:06 C medium
215590.headnode            ...ify_A00015518 edickie         22:05:45 C medium
215591.headnode            ...ify_A00015648 edickie         17:30:18 C medium
215592.headnode            ...ify_A00015759 edickie         17:36:01 C medium
215593.headnode            ...ify_A00015826 edickie         18:04:53 C medium
215594.headnode            ...ify_A00016197 edickie         24:31:51 C medium
215595.headnode            ...ify_A00016720 edickie         01:45:45 C medium
215596.headnode            ...ify_A00016723 edickie         24:47:54 C medium
215597.headnode            ...ify_A00017147 edickie         23:18:23 C medium
215598.headnode            ...ify_A00017294 edickie         21:05:16 C medium
215599.headnode            ...ify_A00018129 edickie         17:55:01 C medium
215600.headnode            ...ify_A00018317 edickie         22:41:47 C medium
215601.headnode            ...ify_A00018335 edickie         23:42:45 C medium
215602.headnode            ...ify_A00018403 edickie         29:44:26 C medium
215603.headnode            ...ify_A00018434 edickie         22:58:22 C medium
215604.headnode            ...ify_A00018553 edickie         21:24:23 C medium
215605.headnode            ...ify_A00018598 edickie         00:00:37 C medium
215606.headnode            ...ify_A00018716 edickie         25:00:27 C medium
215607.headnode            ...ify_A00018979 edickie         01:07:17 C medium
215608.headnode            ...ify_A00019293 edickie         24:42:35 C medium
215609.headnode            ...ify_A00019349 edickie         24:45:35 C medium
215610.headnode            ...ify_A00019750 edickie         23:34:13 C medium
215611.headnode            ...ify_A00019888 edickie         16:55:33 C medium
215612.headnode            ...ify_A00020414 edickie         17:24:33 C medium
215613.headnode            ...ify_A00020416 edickie         17:18:55 C medium
215614.headnode            ...ify_A00020602 edickie         25:44:49 C medium
215615.headnode            ...ify_A00020787 edickie         24:19:33 C medium
215616.headnode            ...ify_A00020805 edickie         25:57:00 C medium
215617.headnode            ...ify_A00020895 edickie         01:48:36 C medium
215618.headnode            ...ify_A00020968 edickie         23:46:53 C medium
215619.headnode            ...ify_A00020984 edickie         25:58:10 C medium
215620.headnode            ...ify_A00021058 edickie         25:05:18 C medium
215621.headnode            ...ify_A00021072 edickie         19:59:16 C medium
215622.headnode            ...ify_A00021081 edickie         00:04:20 C medium
215624.headnode            ...ify_A00021145 edickie         20:17:26 C medium
215626.headnode            ...ify_A00021598 edickie         00:05:21 C medium
215645.headnode            ...ify_A00023132 edickie         01:48:49 C medium
215652.headnode            ...ify_A00023366 edickie         01:39:52 C medium
[edickie@dev01 ZHH]$

so these one's are the list to try to restart..
A00016720
A00018598
A00018979
A00020895
A00021081
A00021598
A00023132
A00023366

```sh
ssh dev01
bids_dir=/KIMEL/tigrlab/external/SchizConnect/COBRE/bids/
bids_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
mkdir -p $sing_home

SIDlist="A00016720
A00018598
A00018979
A00020895
A00021081
A00021598
A00023132
A00023366"

cd ${bids_out}
for SID in  $SIDlist; do
  ls COBRE/out/freesurfer/sub-${SID}
  ls COBRE/out/fmriprep/sub-${SID}/
  ls COBRE/out/fmriprep/sub-${SID}.html
  rm -r COBRE/out/freesurfer/sub-${SID}
  rm -r COBRE/out/fmriprep/sub-${SID}/
  rm COBRE/out/fmriprep/sub-${SID}.html
done
mkdir -p COBRE/work2

cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B ${bids_out}:/out \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    ${ciftify_container} \
    /bids/COBRE/ /out/COBRE/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /out/COBRE/work2 \
    --fs-license /license_file.txt \
    --n_cpus 6  | \
    qsub -V -l walltime=23:00:00,nodes=1:ppn=6 -N ciftify_$SID -j oe -o ${bids_out}/COBRE/logs;
done
```

no functional run?
sub-A00023132 - no functional data
sub-A00023366- no functional data
no fmriprep at all
sub-A00024510 - no T1w data
sub-A00027119 - no T1w data        

A00009656 - no T1w image found      
A00011107 - no T1w image found                      
A00018598 - no T1w image found                


## running this last bit to get QA index pages to write

```sh
ssh dev01
bids_dir=/KIMEL/tigrlab/external/SchizConnect/COBRE/bids/
bids_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
mkdir -p $sing_home


singularity run -H ${sing_home}:/myhome \
  -B ${bids_dir}:/bids \
  -B ${bids_out}:/out \
  -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
  ${ciftify_container} \
  /bids/COBRE/ /out/COBRE/out group \
  --fs-license /license_file.txt

```

A00016720 - fmriprep done
sub-A00020895
sub-A00021081
sub-A00021598
sub-A00024684
sub-A00031186
sub-A00037007

fmriprep not done?
sub-A00018598
sub-A00022592
sub-A00024198
sub-A00025969
sub-A00027755
sub-A00027969
sub-A00031597
sub-A00033812
A00036555
A00036897
A00037619
sub-A00026945 - so many anats..maybe some need to be blacklisted?


```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/external/SchizConnect/COBRE/bids/
bids_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
mkdir -p $sing_home

## rerun the ones that just failed ciftify (but are done fMRIprep)
SIDlist="A00016720
A00020895
A00021081
A00021598
A00024684
A00031186
A00037007"

cd ${bids_out}
for SID in  $SIDlist; do
  rm -r COBRE/out/ciftify/sub-${SID}
done
mkdir -p COBRE/work3

cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B ${bids_out}:/out \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    ${ciftify_container} \
    /bids/COBRE/ /out/COBRE/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /out/COBRE/work2 \
    --fs-license /license_file.txt \
    --n_cpus 1  | \
    qsub -V -l walltime=6:00:00,nodes=1:ppn=1 -N ciftify_$SID -j oe -o ${bids_out}/COBRE/logs;
done

## rerun these ones from the start
SIDlist="A00018598
A00022592
A00024198
A00025969
A00027755
A00027969
A00031597
A00033812
A00036555
A00036897
A00037619
A00026945"

cd ${bids_out}
for SID in  $SIDlist; do
  ls COBRE/out/freesurfer/sub-${SID}
  ls COBRE/out/fmriprep/sub-${SID}/
  ls COBRE/out/fmriprep/sub-${SID}.html
  rm -r COBRE/out/freesurfer/sub-${SID}
  rm -r COBRE/out/fmriprep/sub-${SID}/
  rm COBRE/out/fmriprep/sub-${SID}.html
  rm -r COBRE/out/ciftify/sub-${SID}
done

cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B ${bids_out}:/out \
    -B /quarantine/Freesurfer/6.0.0/freesurfer/license.txt:/license_file.txt \
    ${ciftify_container} \
    /bids/COBRE/ /out/COBRE/out participant \
    --participant_label=$SID \
    --fmriprep-workdir /out/COBRE/work3 \
    --fs-license /license_file.txt \
    --n_cpus 6  | \
    qsub -V -l walltime=23:00:00,nodes=1:ppn=6 -N ciftify_$SID -j oe -o ${bids_out}/COBRE/logs;
done
```

# 2018-11-11 Three participants seem to have timed out for mriqc..will try to restart them..

edickie@dev01 logs]$ grep Resources mriqc_* | grep e=04
mriqc_A00000909.o220742: Resources:           cput=07:38:01,mem=3363512kb,vmem=9039244kb,walltime=04:00:17
mriqc_A00009280.o220755: Resources:           cput=07:50:25,mem=3172820kb,vmem=9267484kb,walltime=04:00:13
mriqc_A00018979.o220799: Resources:           cput=07:50:14,mem=2990788kb,vmem=10049652kb,walltime=04:00:09

```sh
ssh dev02
bids_dir=/KIMEL/tigrlab/external/SchizConnect/COBRE/bids/
bids_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
mkdir -p $sing_home


SIDlist="A00000909
A00009280
A00018979"

cd /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep
mkdir -p COBRE/out COBRE/work COBRE/logs
cd $sing_home
for SID in $SIDlist; do
  subject=$SID
  echo singularity run -H ${sing_home}:/myhome \
    -B ${bids_dir}:/bids \
    -B ${bids_out}:/out \
    /KIMEL/tigrlab/archive/code/containers/MRIQC/poldracklab_mriqc_0.11.0-2018-06-05-1e4ac9792325.img \
    /bids/COBRE/ /out/COBRE/out participant \
    --participant_label=$SID \
    -w /out/COBRE/work \
    --n_procs 2 \
    --no-sub | \
    qsub -V -l walltime=4:00:00,nodes=1:ppn=2 -N mriqc2_$SID -j oe -o ${bids_out}/COBRE/logs;
done
```
