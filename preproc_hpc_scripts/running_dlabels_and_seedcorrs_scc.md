## running the dlabel creation and seedcorr creation on the scc

These should add dlabel files and seed files into the PINT outputs

subject=$1
session=$2
func_base=$3
archive_pipedir=$4
outdir=$5
sing_home=$6
ciftify_container=$7

```sh
ssh dev01
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out
archive_pipedir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4


for preprocfile in `ls ${outputdir}/fmriprep/sub-*/func/sub-*_task-rest_bold_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session="none"
  if [ ! -f ${outputdir}/ciftify_PINT/${subject}/${subject}_task-rest_bold_atlas-pvertexNET_roi-7_fcmap.dscalar.nii ]; then
echo ${cleaning_script} ${subject} ${session} task-rest_bold ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container} # | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N dlabels_${subject}_${session} -j oe -o ${outputdir}/../../ds000030_R1.0.5/logs;
fi
done
```

KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh sub-10171 none task-rest_bold /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh sub-10227 none task-rest_bold /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh sub-10269 none task-rest_bold /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh sub-11122 none task-rest_bold /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh sub-50029 none task-rest_bold /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out /KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home /KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
[edickie@dev01 ~]$
While running:

# 2016-11-11 running all the CAMH data

```sh
ssh dev01

sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

for dataset in "ASDD" "PNSC" "DTI3T" "SPINS" "RTMSWM"; do

archive_pipedir=/KIMEL/tigrlab/scratch/jjeyachandra/test_env/archive/data/${dataset}/pipelines/
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/logs

cd ${sing_home}

cp ${outputdir}/../../ds000030_R1.0.5/out/ciftify_PINT/Yeo7_2011_80verts_roiidx_LUT.txt ${outputdir}/ciftify_PINT/

func_base="task-rest_acq-CMH_run-01_bold"
for preprocfile in `ls ${archive_pipedir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_${func_base}_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))

  if [ ! -f ${outputdir}/ciftify_PINT/${subject}/${session}/${subject}_${session}_${func_base}_atlas-pvertexNET_roi-7_fcmap.dscalar.nii ]; then
echo ${cleaning_script} ${subject} ${session} ${func_base} ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}   | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N dlabels_${subject}_${session} -j oe -o ${logsdir};

fi
done
done
```

# 2016-11-11 running all the ZHH data

```sh
ssh dev02

sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

dataset="ZHH"

archive_pipedir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/out
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/logs

cd ${sing_home}

cp ${outputdir}/../../ds000030_R1.0.5/out/ciftify_PINT/Yeo7_2011_80verts_roiidx_LUT.txt ${outputdir}/ciftify_PINT/

func_base="task-rest_bold"
for preprocfile in `ls ${archive_pipedir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_${func_base}_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
  if [ ! -f ${outputdir}/ciftify_PINT/${subject}/${subject}_${func_base}_atlas-pvertexNET_roi-7_fcmap.dscalar.nii ]; then
echo ${cleaning_script} ${subject} ${session} ${func_base} ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}  | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N dlabels_${subject}_${session} -j oe -o ${logsdir};
fi
done
```

```sh
ssh dev01

sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_PINT_dlabel_seedcorr.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

for dataset in "DTI3T"; do

archive_pipedir=/KIMEL/tigrlab/scratch/jjeyachandra/test_env/archive/data/${dataset}/pipelines/
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/logs

cd ${sing_home}

cp ${outputdir}/../../ds000030_R1.0.5/out/ciftify_PINT/Yeo7_2011_80verts_roiidx_LUT.txt ${outputdir}/ciftify_PINT/

func_base="task-rest_acq-CMH_run-01_bold"
for preprocfile in `ls ${archive_pipedir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_${func_base}_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))

  if [ ! -f ${outputdir}/ciftify_PINT/${subject}/${session}/${subject}_${session}_${func_base}_atlas-pvertexNET_roi-7_fcmap.dscalar.nii ]; then
echo ${cleaning_script} ${subject} ${session} ${func_base} ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container} #  | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N dlabels_${subject}_${session} -j oe -o ${logsdir};

fi
done
done
```
