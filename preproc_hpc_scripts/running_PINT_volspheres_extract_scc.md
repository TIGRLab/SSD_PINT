# Extracting meants from spheres in volume space

## 2019-04-02

If occurs to me that I might not be getting the point across that my method is an
improvement over traditional methods because even the template ROIs are not "traditional".
In other words, they may be an improvement over traditional volume based preprossing and extraction.
So we are adding an extra "comparison" timeseries that are 6mm volume spheres are the location of the
template ROIs.

To extract from these (fairly) we need to also clean a volume representation of the data.

The script to do this is in:
`bin/participant_ciftify_volclean_and_meants.sh`

It has inputs:
```
subject=$1
session=$2
func_base=$3
archive_pipedir=$4
outdir=$5
sing_home=$6
ciftify_container=$7
```

Now to try to run it on all the old data. (note: the local cluster has upgraded to slurm since I started this analysis)

```sh
ssh dev01
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out
archive_pipedir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/logs
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab//scratch/edickie/tigrlab_fmriprep_ciftify_1.3.0.post2-2.3.1-2019-04-04-8ebe3500bebf.img
# ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_ciftify_volclean_and_meants.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

cd ${logsdir}

for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/sub-*_task-rest_*_summary.csv`; do
  subject=$(basename $(dirname ${preprocfile}))
  session="none"
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${subject}_task-rest_bold_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv ]; then
sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=vspheres_${subject}_${session} ${cleaning_script} ${subject} ${session} task-rest_bold ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container} sbatch time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=vspheres_${subject}_${session}
fi
done
```

...and for COBRE (note the container was copied from my scratch to the archive)

```sh
ssh dev01
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/out
archive_pipedir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/COBRE/logs
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.3.0.post2-2.3.1-2019-04-04-8ebe3500bebf.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_ciftify_volclean_and_meants.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

cd ${logsdir}

cp ${outputdir}/../../ds000030_R1.0.5/out/ciftify_meants/templates/tpl-MNI152NLin6Asym_res-02_desc-6mmYeo780_dseg.nii.gz ${outputdir}/ciftify_meants/templates/

for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_task-rest_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_bold_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv ]; then
sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=vspheres_${subject}_${session} ${cleaning_script} ${subject} ${session} task-rest_bold ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}
fi
done

for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_task-rest_run-01_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_run-01_bold_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv ]; then
sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=vspheres_${subject}_${session} ${cleaning_script} ${subject} ${session} task-rest_run-01_bold ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}
fi
done

for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_task-rest_run-02_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_run-02_bold_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv ]; then
sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=vspheres_${subject}_${session} ${cleaning_script} ${subject} ${session} task-rest_run-02_bold ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}
fi
done
```

```sh
ssh dev01
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out
archive_pipedir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/logs
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.3.0.post2-2.3.1-2019-04-04-8ebe3500bebf.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_ciftify_volclean_and_meants.sh

module load singularity/2.5.2
export OMP_NUM_THREADS=4

cd ${logsdir}

cp ${outputdir}/../../ds000030_R1.0.5/out/ciftify_meants/templates/tpl-MNI152NLin6Asym_res-02_desc-6mmYeo780_dseg.nii.gz ${outputdir}/ciftify_meants/templates/

for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_task-rest_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_bold_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv ]; then
sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=vspheres_${subject}_${session} ${cleaning_script} ${subject} ${session} task-rest_bold ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}
fi
done
```

```sh
ssh dev01

module load singularity/2.5.2
export OMP_NUM_THREADS=4

sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.3.0.post2-2.3.1-2019-04-04-8ebe3500bebf.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_ciftify_volclean_and_meants.sh

for dataset in "DTI3T" "SPINS" "RTMSWM" "PNSC" "ASDD"; do

archive_pipedir=/KIMEL/tigrlab/scratch/jjeyachandra/test_env/archive/data/${dataset}/pipelines/
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/logs

cd ${logsdir}

cp ${outputdir}/../../ds000030_R1.0.5/out/ciftify_meants/templates/tpl-MNI152NLin6Asym_res-02_desc-6mmYeo780_dseg.nii.gz ${outputdir}/ciftify_meants/templates/

func_base="task-rest_acq-CMH_run-01_bold"
for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_${func_base}_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))

  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_${func_base}_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv ]; then
sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=vspheres_${subject}_${session} ${cleaning_script} ${subject} ${session} ${func_base} ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}
fi
done
done

```
