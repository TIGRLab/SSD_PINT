# Extracting Schaefer volume and surface versions

## 2023-20-02


A reviewer has asked to see how this all compares to Schaefer - so we will run Schaefer (7 networks, 100 ROI version) in volume and surface and recreate the results to see.

We will use the same cleaned volume and surface files

The script to do this is in:
`bin/participant_ciftify_volclean_and_meants.sh`

It has inputs:
```
subject=$1
session=$2
func_base=$3
outdir=$4
sing_home=$5
ciftify_container=$6
```

Now to try to run it on all the old data. (note: the local cluster has upgraded to slurm since I started this analysis)

```sh
ssh dev01
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out
archive_pipedir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/logs
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
#ciftify_container=/KIMEL/tigrlab//scratch/edickie/tigrlab_fmriprep_ciftify_1.3.0.post2-2.3.1-2019-04-04-8ebe3500bebf.img
# ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.1.0-2018-10-12-dcfba6cc0add.img

ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_v1.3.2-2.3.3-2019-08-16-c0fcb37f1b56.simg

cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/code/bin/participants_Schaefer_meants.sh

# module load singularity
export OMP_NUM_THREADS=4

# Note that the subcortical templates are generated from the main template using the commands in the templates/README.md from this repo..
cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz ${outputdir}/ciftify_meants/templates
cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order.dlabel.nii ${outputdir}/ciftify_meants/templates


cd ${logsdir}

for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/sub-*_task-rest_*_summary.csv`; do
  subject=$(basename $(dirname ${preprocfile}))
  session="none"
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${subject}_task-rest_bold_desc-clean_atlas-Shaefer7N100P_timeseries.csv ]; then
sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=schaefer_${subject}_${session} ${cleaning_script} ${subject} ${session} task-rest_bold ${outputdir} ${sing_home} ${ciftify_container} 
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
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_v1.3.2-2.3.3-2019-08-16-c0fcb37f1b56.simg

cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/code/bin/participants_Schaefer_meants.sh

# module load singularity
export OMP_NUM_THREADS=4

# Note that the subcortical templates are generated from the main template using the commands in the templates/README.md from this repo..
cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz ${outputdir}/ciftify_meants/templates
cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order.dlabel.nii ${outputdir}/ciftify_meants/templates


cd ${logsdir}
for func_base in "task-rest_bold" "task-rest_run-01_bold" "task-rest_run-02_bold"; do
for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_${func_base}_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${subject}_${session}_${func_base}_desc-clean_atlas-Shaefer7N100P_timeseries.csv ]; then
  sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=schaefer_${subject}_${session} ${cleaning_script} ${subject} ${session} ${func_base} ${outputdir} ${sing_home} ${ciftify_container} 
fi
done
done


```

```sh
ssh dev01
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out
archive_pipedir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/logs
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_v1.3.2-2.3.3-2019-08-16-c0fcb37f1b56.simg

cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/code/bin/participants_Schaefer_meants.sh

# module load singularity
export OMP_NUM_THREADS=4

# Note that the subcortical templates are generated from the main template using the commands in the templates/README.md from this repo..
cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz ${outputdir}/ciftify_meants/templates
cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order.dlabel.nii ${outputdir}/ciftify_meants/templates


cd ${logsdir}

for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_task-rest_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))
  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_task-rest_bold_desc-clean_atlas-Shaefer7N100P_timeseries.csv ]; then
    sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=schaefer_${subject}_${session} ${cleaning_script} ${subject} ${session} task-rest_bold ${outputdir} ${sing_home} ${ciftify_container}
fi
done
```

```sh
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/code/bin/participants_Schaefer_meants.sh
sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_v1.3.2-2.3.3-2019-08-16-c0fcb37f1b56.simg

for dataset in "DTI3T" "SPINS" "RTMSWM" "PNSC" "ASDD"; do

archive_pipedir=/KIMEL/tigrlab/scratch/jjeyachandra/test_env/archive/data/${dataset}/pipelines/
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/out
logsdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${dataset}/logs

cd ${logsdir}

cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz ${outputdir}/ciftify_meants/templates
cp /projects/edickie/code/SZ_PINT/data/raw/templates/Schaefer2018_100Parcels_7Networks_order.dlabel.nii ${outputdir}/ciftify_meants/templates

for func_base in "task-rest_acq-CMH_bold" "task-rest_acq-CMH_run-01_bold" "task-rest_acq-CMH_run-02_bold"; do
for preprocfile in `ls ${outputdir}/ciftify_PINT/sub-*/ses-*/sub-*_${func_base}_*_summary.csv`; do
  subject=$(basename $(dirname $(dirname ${preprocfile})))
  session=$(basename $(dirname ${preprocfile}))

  if [ ! -f ${outputdir}/ciftify_meants/${subject}/${session}/${subject}_${session}_${func_base}_desc-clean_atlas-Shaefer7N100P_timeseries.csv ]; then
    echo sbatch --time=00:20:00 --cpus-per-task=4 --nodes=1 --job-name=schaefer_${subject}_${session} ${cleaning_script} ${subject} ${session} ${func_base} ${outputdir} ${sing_home} ${ciftify_container}
fi
done
done
done

```
