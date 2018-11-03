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
echo ${cleaning_script} ${subject} ${session} task-rest_bold ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}  | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N dlabels_${subject}_${session} -j oe -o ${outputdir}/../../ds000030_R1.0.5/logs;
fi
done
```
