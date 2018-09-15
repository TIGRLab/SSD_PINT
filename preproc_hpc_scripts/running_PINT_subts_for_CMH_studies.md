## My commands and notes for running the PINT and post-processing steps for CAMH studies

Started 2018-09-13

Before this..Jerry ran the ciftify container on particpants in his test space..

Now i will try to run PINT.

## rerunning PINT and cleaning for those last ZHH peeps

```sh
ssh dev02

template_out=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out

sing_home=/KIMEL/tigrlab/scratch/edickie/saba_PINT/sing_home
ciftify_container=/KIMEL/tigrlab/archive/code/containers/FMRIPREP_CIFTIFY/tigrlab_fmriprep_ciftify_1.1.2-2.0.9-2018-07-31-d0ccd31e74c5.img
cleaning_script=/KIMEL/tigrlab/projects/edickie/code/SZ_PINT/bin/participant_ciftify_clean_PINT_subcortts.sh

# projectname="ASDD"
# projectname="DTI3T"
projectname="PNSC"
# projectname="RTMSWM"
# projectname="SPINS"
archive_pipedir=/KIMEL/tigrlab/scratch/jjeyachandra/test_env/archive/data/${projectname}/pipelines/
outputdir=/KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${projectname}/out

mkdir /KIMEL/tigrlab/scratch/edickie/saba_PINT/ciftify_fmriprep/${projectname}/logs

## this only needs to be run once to set up the templates
mkdir -p ${outputdir}/ciftify_clean_img ${outputdir}/ciftify_PINT ${outputdir}/ciftify_meants/templates
cp ${template_out}/ciftify_clean_img/*.json ${outputdir}/ciftify_clean_img/
cp ${template_out}/ciftify_PINT/*.csv ${outputdir}/ciftify_PINT/
cp ${template_out}/ciftify_meants/templates/*.nii ${outputdir}/ciftify_meants/templates/

## now actually submitting the shiz
module load singularity/2.5.2
export OMP_NUM_THREADS=4


func_base="task-rest_acq-CMH_run-01_bold"
for preprocfile in `ls ${archive_pipedir}/fmriprep/sub-*/ses-*/func/sub-*_ses-*_${func_base}_space-T1w_preproc.nii.gz`; do
  subject=$(basename $(dirname $(dirname $(dirname ${preprocfile}))))
  session=$(basename $(dirname $(dirname ${preprocfile})))
  if [ ! -f ${outputdir}/ciftify_PINT/${subject}/${session}/${subject}_${session}_${func_base}_desc-clean_bold_summary.csv ]; then
echo ${cleaning_script} ${subject} ${session} ${func_base} ${archive_pipedir} ${outputdir} ${sing_home} ${ciftify_container}  | qsub -V -l walltime=00:20:00,nodes=1:ppn=4 -N pint_${subject}_${session} -j oe -o ${outputdir}/../../${projectname}/logs;
fi
done
```

ran section above with
