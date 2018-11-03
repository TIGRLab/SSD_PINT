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
