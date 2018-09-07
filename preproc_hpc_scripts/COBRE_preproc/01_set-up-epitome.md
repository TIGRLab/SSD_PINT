# Setting up the COBRE data for epitome analyses

Run by Erin W Dickie, November 2017

# 1. we found all the T1w and Resting State Runs

We used the code below
Note: this pauses with a message to "kill with Ctrl-D" when it can't find any images
(I just hit Ctrl-D a bunch o times)

```sh
working_dir=/external/SchizConnect/COBRE

cd ${working_dir}

mkdir -p epitome/COBRE
subjects=`cd data/nii; ls -1d A00*`

for subj in ${subjects}; do
    fmri_data=$(find data/nii/${subj} -name Resting_State*.nii.gz)
    parallel "mkdir -p ${working_dir}/epitome/COBRE/${subj}/REST/SESS01/RUN0{#}" ::: $fmri_data
    parallel "cd ${working_dir}/epitome/COBRE/${subj}/REST/SESS01/RUN0{#}; ln -s ../../../../../../{} {/}" ::: $fmri_data

    t1_data=$(find data/nii/${subj} -name MPRAGE_??_??*.nii.gz)
    parallel "mkdir -p ${working_dir}/epitome/COBRE/${subj}/T1/SESS01/RUN0{#}" ::: $t1_data
    parallel "cd ${working_dir}/epitome/COBRE/${subj}/T1/SESS01/RUN0{#}; ln -s ../../../../../../{} {/}" ::: $t1_data
done
```

## Excluding obviously bad scans pre epitome.

#### For Resting_State

+ I am only excluding scans where the number of volumes is less than 150
+ two subjects fit this criteria

A00020602/RST/SESS01/RUN01/Resting_State_12_01.nii.gz 150
A00020602/RST/SESS01/RUN02/Resting_State_13_01.nii.gz 89
A00020602/RST/SESS01/RUN03/Resting_State_14_01.nii.gz 150
A00034381/RST/SESS01/RUN01/Resting_State_12_01.nii.gz 150
A00034381/RST/SESS01/RUN02/Resting_State_13_01.nii.gz 62

```sh
## this subject I am moving run3 data to run2 after deleting run 2
rm -r A00020602/REST/SESS01/RUN02
mv A00020602/REST/SESS01/RUN03 A00020602/REST/SESS01/RUN02
## this subject I am just deleting run 2
rm -r A00034381/REST/SESS01/RUN02ls
```
#### For T1w Some File have many

+ looking at the slicesdir pics..some of the images from RUN03 on look terrible
+ so I will only keep RUN01 and RUN02 for everyone

```sh
cd /external/SchizConnect/COBRE/epitome/COBRE
rm -r A000*/T1/SESS01/RUN03
rm -r A000*/T1/SESS01/RUN04
rm -r A000*/T1/SESS01/RUN05
rm -r A000*/T1/SESS01/RUN06
```

# AND NOW for the magik!!! Submit epitome jobs

```sh
cd /external/SchizConnect/COBRE/epitome/COBRE

cp /projects/edickie/code/SZ_PINT/COBRE_preproc/epitome* ./


source epitome_env.sh
chmod +x epitome_run.sh
parallel "echo ${PWD}/epitome_run.sh {}" ::: $(ls -1d A0*) | \
  qbatch -b sge -N epitome_COBRE_run2 -
parallel "echo ${PWD}/epitome_run.sh {}" ::: $(ls -1d A0*) | \
  qbatch -b sge -N epitome_COBRE_run3 -

## because the queue was acting up.. and I want a record .. just sending them in serial on frankin
for subid in $(ls -1d A0*); do
  ./epitome_run.sh ${subid}
done | tee run_epitome_serial.log
```
