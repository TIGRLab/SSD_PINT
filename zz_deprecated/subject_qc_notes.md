## Notes issues with CMH subject_list

## anatomical Fails:
-all COBRE ok for anat
COBRE: 
COBRE_A00028806_SESS01 Really bad
ZHH anat issues:
EXP_21762_SESS01
EXP_21768_SESS01
EXP_22050_SESS01
EXP_22119_SESS01
EXP_22181_SESS01
EXP_22361_SESS01
EXP_22549_SESS01
EXP_22590_SESS01
EXP_21735_SESS01
EXP_21434_SESS01
EXP_21480_SESS01
EXP_21506_SESS01
EXP_21605_SESS01
EXP_21638_SESS01
EXP_21719_SESS01
EXP_21730_SESS01
EXP_21756_SESS01
EXP_21830_SESS01
EXP_21849_SESS01
EXP_21875_SESS01
EXP_21966_SESS01
EXP_22136_SESS01
EXP_22184_SESS01
EXP_22300_SESS01

## Some weird DMN conn in qc_fmri
funcfiles=`ls */MNINonLinear/Results/REST_0?/REST_0?.nii.gz`
for funcfile in $funcfiles; do dim4val=`fslval $funcfile dim4`; echo $dim4val $funcfile; done | sort

### all these COBRE runs are not acceptable (although mayber A00018129 could be concatenated?)
46 COBRE_A00014804_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
47 COBRE_A00031478_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
48 COBRE_A00012767_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
53 COBRE_A00037619_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
56 COBRE_A00018129_SESS01/MNINonLinear/Results/REST_02/REST_02.nii.gz
61 COBRE_A00018129_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
58 COBRE_A00028189_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
58 COBRE_A00028402_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
60 COBRE_A00028406_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
61 COBRE_A00023330_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
62 COBRE_A00023143_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
62 COBRE_A00037034_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
64 COBRE_A00031764_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
65 COBRE_A00000456_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
67 COBRE_A00033214_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
67 COBRE_A00035836_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
69 COBRE_A00013216_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
71 COBRE_A00028806_SESS01/MNINonLinear/Results/REST_01/REST_01.nii.gz
48 RTMSWM/hcp/RTMSWM_CMH_WM003_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
48 SPINS/hcp/SPN01_CMH_0033_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
55 RTMSWM/hcp/RTMSWM_CMH_WM002_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
57 RTMSWM/hcp/RTMSWM_CMH_WM015_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
60 RTMSWM/hcp/RTMSWM_CMH_WM019_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
63 DTI3T/hcp/DTI_CMH_S158_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
63 RTMSWM/hcp/RTMSWM_CMH_WM054_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
71 DTI3T/hcp/DTI_CMH_H127_02/MNINonLinear/Results/REST_01/REST_01.nii.gz


### These one..oddly are censored too...but nothing in censored that's over 100 left over..
77 RTMSWM/hcp/RTMSWM_CMH_WM018_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
79 RTMSWM/hcp/RTMSWM_CMH_WM063_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
79 SPINS/hcp/SPN01_CMH_0103_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
80 SPINS/hcp/SPN01_CMH_0050_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
89 SPINS/hcp/SPN01_CMH_0095_01/MNINonLinear/Results/REST_01/REST_01.nii.gz
91 DTI3T/hcp/DTI_CMH_S019_02/MNINonLinear/Results/REST_01/REST_01.nii.gz
95 DTI3T/hcp/DTI_CMH_H020_02/MNINonLinear/Results/REST_01/REST_01.nii.gz
98 RTMSWM/hcp/RTMSWM_CMH_WM013_01/MNINonLinear/Results/REST_01/REST_01.nii.gz


## I search for issued but running the following commands inside hcp/logs

```sh
 grep Resources * | grep walltime=01 ## will find people would timedout
 grep "does not exist" * # will identify error messages for missing inputs
```

## these subjects needed to be excluded because of their fMRI data quality
ASDD_CMH_HEF011_01 -> fMRI failed trsub
RTMSWM_CMH_WM017_01 -> fMRI failed trsub
RTMSWM_CMH_WM047_01 -> fMRI failed trsub

## Files missing where fs2hcp did not finish
DTI_CMH_S055_02
DTI_CMH_H036_02
DTI_CMH_S156_01
DTI_CMH_S171_01
SPN01_CMH_0004_01
PNS_CMH_0023_01


## subjects that timed out..
