## Notes issues with CMH subject_list

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
