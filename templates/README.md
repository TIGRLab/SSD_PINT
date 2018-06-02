all files here copied from /projects/gjacobs/Functional_PALM_PNC_analysis/Labels on 2017-11-21

## 2018-06-01 to separate out the thalamus ROIS

```sh
module load connectome-workbench

wb_command -cifti-separate \
 rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
 COLUMN \
 -volume THALAMUS_LEFT rsn_THALAMUS_LEFT.nii.gz \
 -volume THALAMUS_RIGHT rsn_THALAMUS_RIGHT.nii.gz \
 -volume CEREBELLUM_LEFT rsn_CEREBELLUM_LEFT.nii.gz \
 -volume CEREBELLUM_RIGHT rsn_CEREBELLUM_RIGHT.nii.gz

wb_command -cifti-create-dense-from-template \
 rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
 rsn_yeo-cortex_buckner-cerebellum_MWfix.dlabel.nii \
 -volume CEREBELLUM_LEFT rsn_CEREBELLUM_LEFT.nii.gz \
 -volume CEREBELLUM_RIGHT rsn_CEREBELLUM_RIGHT.nii.gz

 wb_command -cifti-create-dense-from-template \
  rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
  rsn_thalamus_7networks_networks_MWfix.dlabel.nii \
  -volume THALAMUS_LEFT rsn_THALAMUS_LEFT.nii.gz \
  -volume THALAMUS_RIGHT rsn_THALAMUS_RIGHT.nii.gz

rm rsn_THALAMUS_*
rm rsn_CEREBELLUM_*

```
