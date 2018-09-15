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

## making separated ones into the templates repo

```sh
module load connectome-workbench
templates_dir=/scratch/edickie/saba_PINT/ciftify_fmriprep/ZHH/out/ciftify_meants/templates

wb_command -cifti-separate \
 rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
 COLUMN \
 -volume THALAMUS_LEFT rsn_THALAMUS_LEFT.nii.gz \
 -volume THALAMUS_RIGHT rsn_THALAMUS_RIGHT.nii.gz \
 -volume CEREBELLUM_LEFT rsn_CEREBELLUM_LEFT.nii.gz \
 -volume CEREBELLUM_RIGHT rsn_CEREBELLUM_RIGHT.nii.gz \
 -volume ACCUMBENS_LEFT rsn_ACCUMBENS_LEFT.nii.gz \
 -volume ACCUMBENS_RIGHT rsn_ACCUMBENS_RIGHT.nii.gz \
 -volume CAUDATE_LEFT rsn_CAUDATE_LEFT.nii.gz \
 -volume CAUDATE_RIGHT rsn_CAUDATE_RIGHT.nii.gz\
 -volume PUTAMEN_LEFT rsn_PUTAMEN_LEFT.nii.gz \
 -volume PUTAMEN_RIGHT rsn_PUTAMEN_RIGHT.nii.gz

wb_command -cifti-create-dense-from-template \
 rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
 ${templates_dir}/7RSN_roi-Lcerebellum_atlas.dlabel.nii \
 -volume CEREBELLUM_LEFT rsn_CEREBELLUM_LEFT.nii.gz

wb_command -cifti-create-dense-from-template \
  rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
  ${templates_dir}/7RSN_roi-Rcerebellum_atlas.dlabel.nii \
  -volume CEREBELLUM_RIGHT rsn_CEREBELLUM_RIGHT.nii.gz

wb_command -cifti-create-dense-from-template \
  rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
  ${templates_dir}/7RSN_roi-Lthalamus_atlas.dlabel.nii \
  -volume THALAMUS_LEFT rsn_THALAMUS_LEFT.nii.gz

wb_command -cifti-create-dense-from-template \
  rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
  ${templates_dir}/7RSN_roi-Rthalamus_atlas.dlabel.nii \
  -volume THALAMUS_RIGHT rsn_THALAMUS_RIGHT.nii.gz

wb_command -cifti-create-dense-from-template \
  rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
  ${templates_dir}/7RSN_roi-Lstriatum_atlas.dlabel.nii \
 -volume ACCUMBENS_LEFT rsn_ACCUMBENS_LEFT.nii.gz \
 -volume CAUDATE_LEFT rsn_CAUDATE_LEFT.nii.gz \
 -volume PUTAMEN_LEFT rsn_PUTAMEN_LEFT.nii.gz

wb_command -cifti-create-dense-from-template \
  rsn_yeo-cortex_buckner-cerebellum_choi-striatum_thalamus_7networks_networks_MWfix.dlabel.nii \
  ${templates_dir}/7RSN_roi-Rstriatum_atlas.dlabel.nii \
 -volume ACCUMBENS_RIGHT rsn_ACCUMBENS_RIGHT.nii.gz \
 -volume CAUDATE_RIGHT rsn_CAUDATE_RIGHT.nii.gz\
 -volume PUTAMEN_RIGHT rsn_PUTAMEN_RIGHT.nii.gz

rm rsn_THALAMUS_*
rm rsn_CEREBELLUM_*
rm rsn_ACCUMBENS_*
rm rsn_CAUDATE_*
rm rsn_PUTAMEN_*
```
