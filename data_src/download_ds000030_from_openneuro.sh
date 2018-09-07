mkdir /scratch/edickie/saba_PINT/bids_in/ds000030

chuncks='sub10159-10299
sub10304-10388
sub10428-10575
sub10624-10697
sub10704-10893
sub10912-10998
sub11019-11098
sub11104-11156
sub50004-50008
sub50010-50038
sub50043-50059
sub50060-50085'

for chunck_name in ${chuncks}; do
if [ ! -f ds000030_R1.0.5_${chunck_name}.zip ]; then
  wget https://s3.amazonaws.com/openneuro/ds000030/ds000030_R1.0.5/compressed/ds000030_R1.0.5_${chunck_name}.zip
fi
done

for chunck_name in ${chuncks}; do
if [ ! -f ds000030_R1.0.5_derivatives_freesurfer_${chunck_name}.zip ]; then
wget https://s3.amazonaws.com/openneuro/ds000030/ds000030_R1.0.5/compressed/ds000030_R1.0.5_derivatives_freesurfer_${chunck_name}.zip
fi
done
