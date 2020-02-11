# from quarantine
module load matlab/R2014a
module load FSL/5.0.10
module load FIX/1.061
module load R/3.2.4
module load R-extras/3.2.4
module load AFNI/2014.12.16
module load freesurfer/6.0.0
module load connectome-workbench/1.2.3

# load home-spun tools
module load /archive/code/packages.module
export PATH=$PATH:/projects/edickie/code/ciftify/ciftify/bin
export PYTHONPATH=$PYTHONPATH:/projects/edickie/code/ciftify
export CIFTIFY_TEMPLATES=/projects/edickie/code/ciftify/data

# experiment specific
export SCRIPTIT_DATA='/external/SchizConnect/COBRE/epitome'
export SUBJECTS_DIR='/external/SchizConnect/COBRE/freesurfer'
export HCP_DATA='/external/SchizConnect/COBRE/hcp'
