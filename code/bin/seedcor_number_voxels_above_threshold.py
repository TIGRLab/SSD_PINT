#!/usr/bin/env python3
"""
For each network, will count the number of vertices in a seedcorr output that
is above a range of thresholds.

Usage:
  seedcor_number_voxels_above_threshold [options] <output_base_dir> <pheno_qa_pass_csv> <output_csv>

Arguments:
    <output_base_dir>     The base output directorty with the PINT and seedcorr outputs
    <pheno_qa_pass_csv>   The qa pass csv input file (the filepaths are parsed from this)
    <output_csv>          The path to the output file

Options:
  --thres-start start     Starthing threshold [default: 0.3] count vertices above.
  --thres-end end         Ending threshold [default: 0.7] in range.
  --thres-step step       The step value [default: 0.02] in the threshold range
  --debug                 Debug logging in Erin's very verbose style
  --help                  Print help

DETAILS

Written by Erin W Dickie, November 21, 2018
"""
# coding: utf-8

# In[28]:


import ciftify.niio
import os
import seaborn
import numpy as np
import pandas as pd
from docopt import docopt



def count_num_above_threshold(input_array, threshold):
    '''count the number of input array values above a threshold'''
    result = len(np.where( input_array > threshold)[0])
    return result



# In[137]:
def get_dscalar_filename(dataset, subject, session, summary_filename, network, outbase):
    '''create the filename of the dscalar file with the seed connectivity values
    using the other info in the dataframe'''
    fc_filename = summary_filename.replace('_desc-clean_bold_summary.csv',
                                          '_atlas-pvertexNET_roi-{}_fcmap.dscalar.nii'.format(network))
    if session:
        fc_path = os.path.join(outbase, dataset, "out", "ciftify_PINT", subject, session, fc_filename)
    else:
        fc_path = os.path.join(outbase, dataset, "out", "ciftify_PINT", subject, fc_filename)
    return fc_path


def main():
    arguments       = docopt(__doc__)
    outbase = arguments['<output_base_dir>']
    qapass_pheno_file = arguments['<pheno_qa_pass_csv>']
    output_csv_path = arguments['<output_csv>']
    debug = arguments['--debug']
    thresholds = np.arange(float(arguments['--thres-start']),
                        float(arguments['--thres-end']),
                        float(arguments['--thres-step']))

    # example_fcmap = '/scratch/edickie/saba_PINT/ciftify_fmriprep/ds000030_R1.0.5/out/ciftify_PINT/sub-10159/sub-10159_task-rest_bold_atlas-pvertexNET_roi-2_fcmap.dscalar.nii'
    # fc_map_data = ciftify.niio.load_cifti(example_fcmap)
    # thresholds = np.arange(0.3, 0.72, 0.02)
    # test_out = np.zeros(len(thresholds))
    # for i in range(len(thresholds)):
    #     test_out[i] = count_num_above_threshold(fc_map_data, thresholds[i])
    # test_out


    # read in the subject list
    subject_list = pd.read_csv(qapass_pheno_file)
    subdf = subject_list[['dataset','subject',"session", "filename"]]
    subdf6 = pd.concat([subdf]*6)

    ## create a thresholds vector and output template
    threscols = []
    for thres in thresholds:
        threscols.append('gt_{0:.2f}'.format(thres))
    thresdf = pd.DataFrame(columns = threscols, dtype=int, index=range(len(subdf6)))

    ## create a NETWORK column
    net_df = pd.DataFrame({"NETWORK":[2,3,4,5,6,7]*len(subdf)})

    ## put all of these together
    df = pd.concat([subdf6.sort_values('subject').reset_index(), net_df.reset_index(), thresdf], axis = 1)

    ## loop over everything and calculate
    for i in df.index:
        if debug: print(df.subject[i])
        session = df.session[i] if type(df.session[i])==str else None
        fc_path = get_dscalar_filename(df.dataset[i], df.subject[i], session, df.filename[i], df.NETWORK[i], outbase)
        try:
            fc_map_data = ciftify.niio.load_cifti(fc_path)
            for threshold in thresholds:
                df.loc[i,'gt_{0:.2f}'.format(threshold)] = count_num_above_threshold(fc_map_data, threshold)
        except:
            print('{} failed to load'.format(fc_path))

    df.to_csv(output_csv_path)

if __name__ == "__main__":
    main()
