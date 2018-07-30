#!/usr/bin/env python

'''Change the pixdim4 value and xyzt_units value in the header of a nifti file

Usage:
	EditHeader.py <Nifti file> <TR_value> <xyzt_units>

Argument(s):
	<Nifti file>	A nifti file
	<TR_value>	This is the TR that will be substituted in the field "pixdim4" in the header
	<xyzt_units>	This is the value to be substituted in xyzt_units field - it is an integer.

Written by Saba Shahab, July 2018
'''

import os
import sys
import nibabel

if len(sys.argv) < 4:
	print("Usage: EditHeader.py <Nifti file> <TR_value> <xyzt_units>")
	sys.exit()

#Load the image file
n1_img = nibabel.load(sys.argv[1])

# Functions
# This function edits the pixdim4 field in the header
def EditTR(filepath, TR):
	n1_img.header['pixdim'][4] = TR
	nibabel.save(n1_img, filepath)


def EditTRunits(filepath, TRunits):
	n1_img.header['xyzt_units'] = TRunits
	nibabel.save(n1_img, filepath)
	

if __name__ == "__main__":
	#EditTR(sys.argv[1], sys.argv[2])
	EditTRunits(sys.argv[1], sys.argv[3])
	print("Edited the header in the nifti file:", sys.argv[1])

