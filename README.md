Companion code and data to:
Cerliani L., Mennes M., Thomas RM, DI Martino A, Thioux M. Keysers C. (2015) "Increased functional connectivity between subcortical and cortical resting-state networks in autism spectrum disorder"; JAMA Psychiatry, in press.

This distribution provides the data and the code do_FNC_group_differences.m to reproduce the analyses of differences in Functional Network Connectivity (FNC) described in the published manuscript. Please read the commented lines at the beginning of the script - and the disclaimer below - before using it.

The data.zip file needs to be unzipped in order to run the script do_FNC_group_differences.m.
FSL (http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) must be installed

DISCLAIMER
The software is provided 'AS IS'. The public distribution of this script
is UNIQUELY intended to provide an accurate documentation of the analyses 
performed in the published manuscript. The author declines any
responsibility regarding machine malfunctioning or data loss due to the
use of this script in the version provided here or in any modification
of it.

CONTENT

do_FNC_group_differences.m
Main script to reproduce the analyses. Please read the commented lines at the beginning of the script before using it.

data.zip
The data needed for the analyses (needs to be unzipped)

data/dr_stage1
The output of the first (spatial) regression using the output of metaICA as templates

data/list_359_IDs
The subject IDs of the participants included in the analyses.

data/pheno_359.mat
All the phenotypical data provided in the ABIDE document (Phenotypic_V1_0b.csv) available on the ABIDE database

CC_bptf_published.mat
The values of FNC used for the results reported in the paper (see more information inside the matlab script)

addons
Various supporting software. See also the commented lines at the beginning of the main script

metaICA_results
Self-explanatory. See the README inside the folder for more information.

###
Leonardo Cerliani
June 7th, 2015








