IMPORTANT NOTE - The code is compatible with the COBRA 2.0 release circa 2015, and may not be compatible with COBRA 3.0 due to major changes in the COBRA toolbox.
We suggest researchers developing context specific models to use the new version of the function createTissueSpecificModel.m available with COBRA 3.0 in order to ensure the usage of the most up-to-date extraction features

# Context_Specific_Models_from_GeMS
A systematic evaluation of methods for tailoring genome-scale metabolic models
DOI: Ref:

This repository provides a code base developed to benchmark different algorithms for generating context-specific metabolic network models. This contains code for generating models, processed RNA-Seq and genome-wide CRISPR-Cas9 gene knockout screen data, and code to use these data to compare the accuracy of the resulting models. This can be adapted for others to compare new algorithms and different parameter sets.


## Folder « Models and data »
Contains all the data and input models used to extract context specific models for the four cell lines (A375, HL60, KBM7 and K562)

## Folder « Extraction Methods »
Contains the code to run the different extraction algorithms

## Folder « Run analysis »
Code to run the gene essentiality and metabolic capabilities analysis

## Example to run the analysis on a model extracted for A375 using Fastcore:
[GenEss_score,Func_score] = Run_GenEss_Func('FastCore_CB1_A375','A375')
