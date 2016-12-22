# Context_Specific_Models_from_GeMS
A systematic evaluation of methods for tailoring genome-scale metabolic models
DOI: Ref:

## Folder « Models and data »
Contains all the data and input models used to extract context specific models for the four cell lines (A375, HL60, KBM7 and K562)

## Folder « Extraction Methods »
Contains the code to run the different extraction algorithms

## Folder « Run analysis »
Code to run the gene essentiality and metabolic capabilities analysis

Example to run the analysis on a model extracted for A375 using Fastcore:
[GenEss_score,Func_score] = Run_GenEss_Func('FastCore_CB1_A375','A375')
