# MAPK_Experimental_Signature
This repository contains an example of the code to generate a MAPK Experimental Signature from the PlateSeq data deposited at GSE252002.
Briefly, the above cited data, is a PlateSeq drug screening performed on two different PDA cell lines (ASPC1 and PANC1, two plates per line) using 14 different RAF/MEK/ERK Inhibitors.

The aim of the R script, called *MAPK_Exp_Signature_Generation.R*, is to generate a consensus MAPK gene expression signature combining the mechanisms of action of all the tested drugs. 

The input files necessary to run the code are:

* the logTPM normalized counts (*ASPC1_logTPM_Expression_Final.rds* and *PANC1_logTPM_Expression_Final.rds*)
* and the metadata associated to each one of the 4 plates (*ULA-APSC1/PANC1-PS-A/B-_Final.rds*)

