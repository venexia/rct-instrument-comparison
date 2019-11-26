
# Comparison with randomized controlled trials as a standard for evaluating instruments in Mendelian randomization

This repository contains the code associated with the letter by Gill et al, titled "[Comparison with randomized controlled trials as a standard for evaluating instruments in Mendelian randomization](https://doi.org/10.1093/ije/dyz236)". 

This letter was written in response to the article by Walker et al, titled "[Repurposing antihypertensive drugs for the prevention of Alzheimer’s disease: a Mendelian randomization study](https://doi.org/10.1093/ije/dyz155)".

## Scripts

The script 'perform_MR.R' provides the code to calculate Mendelian randomization estimates of antihypertensive drug effects on coronary heart disease and stroke risk using the instruments from Walker et al in R. These estimates are presented in Figure 1 in the letter.

## Data

Three datasets are required to run these analyses. These can be accessed via the following links:

(1) The supplementary tables from Walker et al: [ije-2018-12-1541-File004.xlsx](https://doi.org/10.1093/ije/dyz155)

(2) The European ancestry GWAS results for the outcome 'any stroke' from MEGASTROKE: [MEGASTROKE.1.AS.EUR.out](http://www.megastroke.org/) 

(3) The RCT estimates, provided in this repository as a csv: [rct_estimates.csv](https://github.com/venexia/rct-instrument-comparison/blob/master/rct_estimates.csv)
