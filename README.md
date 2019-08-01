
# Comparison with randomized controlled trials as a standard for evaluating instruments in Mendelian randomization

This repository contains the code associated with the letter by Gill et al, titled "Comparison with randomized controlled trials as a standard for evaluating instruments in Mendelian randomization". This letter was written in response to the article by [Walker et al](https://doi.org/10.1093/ije/dyz155), titled "Repurposing antihypertensive drugs for the prevention of Alzheimer’s disease: a Mendelian randomization study".

## Scripts

The repository consists of the following scripts:

| Script             | Description                                                                                                                                                                                                                                                                                   |
|--------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| perform_MR.R       | Code to calculate Mendelian randomization estimates of antihypertensive drug effects on coronary heart disease and stroke risk using the instruments from [Walker et al](https://doi.org/10.1093/ije/dyz155). These estimates are presented in Figure 1 in the letter.                        |
| calculate_power.R  | Code to calculate the minimum and maximum true effect estimates required to achieve 80% statistical power with a 5% type-I error rate for Mendelian randomization analyses investigating the effect of each considered antihypertensive drug class on coronary heart disease and stroke risk. These estimates are presented in Table 1 in the letter. |
| func_mRnd_binary.R | Function to calculate the power of Mendelian randomization analyses with binary outcomes. Adapted from [mRnd repository](https://github.com/kn3in/mRnd/blob/master/functions.R).                                                                                                              |

## Data

Three datasets are required to run these analyses. These can be accessed via the following links:

(1) The supplementary tables from Walker et al: [ije-2018-12-1541-File004.xlsx](https://doi.org/10.1093/ije/dyz155)

(2) The European ancestry GWAS results for the outcome 'any stroke' from MEGASTROKE: [MEGASTROKE.1.AS.EUR.out](http://www.megastroke.org/) 

(3) The RCT estimates, provided in this repository as a csv: [rct_estimates.csv](https://github.com/venexia/rct-instrument-comparison/blob/master/rct_estimates.csv)
