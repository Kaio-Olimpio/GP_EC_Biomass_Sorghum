# Prediction of biomass sorghum hybrids using environmental covariate-enriched genomic combining ability models in tropical environments

This repository contains the data and scripts necessary to reproduce our analyses. Currently, the referred paper is being reviewed.

## The following files contains: 

-   Data: 
    -    pheno.rds: A .rds file containing the EBLUEs of the 221 hybrids (gen) evaluated at 64 field trials (env). This dataset also contains information on the geographical coordinates (lat and long), sowing and harvest dates (plant.date and harv.date), and experimental precision (heritability - her, coefficient of experimental variation - cv and mean - cv) of these trials. We also provide the parents of each hybrid (Aline and Rline), the standard error of the EBLUEs (std.error), the weights of the first-stage analysis (W), and the year (Ano), type (Ens, where F is for "Final" and P is for "Preliminary") and location (Loc) of each trial. The decoding of the "Loc" function is in Figure 1 of the paper.
    -   geno.rda: A .rda file with the three matrices containing molecular marker information used in the paper. `genoA` is the matrix filtered for the A-lines. `genoR` is the matrix with the same information for the R-lines. `genoF` has the complete information.
    -   envcov.csv: A .csv file containing information on the environmental covariates used in the study. More details about these covariates is in Table 2 of the paper.

-   Supporting:
    -   EC_download.R: A script describing how to download environmental covariates and extract them from rasters
    -   crossvalid.R: A script with two homemade functions designed to expedite the cross-validation process in the paper. More details on how to use the functions can be found at the description located above the function per se.

-   Script: A .R data with the codes used to generate the results of the paper.