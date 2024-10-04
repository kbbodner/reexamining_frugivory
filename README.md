## Overview

This repository contains the data and R code necessary to run the analyses and produce the figures for the 
**Letter**, "Re-examining evidence for birds optimizing fruit size near their geographic limits". This **Letter** is a response to the previously published *Science* article, ["Birds optimize fruit size consumed near their geographic range limits"](https://www.science.org/doi/10.1126/science.adj1856) (Martins et al. 2024).

## Getting Started
The project includes one R script, ***model_evaluation.R*** and one R data file, ***dataframe_gape_scaled.RData***.

 ***model_evaluation.R*** calls in ***dataframe_gape_scaled.RData*** to create and evaluate a "best-performing" model and a null model. The script also contains code to produce residual plots and to perform tests related to model assumptions. 

### Dependencies

* Analysis was performed using R (R Core Team 2023) (version 4.3.1).
* R packages necessary to recreate the analysis include here (Müller 2020) (version 1.0.1), glmmTMB (Brooks et al 2017) (version 1.1.10) and DHARMa (Hartig 2022) (version 0.4.6), MuMIn (Bartoń 2024) (version 1.48.4), parallel (R Core Team 2023) (version 4.3.1) and performance (Lüdecke et al. 2021) (version 0.12.3).

