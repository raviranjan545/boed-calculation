# Calculating Bayesian SIG-Optimal Experimental Designs

This repository contains interactive R scripts for exploring and calculating Bayesian SIG-optimal experimental designs across four different environmental driver models (nutrients, light, toxins and temperature). If you're not familiar with Github, just click [here](https://downgit.github.io/#/home?url=https://github.com/raviranjan545/boed-calculation) to download the entire repository. 

These scripts find experimental designs that maximize the expected Shannon Information Gain (SIG) from the prior to the posterior for non-linear models, and are associated with our [preprint](https://doi.org/10.64898/2026.05.28.728579). They use the Approximate Coordinate Exchange (ACE) algorithm via the `acebayes` package. 

Calculating Bayesian optimal designs can take a long time. We have chosen default values for the number of experimental units, a design space, priors, and additional settings that govern how long and how carefully the algorithm will search for an optimal design (not very long by default). We strongly recommend reading the overview of ACE in the [Supplementary Information](https://www.biorxiv.org/content/10.64898/2026.05.28.728579.supplementary-material) of our preprint and thinking about the settings most applicable to your experiment before running a fine-grained search for the best design.

Note that this repository is *not* meant to reproduce the analysis and figures in the preprint. The code to reproduce all preprint results and figures can be found [here](https://github.com/raviranjan545/boed-paper-reproduction). 

## Getting Started

## 1. Open the R Project

To ensure all file paths work correctly, **you must open this project using the `.Rproj` file.** 1. Download or clone this repository to your computer. 2. Double-click the `.Rproj` file to launch RStudio. 3. This will automatically lock your working directory to the project root, allowing the `here` package to navigate the folders seamlessly.

## 2. Dependencies

The exploratory scripts will automatically attempt to install missing packages when run. However, you can manually install the required dependencies by running:

``` r
install.packages(c("acebayes", "ggplot2", "dplyr", "tidyr", "here", 
                 "parallel", "patchwork", "scales", "GGally", "MASS", "Matrix"))
```

## 3. Repository Structure

- **`*.Rproj`**: The RStudio project file. Always open the project using this file.
- **Exploratory Scripts**:
  - `light_SIG_design_exploration.R` (Eilers-Peeters)
- `toxin_SIG_design_exploration.R` (Log-logistic)
- `temperature_SIG_design_exploration.R` (Norberg)
- `nutrients_SIG_design_exploration.R` (Monod)
- **Support Files**:
  - `stable_SIG_CPP_function.R` & `stable_utility_function.R`: Custom C++ and utility functions sourced by the scripts to prevent underflow issues during log-likelihood calculations.
- `temperature_prior.R`: Generates the correlated multivariate normal prior used specifically for the Norberg temperature model.

## 4. How to Use the Scripts

Each exploratory script is divided into sequential sections:

1.  **Setup:** Loads packages, sets the working directory via `here::i_am()`, and sources the stable calculation files.
2.  **Prior Visualization:** Draws samples from the specified prior distributions (either independent lognormal or correlated multivariate normal) and plots them.
3.  **Design Calculation:** Uses the `acenlm()` and `pacenlm()` functions to search a specified grid for the optimal 5-point design.
4.  **Result Visualization:** Plots the final optimal design.

### Important Note on Parallel Processing (Windows vs. Mac/Linux)

The `acebayes` package relies on fork-based parallel processing (`mclapply`) to rapidly evaluate multiple starting designs via the `pacenlm()` function. **This is not supported on Windows.**

To accommodate all users, the design calculation section of these scripts is split into two blocks: \* **Mac/Linux Users:** Can run the default `pacenlm()` code block. \* **Windows Users:** Must comment out the Mac/Linux block and uncomment the Windows block. The Windows block manually builds a PSOCK cluster (`parLapply`), exports the required packages and environments to the workers and reconstructs the pacenlm output object.

## Further reading

For an overview of Bayesian optimal experimental design for nonlinear regression, see our preprint:

[Ranjan, R., & Thomas, M. K. (2026). Bayesian optimal designs for common single-driver experiments in ecology. *bioRxiv*, 2026-05.](https://www.biorxiv.org/content/10.64898/2026.05.28.728579.abstract)

For a deeper technical explanation of the ACE algorithm, see the original paper:

[Overstall, A. M., Woods, D. C., & Parker, B. M. (2020). Bayesian optimal design for ordinary differential equation models with application in biological science. Journal of the American Statistical Association.](https://www.tandfonline.com/doi/full/10.1080/00401706.2016.1251495)
