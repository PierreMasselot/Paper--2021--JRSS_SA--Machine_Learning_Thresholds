# Machine learning approaches to identify thresholds in a heat-health warning system context
---

R code and data reproducing results in the article:

Masselot P, Chebana F, Campagna C, Lavigne É, Ouarda TBMJ, Gosseling P. Machine learning approaches to identify thresholds in a heat-health warning system context. *Journal of the Royal Statistical Society: Series A (Statistics in Society)*. **In Press**.

## Description

Considers several methods to estimate data-driven thresholds for warning system. The goal of each method is to estimate the threshold directly from the association between the indicators and outcome of interest, looking for a break in the relationship.

Considered methods include:

- Model-based partitioning (MOB)
- Multivariate adaptive regression splines (MARS)
- Patient rule-induction method (PRIM)
- Adaptive index models (AIM)

They are compared in a simulation study and applied to a real-world case of study of heat thresholds in Montréal, Canada. Results suggest that PRIM and MOB give the best thresholds.

## Reproducible R code

The R code and data can be used to replicate the results and adapt the methods to new applications. Scripts are classified into three categories:

- Custom functions used in the code:

  *00_Misc_functions.R* Implements several convenience functions used throughout the analysis
  
  *01_Threshold_functions.R* Wrapper function for each method. Used in simulations
  
- Simulation study:

  *10_Simulations_functions.R* Functions used to generate synthetic data
  
  *11_Simulations.R* Main script implementing the simulations study, calling the generation function and applying each method wrapper
  
  *12_Simulation_plots.R* Uses results from *11_Simulations.R* to produce Figures found in the manuscript
  
- Real-world application:

  *21_Application.R* Script detailing the application of each method on data from Montréal. Produces Figures 2 to 5
  
  *22_Bootstrap.R* Script implementing bootstrap resampling to assess the bias and variability of each method
  
  *23_Plots.R* Uses results from *21_Application.R* and *22_Bootstrap.R* to produce Figures 6 and 7 as well as Table 2

In addition to the main script above, a contributed R package has been created to implement the PRIM method used in the paper. The package includes several methods to perform peeling, pasting and analyze results. It is for now only [on github](https://github.com/PierreMasselot/primr) but I plan to submit to CRAN (hopefully) soon.

Other methods applied and compared here were already implemented in contributed R packages, namely:

- `partykit` for MOB
- `earth` for MARS
- `AIM` for AIM
