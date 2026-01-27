# GHG-depths Analysis

This repository contains code used in a manuscript about the drivers of greenhouse gas concentrations in lakes worldwide.

## Contact

-   Abby Lewis ([lewisa4\@si.edu](mailto:lewisa4@si.edu))
-   Joe Rabaey ([jsr339\@cornell.edu](mailto:jsr339@cornell.edu))

## Data

All data are sourced directly from the GHG-depths dataset (<https://doi.org/10.6073/pasta/2b72b89bbfbb3da0e198f392a9cbad18>; Rabaey et al. 2026). GHG-depths is a geographically extensive dataset of depth-profile CO<sub>2</sub> and CH<sub>4</sub> concentrations, along with environmental drivers, from 522 lakes. These data span 7 continents and 38 countries, with a total of 2,558 discrete sampling events. The GHG-depths dataset is released for noncommercial use only and is licensed under a Creative Commons Attribution 4.0 International License (CC BY 4.0). All publications that use GHG-depths are encouraged to appropriately cite the data and corresponding paper, and to contact and collaborate with the site-specific data providers for expertise in interpreting data when appropriate.

## File structure

All analysis code is located in `./DataAnalysisScripts`. A compiled dataset of surface and bottom water data is created by `./DataAnalysisScripts/01_Surf_bot_sum.Rmd` and output in `./CompiledData`. Other analysis code files can be run in any order. Figures created by these analysis files are archived in `./Figures`.

## 
