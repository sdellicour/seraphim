"seraphim" 2.0 <img src="unix_OS/man/logo_seraphim.png" align="right" alt="" width="200" />
===============

`seraphim` is a R package for studying phylogenetically informed movements. It can for instance be used to investigate the impact of environmental factors on the dispersal history and dynamics of viral lineages, to estimate lineage dispersal statistics, and to map continuous phylogeographic reconstructions.

## What's new in "seraphim" 2.0?
* spatio-temporal information embedded in annotated MCC or posterior trees retrieved from a Bayesian continuous phylogeographic inference can now be extracted using two new functions — `mccTreeExtractions`and `postTreeExtractions` — to generate similar extraction files, on per tree and with each row corresponding to a distinct phylogenetic branch.
* since the initial 2016 [application note](https://academic.oup.com/bioinformatics/article/32/20/3204/2196575?login=true) of the package, the previous `spreadGraphic1` function has been replaced by the `spreadGraphic2` function that lead to the generation of uncertainty polygons that can be saved as continuous vectorial files (e.g. in a shapefile format) instead of in a raster file (see the related tutorial [here](https://github.com/sdellicour/seraphim/blob/master/tutorials/plotting_the_dispersal_history.pdf)).
* the `spreadStatistics` function has been updated to now also include the estimation of isolation-by-distance (IBD) signal metrics (see our 2024 [study](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002914) as well as the dedicated [tutorial](https://github.com/sdellicour/seraphim/blob/master/tutorials/estimating_dispersal_statistics.pdf) for further detail).
* the `spreadFactors` function now focuses on testing the association between environmental factors on the diffusion - instead of the dispersal - velocity of lineages (see our 2025 [study]() as well as the dedicated [tutorial](https://github.com/sdellicour/seraphim/blob/master/tutorials/impact_on_diffusion_velocity.pdf) for further detail).
* the `spreadFactors` function can now also be used to conduct alternative post hoc analyses on the isolation-by-resistance (IBR), i.e. to what extent environmental factors can be associated with a deviation from an IBD pattern (see our 2025 [study]()).
* in addition to the two post hoc approaches implemented in the `spreadFactors` function, the package can now also be used to follow prior-informed landscape phylogeographic approaches to 

## Stay tuned!
If you want to remain informed about last updates or improvements, just send an e-mail to simon.dellicour[at]ulb[dot]be with "seraphim mailing list" in the object.

## Installation
In R, `seraphim` can be installed with the `devtools` package:
```R
install.packages("devtools"); library(devtools)
install_github("sdellicour/seraphim/unix_OS") # (for a Unix OS)
install_github("sdellicour/seraphim/windows") # (for a Windows OS)
```

## References
#### Package references
* Dellicour S, Faria N, Rose R, Lemey P, Pybus OG (_submitted_). SERAPHIM 2.0: an extended toolbox for studying phylogenetically informed movements
* Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. _Bioinformatics_ 32: 3204-3206

### Estimation of dispersal statistics
* Dellicour S, Bastide P, Rocu P, Fargette D, Hardy OJ, Suchard MA, Guindon S, Lemey P (2024). How fast are viruses spreading in the wild? _PLoS Biology_ 22: e3002914
* Dellicour S, Rose R, Pybus OG (2016). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. _BMC Bioinformatics_ 17: 82

### Landscape phylogeographic analyses
Investigating the association between lineage diffusion velocity and environmental factors:
* Dellicour S, Gámbaro F, Jacquot M, Lequime S, Baele G, Gilbert M, Pybus OG, Suchard MA, Lemey P (2025). Comparative performance of viral landscape phylogeography approaches. _Proceedings of the National Academy of Sciences of the USA_, _in press_
* Dellicour S, Rose R, Faria NR, Vieira LFP, Bourhy H, Gilbert M, Lemey P, Pybus OG (2017). Using viral gene sequences to compare and explain the heterogeneous spatial dynamics of virus epidemics. _Molecular Biology & Evolution_ 34: 2563-2571

Investigating the association between lineage dispersal locations and environmental factors:
* Dellicour S, Lequime S, Vrancken B, Gill MS, Bastide P, Gangavarapu K, Matteson NL, Tan Y, du Plessis L, Fisher AA, Nelson MI, Gilbert M, Suchard MA, Andersen KG, Grubaugh ND, Pybus OG, Lemey P (2020). Epidemiological hypothesis testing using a phylogeographic and phylodynamic framework. _Nature Communications_ 11: 5620
* Dellicour S, Troupin C, Jahanbakhsh F, Salama A, Massoudi S, Moghaddam MK, Baele G, Lemey P, Gholami A, Bourhy H (2019). Using phylogeographic approaches to analyse the dispersal history, velocity, and direction of viral lineages – application to rabies virus spread in Iran. _Molecular Ecology_ 28: 4335-4350

Investigating the association between lineage dispersal frequency and environmental factors:
* Dellicour S, Baele G, Dudas G, Faria NR, Pybus OG, Suchard MA, Rambaut A, Lemey P (2018). Phylodynamic assessment of intervention strategies for the West African Ebola virus outbreak. _Nature Communications_ 9: 2222

### Phylogeographic simulators
Simulation of a relaxed random walk (RRW) diffusion process along time-scaled phylogenies:
* 
Simulation of a relaxed random walk (RRW) diffusion process along time-scaled phylogenies:
* 






* 