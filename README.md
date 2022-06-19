seraphim <img src="unix_OS/man/logo_seraphim.png" align="right" alt="" width="200" />
===============

`seraphim` is a R package for investigating the impact of environmental factors on the dispersal history and dynamics of viral lineages. The package can also be used to estimate dispersal statistics and mapping continuous phylogeographic trees.

## Installation
In R, `seraphim` can be installed with the [`devtools`](https://github.com/hadley/devtools) package:
```R
install.packages("devtools"); library(devtools)
install_github("sdellicour/seraphim/unix_OS") # (for a Unix OS)
install_github("sdellicour/seraphim/unix_OS") # (for a Windows OS)
```

Using `seraphim` with the command line version of `Circuitscape` can sometimes be challenging as the later requires specific Python 2 settings. To facilitate its installation and subsequent application, users can now use `Conda` to install the appropriate environment specified in the `seraphim.yaml` file. This can be done in a `Terminal` with the following command:
```
sudo conda env create -f seraphim.yaml --force
```
Once installed, the `Conda` environment can be activated as follows:
```
conda activate seraphim
```
Within this environment, the command line version of `Circuitscape` running on Python 2 can then be installed as follows:
```
sudo -H pip install -U --upgrade circuitscape
```
Always within the dedicated `Conda` environment, `seraphim` can be installed in R as described above:
```R
library(devtools)
install_github("sdellicour/seraphim/unix_OS")
```

## References
* Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. _Bioinformatics_ 32: 3204-3206.
* Dellicour S, Rose R, Pybus OG (2016). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. _BMC Bioinformatics_ 17: 82.
