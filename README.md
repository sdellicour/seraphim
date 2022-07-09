seraphim <img src="unix_OS/man/logo_seraphim.png" align="right" alt="" width="200" />
===============

`seraphim` is a R package for investigating the impact of environmental factors on the dispersal history and dynamics of viral lineages. The package can also be used to estimate dispersal statistics and mapping continuous phylogeographic trees.

## References
* Dellicour S, Rose R, Faria N, Lemey P, Pybus OG (2016). SERAPHIM: studying environmental rasters and phylogenetically-informed movements. _Bioinformatics_ 32: 3204-3206.
* Dellicour S, Rose R, Pybus OG (2016). Explaining the geographic spread of emerging epidemics: a framework for comparing viral phylogenies and environmental landscape data. _BMC Bioinformatics_ 17: 82.

## Stay tuned!
If you want to remain informed about last updates or improvements, just send an e-mail to simon.dellicour[at]ulb[dot]be with "seraphim mailing list" in the object.

## Installation
In R, `seraphim` can be installed with the `devtools` package:
```R
install.packages("devtools"); library(devtools)
install_github("sdellicour/seraphim/unix_OS") # (for a Unix OS)
install_github("sdellicour/seraphim/windows") # (for a Windows OS)
```

Using `seraphim` with the command line version of the program `Circuitscape` can sometimes be challenging as the later requires specific Python 2 settings. To facilitate its installation and subsequent application, users can now use `Conda` to install the appropriate environment specified in the `seraphim.yaml` file. `Conda` is an open source package/environment management system that can be installed following the instructions available [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once `Conda` is installed, the installation of the `seraphim` environment can be performed in a `Terminal` using the following command:
```
sudo conda env create -f seraphim.yaml --force
```
And the `Conda` environment can then be activated as follows:
```
conda activate seraphim
```
Within this environment, the command line version of `Circuitscape` running on Python 2 can subsequently be installed as follows:
```
sudo -H pip install -U --upgrade circuitscape
```
Always within the dedicated `Conda` environment, `seraphim` can finally be installed in R as described above:
```R
library(devtools)
install_github("sdellicour/seraphim/unix_OS")
```
