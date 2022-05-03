
## Introduction

This set of function will allow you to process and analyse data from one or several .tsv files from ecotaxa.
This/these .tsv file(s) should be an export of a project containing images taken using the Planktoscope.
It will produce a first basic analysis of the diversity in your samples as well as a comparison between samples.

## Installation

You can install planktoscopeR from github with:

```
install.packages("devtools")
devtools::install_github("nhenry50/planktoscopeR")

```

## Usage

To have an example of how to use the functions, have a look at
`vignette("basic_analysis")`
Feel free to add any other tools that might be of interest to interpret classified
Planktoscope data.


## Definitions

- Biovolume : volume of a cell in mm^3 per cubic metre of water. The summed biovolume allows to represent proportion of a group of cells in a volume of water. 
- NBSS : Normalized Biovolume Size Spectra

