---
title: "README"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

This set of function will allow you to process and analyse data from a .tsv file from ecotaxa.
This .tsv file should be an export of a project containing images taken using the Planktoscope.
TODO define which export mode to use
      .tsv file containing different samples? should we merge tsv files from unique samples?

## Functions

- `import.table` is mandatory to use the functions provided in this package. This function will allow you to process the original tsv file and copy the data in different 'data.table' to be then analysed by specific functions

- `ellipsoid.biovol` calculates the biovolume based on the ellipsoid approximation of the volume of a particle. It uses the major and minor axis of the image as well as the pixel size. units are in [mm^3]

- `norm_biovol` calculates the normalized biovolume of a particle. This means that it takes into account the volume filtered by the sampling net. units are in [mm^3/m^3]

- `summed.biovol` computes the sum of biolovolumes for one sample TODO should it sum over biovol or norm_biovol

- `relative.abundance` computes the relative abundance of each taxonomic class. This is done by dividing the absolute count of a certain category by the total number of living objects found in the sample. The result is given in percentage.

- the `NBSS` functions allows to create a data.table containing the Normalized biovolume size spectra. This data.table can then be used to plot the spectra vs the normalized biovolumes.


