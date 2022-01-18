
## Introduction

This set of function will allow you to process and analyse data from a .tsv file from ecotaxa.
This .tsv file should be an export of a project containing images taken using the Planktoscope.
TODO define which export mode to use
      .tsv file containing different samples? should we merge tsv files from unique samples?


## TODO
- allow user to chose the subset of images to process
- sur fichier à part : faire PCA à partir de biovolumes 
- Roxygen
-

## Definitions

- Biovolume : volume of a cell in mm^3 per cubic metre of water. The summed biovolume allows to represent proportion of a group of cells in a volume of water. 
- NBSS : Normalized Biovolume Size Spectra

## Functions

- `import.table` is mandatory to use the functions provided in this package. This function will allow you to process the original tsv file and copy the data in different 'data.table' to be then analysed by specific functions. Within this function, the functions from biovolume and relative abundance are used to add these informations to the tables

- `ellipsoid.vol` calculates the biovolume based on the ellipsoid approximation of the volume of a particle. It uses the major and minor axis of the image as well as the pixel size. units are in [mm^3]

- `biovolume` calculates the normalized biovolume of a particle. This means that it takes into account the volume filtered by the sampling net. units are in [mm^3/m^3]

- `summed.biovol` computes the sum of biovolumes over one sample 

- `relative.abundance` computes the relative abundance of each taxonomic class. This is done by dividing the absolute count of a certain category by the total number of objects found in the sample. The result is given in percentage.

- the `NBSS` functions allows to create a data.table containing the Normalized biovolume size spectra. This data.table can then be used to plot the spectra vs the normalized biovolumes.


