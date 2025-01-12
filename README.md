# Truncate and add Methyl caps to a receptor-ligand structure

## Description

This is a simple script to truncate a protein/receptor structured based on a distance criterium with respect to the present ligand.

## Getting Started

### Dependencies

As a prerequisite, the scripts relies on [AmberTools](https://ambermd.org/AmberTools.php) (`cpptraj`) to process the `.pdb` structure.

For the rest, the script is only reliant on commands, e.g. `awk` and `sed`, for text processing and manipulation which should be available in (most of the) Linux flavours 

###PYMOL specific

This version of the script requires `pymol`. `pymol` ([free version](https://anaconda.org/conda-forge/pymol-open-source) can easily be installed via `conda`:
```
$ conda install -c conda-forge pymol-open-source 
```

> [!IMPORTANT]  
> The script currently contains the activation of a specific `conda` environment where I installed `pymol`. Be sure to change this according to how/where you installed `pymol`.

### Executing program

The script is simply run from the CLI:
```
$ bash TruncAndCap.sh
```

This should produce the truncated and capped complex. Please note that the caps will have missing heavy atoms. These are added automatically when loading the `.pdb` in `leap` to prepare your parameter and topology files.

The script is currently written in a (very) static way. In several places i added a `MODIFY` comment where you might want to change some things to work on your structures. Pay specific attention to file names and paths, the ligand residue number and the distance threshold you want to set.

## Help

Feel free to reach out if you have trouble running the script.

Everyone is welcomed to branch and/or enhance the script as they feel.

## Versions

See `modifications.txt` for update history

## Authors

Charles-Alexandre Mattelaer
