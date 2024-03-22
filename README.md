# DNA flexibility modules transcription factor binding to nucleosomes

## Prerequisite Softwares

Install the following softwares to run analyses:

1.  PEAR

Go to [link](https://cme.h-its.org/exelixis/web/software/pear/) for instructions of downloading and installing PEAR software. Update path to pear executable in `analyses/config.R`.

2.  DNAcycP

Go to directory `DNAcycP/` for instructions of building necessary dependencies and conda environment. Update `DNAcycP/run_cycp_tmp.sh` with environment name and slurm running variables.

3.  Cyclizability Prediction Website

Go to directory `Cyclizability-Prediction-Website/` and build a conda environment with the dependencies listed in `Cyclizability-Prediction-Website/requirements.txt`. 

## Prerequisite R Dependencies

ShortRead dplyr stringr glue SELEX rstatix matrixStats

## Reproducing the analyses

See `analysis/README.md`

## Reproducing the figures


## License

[MIT](https://choosealicense.com/licenses/mit/)