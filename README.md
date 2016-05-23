# A Pipeline for Predixcan 

To know about PrediXcan, please refer the [link](https://github.com/hakyim/PrediXcan).

## Key things to know before using this script

The script is only compatible with

1. the imputation dosage files produced by ricopilli pipeline (link) and it assumes a file structure derived from ricopilli pipeline. 

2. a linux environment with SLURM job sheduler


## This script does four things

1. Converts genotype dosage format ([refer here](http://pngu.mgh.harvard.edu/~purcell/plink/dosage.shtml) for more details ) to allele dosage format required by Predixcan

2. Predicts the gene expression using Predixcan [software](https://github.com/hakyim/PrediXcan).

3. Runs an association analysis between the predicted expression data and a phenotype.

4. Creates a manhattan plot of the results obtained from the step 3.

## Requirements

1. Python ≥ 2.7 and numpy package.
    
2. R and packages(data.table, doParallel, foreach, qqman)
Note: the script will install the packages if not installed already 

3. Linux environment

## Installation

Clone the repository and run the install.sh script once. 

```
git clone https://github.com/veera-dr/Predixcan.ipsych.git 
cd predixcan.ipsych
sh install.sh
```

The install.sh script appends an alias to your bashrc file, so you can call the script by `predixcan`

You can test if installed by

```
predixcan -h 
```
## Usage instructions

The help info from the script is detailed and you can jump start from there quickly.

# Extraction 

Usage:

```
predixcan --extract --chunkdir <value> --outdir <value> 
```

| Argument    | Value                                                      |
|-------------|------------------------------------------------------------|
| `--outdir`  | full path to the folder where you need the outputs         |
| `--chunkdir`| full path to the 'qc1' folder containing the dosage chunks |

*Note: the script by default will look for the 'info' folder in the same directory where 'qc1' folder is located. This is for quality filtering based on infoscore thresholds. You can skip this filtering step with `--noinfo` argument*

#### Details
This will extract a subset of hapmap variants (predixcan requires only subset of variants for prediction) from the dosage files and will convert them to allele dosage format required by predixcan software. In addition, the script also filters out the variants with info score < 0.60 (you can change this with --infoscore argument) before conversion. The converted dosage files will be placed in a folder named "dosages" inside the outdir you specify.


# Prediction

Usage:

```
predixcan --predict --outdir <value> --weights <value> 
```
| Argument    | Value                                                          |
|-------------|----------------------------------------------------------------|
| `--outdir`  | full path to the folder where you need the outputs             |
| `--weights` | full path to the 'weights' folder containing the '.db' files   |

Note: Be default the script looks for folder "dosages" inside the output folder (outdir) that you specify. You can override this with argument `--dosagefolder`. Inside dosage folder there should be a file with name 'samples.txt', containing FIDs and IIDs in the same order as in the dosage files. One of the fam files inside qc1 folder is copied to this folder during the extraction step.

### Details
The script will use the Predixcan.py ([link](https://github.com/hakyim/PrediXcan)) script to predict the gene expressions using the genotype data. The `--weights` argument requires a _path to a folder (not a file)_  containing pre calculated weights for various tissues. The files have `.db` extension and more details can be found in Predixcan page ([link](https://github.com/hakyim/PrediXcan)). 

# Association

Usage:

```
predixcan --assoc --outdir <value> --pheno <value> 
```

| Argument    | Value                                                          |
|-------------|----------------------------------------------------------------|
| `--outdir`  | full path to the folder where you need the outputs             |
| `--pheno`   | full path to the pheno file (same as plink format)             |

Note: the script by default looks for a folder 'expression' inside the output folder(outdir) that you specify. The expression folder should contain the predicted expressions (one to many), which are the output of the prediction step. You can override this with `--expression` argument pointing to a folder with predicted expressions located elsewhere.

### Details
The script will run an association analysis using the PrediXcanAssociation.R Rscript. The Rscript is modified from the actual version. By default the Rscript will run a logistic regression analysis with the phenotype. The phenotype should be coded as 0 and 1 (*not 1 and 2*). In case, if you are using a continuous parameter as a phenotype, just add `--linear` tag in the command.

# Manhattan plot

Usage:

```
predixcan --manhattan --outdir <value>
```
| Argument    | Value                                                          |
|-------------|----------------------------------------------------------------|
| `--outdir`  | full path to the folder where you need the outputs             |

Note: The script by default looks for a folder 'results' inside the output folder (outdir) that you specify. You can override this using argument `--results` with the path to a _folder (not a file)_ containing the results.
