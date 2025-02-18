# Introduction

SumHEM is a heteroscedastic effects model designed to estimate genome-wide SNP effects using GWAS summary statistics. It is a summary-level version of bigRR (https://rdrr.io/rforge/bigRR/man/bigRR-package.html). Using a fast generalized ridge regression (GRR) algorithm under the assumption of uneven variance, SumHEM effectively differentiates the contribution of more or less important variants to phenotypic variation. 

Simulation studies show that SumHEM outperforms LDpred2 (https://github.com/privefl/bigsnpr) in terms of out-of-sample prediction, heritability inference, and estimation consistency. Real data validation further supports these findings, highlighting the effectiveness of SumHEM, especially for highly polygenic traits. 

# System requirements
SumHEM is a easy function to use and requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
This package is supported for macOS and Linux.

## R requirements 
R dependencies: 
```
# Install packages in R
install.packages("dplyr")
install.packages("Matrix")
install.packages("parallel")
install.packages("bigsparser")
install.packages("bigsnpr")
```

# Installation

Install the SumHEM package in R: 
```
# install.packages("remotes")
remotes::install_github("yueyao0/SumHEM")
library(SumHEM)
```

Alternatively, clone the repository and load the SumHEM function: 
```
git clone https://github.com/yueyao0/SumHEM.git
source("SumHEM/R/SumHEM.R")
```

# Quick vignette

You will need the GWAS summary statistics and heritability (optional) of a trait, as well as a linkage disequilibrium (LD) reference. The LD files for the European ancestry population from the UK Biobank, authored by Florian Privé, can be downloaded here: https://figshare.com/articles/dataset/European_LD_reference/13034123. Once downloaded, place the LD `rds` files in a folder and name each file as `LD_chr[chromosome number].rds`. 

For a detailed, step-by-step tutorial on how to conduct your analysis using SumHEM, please refer to our comprehensive guide available below. 

```
# Load example GWAS summary statistics for standing height
df_gwas <- readRDS("SumHEM/data/gwas.height.rds")
# Load SNP information of LD reference panel
df_map <- readRDS("SumHEM/data/map.rds")

# Setting for parameters
ws <- 500 # the number of SNPs in a window for estimation
NCORES <- 10 # the number of cores using in parallel computation
LD_path <- "SumHEM/LD_EUR" # path of folder with LD files
h2input <- 0.4780 # heritability estimated by an external software (optional)

# Estimate SNP effects via SumHEM
df_est <- SumHEM(
  df_gwas = df_gwas, df_map = df_map,
  h2input = h2input, # if missing, it will be computed during the process
  ws = ws, NCORES = NCORES, LD_path = LD_path
)
```

During the process of SumHEM estimation, there will be some information displayed on the screen as a log. 

```
# Analysis started on Wed Feb 19 03:37:24 2025
# Estimating LDSC heritability ...        Time 0.0779950896898905 min 
# LDSC estimates of heritability = 0.47804031
# The number of matched SNPs = 1050026
# Estimating genetic effects ...
# Estimating in chromosome 1 ...  Time 0.32837632894516 min 
# Estimating in chromosome 2 ...  Time 0.32187077999115 min 
# Estimating in chromosome 3 ...  Time 0.284389293193817 min 
# Estimating in chromosome 4 ...  Time 0.230345896879832 min 
# Estimating in chromosome 5 ...  Time 0.24179284175237 min 
# Estimating in chromosome 6 ...  Time 0.313477178414663 min 
# Estimating in chromosome 7 ...  Time 0.210832905769348 min 
# Estimating in chromosome 8 ...  Time 0.204874690373739 min 
# Estimating in chromosome 9 ...  Time 0.168351852893829 min 
# Estimating in chromosome 10 ... Time 0.211142031351725 min 
# Estimating in chromosome 11 ... Time 0.200521258513133 min 
# Estimating in chromosome 12 ... Time 0.186887427171071 min 
# Estimating in chromosome 13 ... Time 0.140968612829844 min 
# Estimating in chromosome 14 ... Time 0.124034980932872 min 
# Estimating in chromosome 15 ... Time 0.108550055821737 min 
# Estimating in chromosome 16 ... Time 0.109349226951599 min 
# Estimating in chromosome 17 ... Time 0.10762494802475 min 
# Estimating in chromosome 18 ... Time 0.108406146367391 min 
# Estimating in chromosome 19 ... Time 0.0756011247634888 min 
# Estimating in chromosome 20 ... Time 0.099292512734731 min 
# Estimating in chromosome 21 ... Time 0.0587003747622172 min 
# Estimating in chromosome 22 ... Time 0.0610900640487671 min 
# SumHEM estimation is DONE with Time 3.89806080261866 min 
# Analysis finished on Wed Feb 19 03:41:24 2025
```

Here is an overview of the results from SumHEM function. 

```
# View the results
str(df_est)
# 'data.frame':   1050026 obs. of  24 variables:
#  $ rsid         : chr  "rs3131969" "rs12562034" "rs4040617" "rs4970383" ...
#  $ chr          : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ pos          : int  754182 768448 779322 838555 846808 853954 854250 864938 870645 873558 ...
#  $ beta         : num  0.009541 -0.001564 -0.009661 -0.00464 0.000286 ...
#  $ se           : num  0.00466 0.00508 0.00467 0.00363 0.00393 ...
#  $ N            : int  193785 193785 193785 193785 193785 193785 193785 193785 193785 193785 ...
#  $ scale        : num  2.05 2.24 2.05 1.6 1.73 ...
#  $ sbeta        : num  0.004649 -0.000699 -0.004703 -0.002903 0.000166 ...
#  $ h2input      : num  0.478 0.478 0.478 0.478 0.478 ...
#  $ M            : int  1050026 1050026 1050026 1050026 1050026 1050026 1050026 1050026 1050026 1050026 ...
#  $ winSize      : int  500 500 500 500 500 500 500 500 500 500 ...
#  $ winID        : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ sig2e        : num  1 1 1 1 1 ...
#  $ sig2b_SumRR  : num  4.55e-07 4.55e-07 4.55e-07 4.55e-07 4.55e-07 ...
#  $ lambda_SumRR : num  2196022 2196022 2196022 2196022 2196022 ...
#  $ b_SumRR      : num  3.58e-04 -7.08e-05 -3.61e-04 -2.67e-04 -3.16e-05 ...
#  $ Rb_SumRR     : num  0.000595 0.000103 -0.000609 0.000122 0.000524 ...
#  $ h2_SumRR     : num  0.000242 0.000242 0.000242 0.000242 0.000242 ...
#  $ hv           : num  0.926 0.92 0.926 0.929 0.937 ...
#  $ sig2b_SumHEM : num  1.72e-06 6.29e-08 1.75e-06 1.01e-06 1.58e-08 ...
#  $ lambda_SumHEM: num  581443 15882144 570044 990728 63434747 ...
#  $ b_SumHEM     : num  9.63e-04 -1.17e-05 -9.92e-04 -5.49e-04 -2.02e-06 ...
#  $ Rb_SumHEM    : num  1.76e-03 2.56e-04 -1.78e-03 -9.59e-05 8.28e-04 ...
#  $ h2_SumHEM    : num  0.000526 0.000526 0.000526 0.000526 0.000526 ...
#  $ h2var_SumHEM : num  8.38e-07 8.38e-07 8.38e-07 8.38e-07 8.38e-07 ...

# Summaries the inferred heritability
df <- df_est %>% distinct(chr,winID, .keep_all = TRUE)
sum(df$h2_SumRR) # inferred heritability of SumRR
# [1] 0.2792048
sum(df$h2_SumHEM) # inferred heritability of SumHEM
# [1] 0.6565036
sqrt(sum(df$h2var_SumHEM)) # se of SumHEM-inferred heritability
# [1] 0.03699596
```

# Input format
The essential input files of SumHEM include the LD reference and the GWAS summary statistics of a trait. 

## LD reference

The genome-wide LD matrix should be split into separate `rds` files for each chromosome. These files should be named as `LD_chr[chromosome number].rds` and placed in the same folder. 

The SNP information in the LD reference should include the following columns with the same column names. If you do not provide the trait heritability for the analysis, please include column `ld` in the data frame, and SumHEM will estimate the value during the process. If you provide the heritability, the column `ld` is optional. 

| Column name | Description | Class in R |
| :--- | :--- | :--- |
| rsid | SNP ID | character |
| chr | Chromosome number | integer |
| pos| Position of SNP | integer |
| ld (optional) | LD score for LDSC heritability model | numeric |

We provide the genome-wide SNP information of LD reference (downloaded from https://figshare.com/articles/dataset/European_LD_reference/13034123) for 1,054,330 HapMap3 variants based on 362,320 European individuals of the UK Biobank. 

```
str(df_map)
# 'data.frame':   1054330 obs. of  4 variables:
#  $ rsid: chr  "rs3131972" "rs3131969" "rs1048488" "rs12562034" ...
#  $ chr : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ pos : int  752721 754182 760912 768448 779322 838555 846808 853954 854250 864938 ...
#  $ ld  : num  3.69 3.73 3.69 1.4 3.68 ...
```

## GWAS summary statistics

It is crucial to format the summary statistics correctly before your analysis with SumHEM. Your input GWAS data file should include the columns as below and with the same column names. Please make sure that the effect allels in GWAS data are matched with those in the reference. 

| Column name | Description | Class in R |
| :--- | :--- | :--- |
| rsid | SNP ID | character |
| chr | Chromosome number | integer |
| pos | Position of SNP | integer |
| beta | Estimated marginal effect | numeric |
| se | Standard error of the estimated marginal effect | numeric |
| N | Sample size | integer |

We provide an example GWAS summary statistics of standing height from the UK Biobank. 

```
str(df_gwas)
# 'data.frame':   1050026 obs. of  6 variables:
#  $ rsid: chr  "rs3131969" "rs12562034" "rs4040617" "rs4970383" ...
#  $ chr : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ pos : int  754182 768448 779322 838555 846808 853954 854250 864938 870645 873558 ...
#  $ beta: num  0.009541 -0.001564 -0.009661 -0.00464 0.000286 ...
#  $ se  : num  0.00466 0.00508 0.00467 0.00363 0.00393 ...
#  $ N   : int  193785 193785 193785 193785 193785 193785 193785 193785 193785 193785 ...
```

# Parameters of SumHEM function

Here are some explanation of parameters in SumHEM function. 

| Parameter | Default value | Description | Class in R |
| :--- | :--- | :--- | :--- |
| df_gwas | "" | A data frame of GWAS summary statistic of trait | data.frame |
| df_map | "" | A data frame of SNP information of LD refrence | data.frame |
| h2input (optional) | "" | Heritability of trait | numeric |
| ws | 500 | The number of SNPs in a window for estimation | integer |
| NCORES | 1 | The number of cores using in parallel computation | integer |
| LD_path | "" | A path of folder with LD files | character |

# Variables in SumHEM output

After analysis, the SumHEM function will return a data frame of estimation results. Here is the explanation of each column. 

| Column name | Description |
| :--- | :--- |
| rsid | SNP ID |
| chr | Chromosome number |
| pos | Position of SNP |
| beta | Estimated marginal effect in GWAS |
| se | Standard error of the estimated marginal effect in GWAS |
| N | Sample size in GWAS |
| scale | sqrt(beta^2+N*se^2), the scale used for normalization |
| sbeta | beta/scale |
| h2input | Input value of parameter h2input in function |
| M | The number of SNP in whole-genome |
| winSize | The number of SNP in the estimation window including the SNP |
| winID | ID of the estimation window including the SNP |
| sig2e | The value of (sigma_e)^2 in the estimation window |
| sig2b_SumRR | The value of (sigma_b)^2 for SumRR in the estimation window |
| lambda_SumRR | The value of shrinkage parameter lambda for SumRR in the estimation window |
| b_SumRR | SumRR estimates for joint effect of the SNP |
| Rb_SumRR | LD %*% b_SumRR in the estimation window |
| h2_SumRR | The inferred heritability of SumRR in the estimation window |
| hv | Hat values used in the estimation window |
| sig2b_SumHEM | The value of (sigma_b)^2 for SumHEM in the estimation window |
| lambda_SumHEM | The value of shrinkage parameter lambda for SumHEM in the estimation window |
| b_SumHEM | SumHEM estimates for joint effect of the SNP |
| Rb_SumHEM | LD %*% b_SumHEM in the estimation window |
| h2_SumHEM | The inferred heritability of SumHEM in the estimation window |
| h2var_SumHEM | The variance of SumHEM-inferred heritability in the estimation window |

# Citation
If you use the SumHEM software, please cite: 

XXX

Xia Shen, Moudud Alam, Freddy Fikse, Lars Rönnegård, A Novel Generalized Ridge Regression Method for Quantitative Genetics, Genetics, Volume 193, Issue 4, 1 April 2013, Pages 1255–1268, https://doi.org/10.1534/genetics.112.146720

# License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

# Contact
Report bugs by opening a new issue on this GitHub page.

Send email to the authors: Yue Yao or Xia Shen.

# Future plans
This repository will be integrated into SumHEM software in the future.






