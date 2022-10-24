# MissHybrid, Imputation of Missing data with Hybrid Approach

## Introduction

Missing data is present in broad area of scientific researches. Imputation of missing data is a challenging task in many biomedical studies although much needed. Both parametric and non-parametric imputation methods have been previously proposed to solve imputation problem. Here we combined some of the popular imputation methods into a two-phase hybrid imputation framework. Through evaluation of imputing performance on both simulated and real data, we showed the hybrid imputation framework has performance advantage over most the popular methods. Specifically, we showed that our imputation framework works best for imputing high-dimensional genomics data such as genome-wide gene expression.

## Key features
1. combine the advantages of both parametric and non-parametric imputation method to fully capture the complex structures among variables in each dataset
2. Improved imputation performance in high dimensional setting
3. Integrate popular methods such as KNN, softImpute and missMDA into the hybrid framework and implemented as an easy-to-use R package

## Install GitHub Version
To install `missHybrid` directly from GitHub, run

```r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("argossy/missHybrid@main")








```

The package includes reference manual, sample data and a Vignette.

## Basic Usage

```r
library(missHybrid)

## simulate data
n = 100
p = 10
dat = matrix(rnorm(1000), 100,10)

## 10% missing
mis_rate = 0.1
for(j in 1:p){
  idx_mis = sample(1:n, n * mis_rate, replace = FALSE)
  dat[idx_mis, j ] = NA
}

# examples of two-phase imputation

# Default method, phase 1: missMDA; phase 2: missForest
Ximp = missHyb(X = dat)

# Specify two phase method, phase 1: missMDA; phase 2: missForest
Ximp = missHyb(X = dat, p1 = 'MDA', p2 = 'MF') # or missHyb(dat, 'MDA', 'MF')

# Specify two phase method, phase 1: missForest; phase 2: MDA
Ximp = missHyb(X = dat, p1 = 'MF', p2 = 'MDA')

# Specify two phase method, phase 1: KNN; phase 2: missMDA
Ximp = missHyb(X = dat, p1 = 'KNN', p2 = 'MDA')




```

## Citation
K Xia, B Zhao, J Gilmore, F Zou. MissHybrid, Imputation of Missing data with Hybrid Approach. (In Preparation)

## Contact
Kai Xia: kxia@med.unc.edu

