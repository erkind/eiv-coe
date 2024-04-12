# TITLE AND DESCRIPTION

Public repository for reproducible results of the paper 'Higher-Moment Estimation of Industry-level Cost of Equity Capital with Errors-in-variables'

# SYSTEM

All scripts have been run on an RStudio 2022.02.3+492 Prairie Trillium Release using R version 4.2.1 installed on a Windows 10 Intel(R) Core(TM) processor i5-4310U CPU 2.00GHz 2.60GHz 64bits x64 system

# REQUIRED PACKAGES

- car
- dplyr
- fbasics
- lmtest
- sandwich
- skedastic

# USAGE

## data-eiv-coe.csv

This is the csv data file used in the paper.

## g-coe-data-eiv.R

The script loads the raw data uploaded on the same repository as a csv file. It defines a subsample as used in the paper. It assigns the variable 'R' to represent the matrix of industry portfolio excess returns over the risk-free rate, and 'X' to denote the matrix of systematic risk factors.

## g-coe-descStats.R

This script begins by importing the 'g-coe-data-eiv.R' file into the session. It then computes various summary statistics and performs normality tests on both dependent and independent variables. The resulting outcomes are stored in two distinct data frames, which can be exported as CSV files using commonly available R methods.

## g-coe-regs.R

Desciption follows...

