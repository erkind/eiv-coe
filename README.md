TITLE
Higher-Moment Estimation of Industry-level Cost of Equity Capital with Errors-in-variables

DESCRIPTION
Public repository for reproducible results of the paper

SYSTEM
All scripts have been run on an RStudio 2022.02.3+492 Prairie Trillium Release using R version 4.2.1 installed on a Windows 10 Intel(R) Core(TM) processor i5-4310U CPU 2.00GHz 2.60GHz 64bits x64 system

REQUIRED PACKAGES
car
dplyr
fbasics
lmtest
sandwich
skedastic

USAGE
g-coe-data-eiv.R
The file loads the original dataset and sets the subsample employed in the paper. It declares "R" the matrix of industry portfolio excess returns over the risk-free rate and "X" the matrix of systematic risk factors.





