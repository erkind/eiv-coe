# TITLE AND DESCRIPTION

Public repository for reproducible results of the paper 'Higher-Moment Estimation of Industry-level Cost of Equity Capital with Errors-in-variables'

# SYSTEM

All scripts have been run on an RStudio 2022.02.3+492 Prairie Trillium Release using R version 4.2.1 installed on a Windows 10 Intel(R) Core(TM) processor i5-4310U CPU 2.00GHz 2.60GHz 64bits x64 system.

# REQUIRED PACKAGES

-   car
-   dplyr
-   fBasics
-   lmtest
-   sandwich
-   skedastic

# USAGE

Below we describe each file we employed in the paper. Some of them are optional and the user needs not to run the codes within to obtain the results.

## data-eiv-coe.csv (optional)

This is the csv data file used in the paper. The user does not need to load this file, it is automatically loaded in the next script. Raw data on industry portfolios as well as systematic risk factors have been obtained from Kenneth F. French [online library](http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html).

## g-coe-data-eiv.R (optional)

The script is optional to the user.

Its function is to load the raw data uploaded on the same repository as a csv file, then defines a subsample that is used in the paper. It assigns the variable 'R.mat' to represent the matrix of industry portfolio excess returns over the risk-free rate, and 'F.mat' to denote the matrix of systematic risk factors.

## g-f-HM-Zmat.R (optional)

This script is optional to the user.

This function, generate.Z.mat, calculates the Z-matrix of higher-moments instruments based on provided input variables. This function takes parameters such as the matrix of independent variables (X.mat), an optional vector of dependent variables (Y.vec), sample size (n), number of predictors (k), and a list of higher-moments (iv.list). It constructs the Z-matrix by computing various higher-moment instruments, such as quadratic and cubic terms, and combines them into the Z-matrix, which is then returned.

## g-coe-descStats.R

The user can run this script in a new R session.

The script begins by importing the data by loading the previous script 'g-coe-data-eiv.R'. It computes various summary statistics and performs normality tests on both dependent and independent variables. The resulting outcomes are stored in two distinct data frames, which can be exported as CSV files using commonly available R methods.

## g-coe-regs.R

The user can run this script in a new R session.

The script implements the OLS and Higher-moments estimations of the regressions of industry portfolio excess returns on Fama-French risk factors. We consider for specification where 'M1' refers to the CAPM, 'M3' to Fama-French three-factor model, 'M4' to Fama-French and Carhart's four-factor model, and 'M5' to Fama-French five-factor model. For each specification, OLS and HM moments regressions are run using a 'for' loop and the output is collected in a separate data frame. HM estimations require the computation of the residuals of the regression of risk factors on the Z matrix of instruments. This is performed in the sections 'Z matrix' and 'regs F on Z'. Selected instruments include z0, z1 and z4 (technical details and discussion provided in the paper). HM regression output also includes the Durbin-Wu-Hausman test results to test the null of no errors-in-variables (EIV) as well as the test for overidentification restrictions. The section 'weak instruments' controls for the problem of weak instruments and reports the results for the F-statistic from the regression of each risk factor on instrumental variables following Stock and Yogo (2005) and Olea and Pflueger (2013). Finally, the section 'instrument exogeneity' controls for the presence of exogenous instruments. This is performed by checking the goodness-of-the-fit statistics of the regressions of the fitted HM regression residuals on each instrument set.

# Contact

For any questions or feedback, feel free to contact us at [erkin.d\@u-pec.fr](mailto:erkin.d@u-pec.fr){.email} or [marc.desban\@u-pec.fr](mailto:marc.desban@u-pec.fr){.email}.
