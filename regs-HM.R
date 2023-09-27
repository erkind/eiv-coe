# HM estimation of factor models for ind portfolios

# data ------------------------------------------------------------------------
# # loading from hard drive, mind customizing the path
# source(paste("C:/Users/68596/Dropbox/eiv-coe/",
#              "scripts/data-eiv-coe.R",
#              sep = "")
#        )
# loading from GitHub
source(paste("http://raw.githubusercontent.com/erkind/",
             "eiv-coe/main/data-eiv-coe.R",
             sep = "")
       )
dim(data); colnames(data)[c(1:5,60:64)]
dim(R); colnames(R)
dim(F1); colnames(F1)
dim(F3); colnames(F3)
dim(F5); colnames(F5)
