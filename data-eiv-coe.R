# eiv-coe data
#
# # loading from hard drive, mind customizing the path
# data <- read.csv(
#    paste(
#       "c:/users/68596/",
#       "Dropbox/eiv-coe/data/",
#       "data-eiv-coe.csv",
#       sep = ""),
#    header = TRUE)
# loading from GitHub
data <- read.csv(
   paste(
      "http://raw.githubusercontent.com/erkind/",
      "eiv-coe/main/data-eiv-coe.csv",
      sep = ""),
   header = TRUE)
dim(data)
colnames(data)
head(data[c(1:3,51:53,57:59)], n = 2)
tail(data[c(1:3,51:53,57:59)], n = 2)
#
# set sample period
t.start<- as.numeric(which(data$date == 199011))
t.end  <- as.numeric(which(data$date == 202202))
#
# excess returns on industry portfolios
R <- as.matrix(data[c(t.start:t.end),2:45], # exclude financials
               nrow = t.end - t.start + 1,
               ncol = 44)
R <- R - data[c(t.start:t.end), 56] # subtract RF
dim(R); colnames(R)
#
# factors
F1 <- as.matrix(data[c(t.start:t.end),51],
                nrow = t.end - t.start + 1,
                ncol = 1)
colnames(F1) <- colnames(data)[51]
F3 <- as.matrix(data[c(t.start:t.end),51:53],
                nrow = t.end - t.start + 1,
                ncol = 3)
F4 <- as.matrix(data[c(t.start:t.end),c(51:53, 61)],
                nrow = t.end - t.start + 1,
                ncol = 4)
F5 <- as.matrix(data[c(t.start:t.end),51:55],
                nrow = t.end - t.start + 1,
                ncol = 5)
#
# counters
nb.obs <- as.numeric(dim(R)[1])
nb.port<- as.numeric(dim(R)[2])
