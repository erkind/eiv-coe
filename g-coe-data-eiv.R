# eiv-coe data
#
# loading raw data
data <- read.csv(paste("http://raw.githubusercontent.com/erkind/",
                       "eiv-coe/main/data-eiv-coe.csv",
                       sep = ""),
                 header = TRUE, sep = ",")
#
# checking data
dim(data)
colnames(data)
head(data[c(1:3,51:53,57:59)], n = 2)
tail(data[c(1:3,51:53,57:59)], n = 2)
#
# set sample period to 199011-202202
t.start<- as.numeric(which(data$date == 199011))
t.end  <- as.numeric(which(data$date == 202202))
#
# vector of risk-free rate
RF <- as.matrix(data[c(t.start:t.end), 56],
                nrow = t.end - t.start + 1,
                ncol = 1)
colnames(RF) <- c("RF")
dim(RF)
colnames(RF)
#
# matrix of industry portfolio excess returns
R.mat <- as.matrix(data[c(t.start:t.end),2:45],
                   nrow = t.end - t.start + 1,
                   ncol = 44)
R.mat <- R.mat - data[c(t.start:t.end), 56]
dim(R.mat)
colnames(R.mat)
#
# matrix of risk factors
F.mat <- as.matrix(data[c(t.start:t.end), c(51:55,61)],
                   nrow = t.end - t.start + 1,
                   ncol = 6)
dim(F.mat)
colnames(F.mat)
#
# remove unnecessary objects
rm("t.start", "t.end")