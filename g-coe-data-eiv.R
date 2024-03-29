# eiv-coe data
#
# loading raw csv file
data <- read.csv(paste("C:/Users/68596/Dropbox/",
                       "eiv-coe/data/data-eiv-coe.csv", sep = ""),
                 header = TRUE, sep = ",")
# data <- read.csv(paste("http://raw.githubusercontent.com/erkind/",
#                        "eiv-coe/main/data-eiv-coe.csv",
#                        sep = ""),
#                  header = TRUE, sep = ",")
dim(data)
colnames(data)
head(data[c(1:3,51:53,57:59)], n = 2)
tail(data[c(1:3,51:53,57:59)], n = 2)
#
# set sample period
t.start<- as.numeric(which(data$date == 199011))
t.end  <- as.numeric(which(data$date == 202202))
#
# risk-free rate
RF <- as.matrix(data[c(t.start:t.end), 56],
                nrow = t.end - t.start + 1,
                ncol = 1)
colnames(RF) <- c("RF")
dim(RF); colnames(RF)
#
# excess returns on industry portfolios
R <- as.matrix(data[c(t.start:t.end),2:45], # exclude financials
               nrow = t.end - t.start + 1,
               ncol = 44)
R <- R - data[c(t.start:t.end), 56] # subtract RF
dim(R); colnames(R)
#
# factors
X <- as.matrix(data[c(t.start:t.end), c(51:55,61)],
               nrow = t.end - t.start + 1,
               ncol = 6)
dim(X); colnames(X)
#
rm("t.start", "t.end")