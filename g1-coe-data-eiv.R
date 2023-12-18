# eiv-coe data
#
# loading raw csv file
data <- read.csv(paste("http://raw.githubusercontent.com/erkind/",
                       "eiv-coe/main/data-eiv-coe.csv",
                       sep = ""),
                 header = TRUE, sep = ",")
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
#
# excess returns on industry portfolios
R <- as.matrix(data[c(t.start:t.end),2:45], # exclude financials
               nrow = t.end - t.start + 1,
               ncol = 44)
R <- R - data[c(t.start:t.end), 56] # subtract RF
#
# factors
X <- as.matrix(data[c(t.start:t.end), c(51:55,61)],
               nrow = t.end - t.start + 1,
               ncol = 6)
#
rm("t.start", "t.end")