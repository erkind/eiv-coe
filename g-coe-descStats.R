# descriptive statistics
#
# call data file
source(paste("http://raw.githubusercontent.com/erkind/",
             "eiv-coe/main/g-coe-data-eiv.R",
             sep = ""))
#
# required
library(fBasics)
#
# summary stats for ind portfolios ex returns
descstats.R <- data.frame(
   series = character(0),
   nb.obs = numeric(0),
   min.val = numeric(0),
   Q1.val = numeric(0), Q2.val = numeric(0), Q3.val = numeric(0),
   max.val = numeric(0),
   avg = numeric(0),
   st.dev = numeric(0),
   sk = numeric(0), t.sk = numeric(0),
   kurt = numeric(0), t.kurt = numeric(0),
   ks.stat = numeric(0), ks.pval = numeric(0),
   jb.stat = numeric(0), jb.pval = numeric(0),
   sw.stat = numeric(0), sw.pval = numeric(0)
   )
for(i in 1:length(colnames(R))){
   # labels
   industry= colnames(R)[i]
   # descriptive stats
   nb.obs  = as.numeric(length(R[,i]))
   min.val = min(R[,i])
   min.val = as.numeric(quantile(R[,i])[1])
   Q1.val  = as.numeric(quantile(R[,i])[2])
   Q2.val  = as.numeric(quantile(R[,i])[3])
   Q3.val  = as.numeric(quantile(R[,i])[4])
   max.val = as.numeric(quantile(R[,i])[5])
   avg.val = mean(R[,i])
   st.dev  = sd(R[,i])
   sk      = skewness(R[,i])
   t.sk    = sk / sqrt(6 / nb.obs)
   kurt    = kurtosis(R[,i])
   t.kurt  = kurt / sqrt(24 / nb.obs)
   ks.stat = as.numeric(
      ks.test(R[,i], "pnorm", mean(R[,i]), sd = sd(R[,i]),
              alternative = c("two.sided"))$statistic)
   ks.pval = as.numeric(
      ks.test(R[,i], "pnorm", mean(R[,i]), sd = sd(R[,i]),
              alternative = c("two.sided"))$p.value)
   jb.stat = as.numeric(jarqueberaTest(R[,i])@test$statistic)
   jb.pval = as.numeric(jarqueberaTest(R[,i])@test$p.value)
   sw.stat = as.numeric(normalTest(R[,i], method = c("sw"))@test$statistic)
   sw.pval = as.numeric(normalTest(R[,i], method = c("sw"))@test$p.value)
   # append results
   descstats.R <- rbind(
      descstats.R,
      data.frame(
         series  = industry,
         nb.obs  = nb.obs,
         min     = min.val,
         Q1      = Q1.val,
         Q2      = Q2.val,
         Q3      = Q3.val,
         max     = max.val,
         avg     = avg.val,
         std     = st.dev,
         skew    = sk,
         t.skew  = t.sk,
         kurt    = kurt,
         t.kurt  = t.kurt,
         ks.stat = ks.stat,
         ks.pval = ks.pval,
         jb.stat = jb.stat,
         jb.pval = jb.pval,
         sw.stat = sw.stat,
         sw.pval = sw.pval
      )
   )
}
#
# summary stats for risk factors
descstats.X <- data.frame(
   series = character(0),
   nb.obs = numeric(0),
   min.val = numeric(0),
   Q1.val = numeric(0), Q2.val = numeric(0), Q3.val = numeric(0),
   max.val = numeric(0),
   avg = numeric(0),
   st.dev = numeric(0),
   sk = numeric(0), t.sk = numeric(0),
   kurt = numeric(0), t.kurt = numeric(0),
   ks.stat = numeric(0), ks.pval = numeric(0),
   jb.stat = numeric(0), jb.pval = numeric(0),
   sw.stat = numeric(0), sw.pval = numeric(0)
)
for(i in 1:length(colnames(X))){
   # labels
   factor  = colnames(X)[i]
   # descriptive stats
   nb.obs  = as.numeric(length(X[,i]))
   min.val = min(X[,i])
   min.val = as.numeric(quantile(X[,i])[1])
   Q1.val  = as.numeric(quantile(X[,i])[2])
   Q2.val  = as.numeric(quantile(X[,i])[3])
   Q3.val  = as.numeric(quantile(X[,i])[4])
   max.val = as.numeric(quantile(X[,i])[5])
   avg.val = mean(X[,i])
   st.dev  = sd(X[,i])
   sk      = skewness(X[,i])
   t.sk    = sk / sqrt(6 / nb.obs)
   kurt    = kurtosis(X[,i])
   t.kurt  = kurt / sqrt(24 / nb.obs)
   ks.stat = as.numeric(
      ks.test(X[,i], "pnorm", mean(X[,i]), sd = sd(X[,i]),
              alternative = c("two.sided"))$statistic)
   ks.pval = as.numeric(
      ks.test(X[,i], "pnorm", mean(X[,i]), sd = sd(X[,i]),
              alternative = c("two.sided"))$p.value)
   jb.stat = as.numeric(jarqueberaTest(X[,i])@test$statistic)
   jb.pval = as.numeric(jarqueberaTest(X[,i])@test$p.value)
   sw.stat = as.numeric(normalTest(X[,i], method = c("sw"))@test$statistic)
   sw.pval = as.numeric(normalTest(X[,i], method = c("sw"))@test$p.value)
   # append results
   descstats.X <- rbind(
      descstats.X,
      data.frame(
         series  = factor,
         nb.obs  = nb.obs,
         min     = min.val,
         Q1      = Q1.val,
         Q2      = Q2.val,
         Q3      = Q3.val,
         max     = max.val,
         avg     = avg.val,
         std     = st.dev,
         skew    = sk,
         t.skew  = t.sk,
         kurt    = kurt,
         t.kurt  = t.kurt,
         ks.stat = ks.stat,
         ks.pval = ks.pval,
         jb.stat = jb.stat,
         jb.pval = jb.pval,
         sw.stat = sw.stat,
         sw.pval = sw.pval
      )
   )
}
#
# output
View(descstats.R)
View(descstats.X)