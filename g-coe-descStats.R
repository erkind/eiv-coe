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
descstats.R.mat <- data.frame(
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
for(i in 1:length(colnames(R.mat))){
   # labels
   industry= colnames(R.mat)[i]
   # descriptive stats
   nb.obs  = as.numeric(length(R.mat[,i]))
   min.val = min(R[,i])
   min.val = as.numeric(quantile(R.mat[,i])[1])
   Q1.val  = as.numeric(quantile(R.mat[,i])[2])
   Q2.val  = as.numeric(quantile(R.mat[,i])[3])
   Q3.val  = as.numeric(quantile(R.mat[,i])[4])
   max.val = as.numeric(quantile(R.mat[,i])[5])
   avg.val = mean(R.mat[,i])
   st.dev  = sd(R.mat[,i])
   sk      = skewness(R.mat[,i])
   t.sk    = sk / sqrt(6 / nb.obs)
   kurt    = kurtosis(R.mat[,i])
   t.kurt  = kurt / sqrt(24 / nb.obs)
   ks.stat = as.numeric(
      ks.test(R.mat[,i], "pnorm", mean(R.mat[,i]), sd = sd(R.mat[,i]),
              alternative = c("two.sided"))$statistic)
   ks.pval = as.numeric(
      ks.test(R.mat[,i], "pnorm", mean(R.mat[,i]), sd = sd(R.mat[,i]),
              alternative = c("two.sided"))$p.value)
   jb.stat = as.numeric(jarqueberaTest(R.mat[,i])@test$statistic)
   jb.pval = as.numeric(jarqueberaTest(R.mat[,i])@test$p.value)
   sw.stat = as.numeric(normalTest(R.mat[,i], method = c("sw"))@test$statistic)
   sw.pval = as.numeric(normalTest(R.mat[,i], method = c("sw"))@test$p.value)
   # append results
   descstats.R.mat <- rbind(
      descstats.R.mat,
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
descstats.F.mat <- data.frame(
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
for(i in 1:length(colnames(F.mat))){
   # labels
   factor  = colnames(F.mat)[i]
   # descriptive stats
   nb.obs  = as.numeric(length(F.mat[,i]))
   min.val = min(F.mat[,i])
   min.val = as.numeric(quantile(F.mat[,i])[1])
   Q1.val  = as.numeric(quantile(F.mat[,i])[2])
   Q2.val  = as.numeric(quantile(F.mat[,i])[3])
   Q3.val  = as.numeric(quantile(F.mat[,i])[4])
   max.val = as.numeric(quantile(F.mat[,i])[5])
   avg.val = mean(F.mat[,i])
   st.dev  = sd(F.mat[,i])
   sk      = skewness(F.mat[,i])
   t.sk    = sk / sqrt(6 / nb.obs)
   kurt    = kurtosis(F.mat[,i])
   t.kurt  = kurt / sqrt(24 / nb.obs)
   ks.stat = as.numeric(
      ks.test(F.mat[,i], "pnorm", mean(F.mat[,i]), sd = sd(F.mat[,i]),
              alternative = c("two.sided"))$statistic)
   ks.pval = as.numeric(
      ks.test(F.mat[,i], "pnorm", mean(F.mat[,i]), sd = sd(F.mat[,i]),
              alternative = c("two.sided"))$p.value)
   jb.stat = as.numeric(jarqueberaTest(F.mat[,i])@test$statistic)
   jb.pval = as.numeric(jarqueberaTest(F.mat[,i])@test$p.value)
   sw.stat = as.numeric(normalTest(F.mat[,i], method = c("sw"))@test$statistic)
   sw.pval = as.numeric(normalTest(F.mat[,i], method = c("sw"))@test$p.value)
   # append results
   descstats.F.mat <- rbind(
      descstats.F.mat,
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
View(descstats.R.mat)
View(descstats.F.mat)