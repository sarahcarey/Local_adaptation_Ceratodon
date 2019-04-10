
### Shapiro-Wilks test of normal distribution of data

## analyses for McDaniel et al done with R version 3.5.3


### load in data
combined <- na.omit(read.csv("data/CombinedENA.csv"))
combined$Individual <- as.character(combined$Individual)
combined$Clone <- as.character(combined$Clone)

maledata <- na.omit(read.csv("data/MaleENA.csv"))
maledata$Individual <- as.character(maledata$Individual)
maledata$Clone <- as.character(maledata$Clone)

femaledata <- na.omit(read.csv("data/FemaleENA.csv"))
femaledata$Individual <- as.character(femaledata$Individual)
femaledata$Clone <- as.character(femaledata$Clone)
femaledata$Sex <- as.factor(femaledata$Sex)


# sexes combined
swDtG <- shapiro.test(combined$DtG)
swDtG
#W = 0.81434, p-value < 2.2e-16
swDt10G <- shapiro.test(combined$Dt10G)
swDt10G
#W = 0.88146, p-value < 2.2e-16
swPS <- shapiro.test(combined$PS)
swPS
#W = 0.98274, p-value = 2.262e-09
swLL <- shapiro.test(combined$LL)
swLL
#W = 0.96415, p-value = 8.369e-15


# females
swDtG <- shapiro.test(femaledata$DtG)
swDtG
#W = 0.84412, p-value < 2.2e-16
swDt10G <- shapiro.test(femaledata$Dt10G)
swDt10G
#W = 0.86739, p-value < 2.2e-16
swPS <- shapiro.test(femaledata$PS)
swPS
#W = 0.97905, p-value = 3.685e-07
swLL <- shapiro.test(femaledata$LL)
swLL
#W = 0.97295, p-value = 1.282e-08


# males
swDtG <- shapiro.test(maledata$DtG)
swDtG
#W = 0.7602, p-value < 2.2e-16
swDt10G <- shapiro.test(maledata$Dt10G)
swDt10G
#W = 0.87752, p-value < 2.2e-16
swPS <- shapiro.test(maledata$PS)
swPS
#W = 0.99478, p-value = 0.1676
swLL <- shapiro.test(maledata$LL)
swLL
#W = 0.93714, p-value = 2.645e-12
swNP <- shapiro.test(maledata$NP)
swNP
#W = 0.964, p-value = 1.287e-08
swDtP <- shapiro.test(maledata$DtP)
swDtP
#W = 0.73078, p-value < 2.2e-16



