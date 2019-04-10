
### Box-Cox transformations

## analyses for McDaniel et al done with R version 3.5.3
## MASS version 7.3-51.4


install.packages("MASS")
library("MASS")

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


################################
# sexes combined
BCdtg <- boxcox(DtG ~ Population, 
                plotit = TRUE, data = combined)
with(BCdtg, x[which.max(y)])
# result: -1.6
# will use: -1.75

BCdt10g <- boxcox(Dt10G ~ Population, 
                plotit = TRUE, data = combined)
with(BCdt10g, x[which.max(y)])
# result: -0.1818182
# will use: -0.25

BCps <- boxcox(PS ~ Population, 
                plotit = TRUE, data = combined)
with(BCps, x[which.max(y)])
# result: 0.7474747
# will use: 0.75

BCll <- boxcox(LL ~ Population, 
                plotit = TRUE, data = combined)
with(BCll, x[which.max(y)])
# result: 0.22222222
# will use: 0.25

#########################################
# males

BCdtgM <- boxcox(DtG ~ Population, 
                plotit = TRUE, data = maledata)
with(BCdtgM, x[which.max(y)])
# result: -2
# will use: -2

BCdt10gM <- boxcox(Dt10G ~ Population, 
                plotit = TRUE, data = maledata)
with(BCdt10gM, x[which.max(y)])
# result: -0.3838384
# will use: -0.5

BCpsM <- boxcox(PS ~ Population, 
                plotit = TRUE, data = maledata)
with(BCpsM, x[which.max(y)])
# result: 1.030303
# will use: 1

BCllM <- boxcox(LL ~ Population, 
                plotit = TRUE, data = maledata)
with(BCllM, x[which.max(y)])
# result: 0.1010101
# will use: 0

###make the 0s become NA for NP just for the boxcox test
maledata$NP[maledata$NP==0] <- NA
BCnp <- boxcox(NP ~ Population, 
               plotit = TRUE, data = maledata)
with(BCnp, x[which.max(y)])
# result: 0.5858586
# will use: 0.5

BCdtp <- boxcox(DtP ~ Population, 
                plotit = TRUE, data = maledata)
with(BCdtp, x[which.max(y)])
# result: 2
# will use: 2

#########################################
# females

BCdtgF <- boxcox(DtG ~ Population, 
                plotit = TRUE, data = femaledata)
with(BCdtgF, x[which.max(y)])
# result: -1.272727
# will use: -1.25

BCdt10gF <- boxcox(Dt10G ~ Population, 
                plotit = TRUE, data = femaledata)
with(BCdt10gF, x[which.max(y)])
# result: -0.22222222
# will use: -0.25

BCpsF <- boxcox(PS ~ Population, 
                plotit = TRUE, data = femaledata)
with(BCpsF, x[which.max(y)])

# result: 0.707070
# will use: 0.75

BCllF <- boxcox(LL ~ Population, 
                plotit = TRUE, data = femaledata)
with(BCllF, x[which.max(y)])
# result: 0.2222222
# will use: 0.25

################ transformations to use on the data based on boxcox test ###################

## males
maledatatransDtG <- ((maledata$DtG)^(-2))
maledatatransDt10G <- ((maledata$Dt10G)^(-0.5)) 
maledatatransPS <- (maledata$PS)
maledatatransLL <- log(maledata$LL)
maledatatransDtP <- ((maledata$DtP)^(2))
maledatatransNP <- ((maledata$NP)^(0.5))


## females
femaledatatransDtG <- ((femaledata$DtG)^(-1.25))
femaledatatransDt10G <- ((femaledata$Dt10G)^(-0.25))
femaledatatransPS <- ((femaledata$PS)^(0.75))
femaledatatransLL <- ((femaledata$LL)^(0.25))

## combined
combinedtransDtG <- ((combined$DtG)^(-1.75))
combinedtransDt10G <- ((combined$Dt10G)^(-0.25))
combinedtransPS <- ((combined$PS)^(0.75))
combinedtransLL <- ((combined$LL)^(0.25))

