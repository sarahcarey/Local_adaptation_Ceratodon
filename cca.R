
# Canonical correspondance analysis

## analyses for McDaniel et al done with R version 3.5.3
## vegan version 2.5-4

install.packages("vegan")
library("vegan") 

#### load in data

maledata <- na.omit(read.csv("data/MaleENA.csv"))
maledata$Individual <- as.character(maledata$Individual)
maledata$Clone <- as.character(maledata$Clone)

femaledata <- na.omit(read.csv("data/FemaleENA.csv"))
femaledata$Individual <- as.character(femaledata$Individual)
femaledata$Clone <- as.character(femaledata$Clone)
femaledata$Sex <- as.factor(femaledata$Sex)

##### transform data

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


### write new data frame for transformed data 
malePopulation <- as.character(maledata$Population)
maleIndividual <- as.character(maledata$Individual)
maleTransALL <- as.data.frame(cbind(maledatatransDtG,maledatatransDt10G,maledatatransPS,maledatatransLL,maledatatransDtP,maledatatransNP,malePopulation,maleIndividual))
maleTransALL$maledatatransDtG <- as.numeric(as.character(maleTransALL$maledatatransDtG))
maleTransALL$maledatatransDt10G <- as.numeric(as.character(maleTransALL$maledatatransDt10G))
maleTransALL$maledatatransPS <- as.numeric(as.character(maleTransALL$maledatatransPS))
maleTransALL$maledatatransLL <- as.numeric(as.character(maleTransALL$maledatatransLL))
maleTransALL$maledatatransNP <- as.numeric(as.character(maleTransALL$maledatatransNP))
maleTransALL$maledatatransDtP <- as.numeric(as.character(maleTransALL$maledatatransDtP))

femalePopulation <- as.character(femaledata$Population)
femaleIndividual <- as.character(femaledata$Individual)
femaleTransALL <- as.data.frame(cbind(femaledatatransDtG,femaledatatransDt10G,femaledatatransPS,femaledatatransLL,femalePopulation,femaleIndividual))
femaleTransALL$femaledatatransDtG <- as.numeric(as.character(femaleTransALL$femaledatatransDtG))
femaleTransALL$femaledatatransDt10G <- as.numeric(as.character(femaleTransALL$femaledatatransDt10G))
femaleTransALL$femaledatatransPS <- as.numeric(as.character(femaleTransALL$femaledatatransPS))
femaleTransALL$femaledatatransLL <- as.numeric(as.character(femaleTransALL$femaledatatransLL))


### averaging the traits for cca 
maleTransTraitsAvg <- aggregate(maleTransALL, by=list(maleTransALL$malePopulation, maleTransALL$maleIndividual), FUN=mean)
femaleTransTraitsAvg <- aggregate(femaleTransALL, by=list(femaleTransALL$femalePopulation, femaleTransALL$femaleIndividual), FUN=mean)

maleTransTraitsAvgSort <- maleTransTraitsAvg[order(maleTransTraitsAvg$Group.1),,]
femaleTransTraitsAvgSort <- femaleTransTraitsAvg[order(femaleTransTraitsAvg$Group.1),,] 


## final data for use in constrained ordination

ordinationMaleTransTraits <- cbind(maleTransTraitsAvgSort[3:8])
colnames(ordinationMaleTransTraits) <- c("DtG", "Dt10G", "PS", "LL","DtP","NP")
ordinationMaleTransTraits<- ordinationMaleTransTraits[ordinationMaleTransTraits$NP >= 0, ]

ordinationFemaleTransTraits <- cbind(femaleTransTraitsAvgSort[3:6])
colnames(ordinationFemaleTransTraits) <- c("DtG", "Dt10G", "PS", "LL")

#### environmental variables for use in cca
TransENVM <- read.csv("data/ENA_TransAvg_ENV_M.csv")
TransENVF <- read.csv("data/ENA_TransAvg_ENV_F.csv")

### environmental variables including population for use in graphing 
TransENVMwPop <- read.csv("data/ENA_TransAvgPop_ENV_M.csv") 
TransENVFwPop <- read.csv("data/ENA_TransAvgwPop_ENV_F.csv") 


##############################################################


## males 

##testing for any significance of env vars
ENAordMALL <- cca(formula = ordinationMaleTransTraits ~ ., data = TransENVM) 
anova(ENAordMALL)

### determining the best environmental variables to use in cca 
mod1 <- cca(ordinationMaleTransTraits ~ ., TransENVM)
mod0 <- cca(ordinationMaleTransTraits ~ 1, TransENVM)
mod <- step(mod0, scope=formula(mod1), test="perm")
mod
#cca(formula = ordinationMaleTransTraits ~ June + April, data = TransENVM)


ENAordM <- cca(formula = ordinationMaleTransTraits ~ June + May + April + Latitude, data = TransENVM)
ENAordM
summary(ENAordM)
anova(ENAordM)
anova(ENAordM, by = "term")
anova(ENAordM, by = "mar")
anova(ENAordM, by = "axis")


### male cca plot

colvec <- c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
            "deeppink")
popvec <- c("A","D","E","N","R","S","W")

png("figures/male_CCA_2019.png", width = 8, height = 6, units = 'in', res = 300)

par(mar=c(5.1,5.1,4.1,2.1))
ENAordMALL <- cca(formula = ordinationMaleTransTraits ~ ., data = TransENVM) 
plot(ENAordMALL,xlim=c(-10,10), ylim=c(-10,10),type="n", main="Males", 
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
with(TransENVMwPop, points(ENAordMALL, display="sites", pch=21, col="gray86", 
                           bg="gray86"))
with(TransENVMwPop, ordiellipse(ENAordMALL, groups=Population,draw="polygon",
                                col=colvec, kind="sd",lwd=0.1,alpha=50))
with(TransENVMwPop, ordiellipse(ENAordMALL, groups=Population,draw="line",
                                col=colvec,kind="sd",lwd=2,alpha=50))
ef <- envfit(ENAordMALL,TransENVM,perm=999)
plot(ef,col="black",lwd=10, cex=1.5)
with(ef, legend("topright", legend = popvec, bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

dev.off()



## females

##testing for any significance of env vars
ENAordFALL <- cca(formula = ordinationFemaleTransTraits ~ ., data = TransENVF) 
anova(ENAordFALL)

#### because the cca with all traits is not significant we should
## not proceed with model selection for females 


mod1 <- cca(ordinationFemaleTransTraits ~ ., TransENVF)
mod0 <- cca(ordinationFemaleTransTraits ~ 1, TransENVF)
mod <- step(mod0, scope=formula(mod1), test="perm")
mod
#Call: cca(formula = ordinationFemaleTransTraits ~ 1, data = TransENVF)


ENAordF <- cca(formula = ordinationFemaleTransTraits ~ 1, data = TransENVF) 
ENAordF
summary(ENAordF)

anova(ENAordF, by = "term")
anova(ENAordF, by = "mar")
anova(ENAordF, by = "axis")


### female cca plot

colvec <- c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
            "deeppink")

popvec <- c("A","D","E","N","R","S","W")


png("figures/female_CCA_2019.png", width = 8, height = 6, units = 'in', res = 300)

par(mar=c(5.1,5.1,4.1,2.1))
ENAordFALL <- cca(formula = ordinationFemaleTransTraits ~ ., data = TransENVF) 
plot(ENAordFALL,xlim=c(-10,10), ylim=c(-10,10),type="n", main="Females", 
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
with(TransENVFwPop, points(ENAordFALL, display="sites", pch=21, col="gray86", 
                           bg="gray86"))
with(TransENVFwPop, ordiellipse(ENAordFALL, groups=Population,draw="polygon",
                                col=colvec, kind="sd",lwd=0.1,alpha=50))
with(TransENVFwPop, ordiellipse(ENAordFALL, groups=Population,draw="line",
                                col=colvec,kind="sd",lwd=2,alpha=50))
ef <- envfit(ENAordFALL,TransENVF,perm=999)
plot(ef,col="black",lwd=10, cex=1.5)
with(ef, legend("topright", legend = popvec, bty = "n",
                col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

dev.off()








#### cca of males using only shared traits


ordinationMaleTransTraits <- cbind(maleTransTraitsAvgSort[3:6])
colnames(ordinationMaleTransTraits) <- c("DtG", "Dt10G", "PS", "LL")
TransENVM <- read.csv("data/ENA_TransAvg_ENV_M.csv")
TransENVMwPop <- read.csv("data/ENA_TransAvgPop_ENV_M.csv") 

ENAordMALL <- cca(formula = ordinationMaleTransTraits ~ ., data = TransENVM) 
anova(ENAordMALL)

anova(ENAordMALL, by = "term")
anova(ENAordMALL, by = "mar")
anova(ENAordMALL, by = "axis")


### determining the best environmental variables to use in cca 
mod1 <- cca(ordinationMaleTransTraits ~ ., TransENVM)
mod0 <- cca(ordinationMaleTransTraits ~ 1, TransENVM)
mod <- step(mod0, scope=formula(mod1), test="perm")
mod
#cca(formula = ordinationMaleTransTraits ~ June + April, data = TransENVM)


ENAordM <- cca(formula = ordinationMaleTransTraits ~ June + April, data = TransENVM)
ENAordM
summary(ENAordM)
anova(ENAordM)
anova(ENAordM, by = "term")
anova(ENAordM, by = "mar")
anova(ENAordM, by = "axis")

colvec <- c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
            "deeppink")

popvec <- c("A","D","E","N","R","S","W")

png("figures/malewodevo_CCA_2019.png", width = 8, height = 6, units = 'in', res = 300)

par(mar=c(5.1,5.1,4.1,2.1))
ENAordMALL <- cca(formula = ordinationMaleTransTraits ~ ., data = TransENVM) 
plot(ENAordMALL,xlim=c(-10,10), ylim=c(-10,10),type="n", main="Males", 
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
with(TransENVMwPop, points(ENAordMALL, display="sites", pch=21, col="gray86", 
                           bg="gray86"))
with(TransENVMwPop, ordiellipse(ENAordMALL, groups=Population,draw="polygon",
                                col=colvec, kind="sd",lwd=0.1,alpha=50))
with(TransENVMwPop, ordiellipse(ENAordMALL, groups=Population,draw="line",
                                col=colvec,kind="sd",lwd=2,alpha=50))
ef <- envfit(ENAordMALL,TransENVM,perm=999)
plot(ef,col="black",lwd=10, cex=1.5)
with(ef, legend("topright", legend = popvec, bty = "n",
                col = colvec, pch = 21, pt.bg = colvec, cex = 1.5))

dev.off()


