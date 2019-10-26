
## Least square means

## analyses for McDaniel et al done with R version 3.5.3
## lsmeans version 2.30

install.packages("lsmeans")
library("lsmeans")

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
maleTransALL$maledatatransDtP <- as.numeric(as.character(maleTransALL$maledatatransDtP))
maleTransALL$maledatatransNP <- as.numeric(as.character(maleTransALL$maledatatransNP))

femalePopulation <- as.character(femaledata$Population)
femaleIndividual <- as.character(femaledata$Individual)
femaleTransALL <- as.data.frame(cbind(femaledatatransDtG,femaledatatransDt10G,femaledatatransPS,femaledatatransLL,femalePopulation,femaleIndividual))
femaleTransALL$femaledatatransDtG <- as.numeric(as.character(femaleTransALL$femaledatatransDtG))
femaleTransALL$femaledatatransDt10G <- as.numeric(as.character(femaleTransALL$femaledatatransDt10G))
femaleTransALL$femaledatatransPS <- as.numeric(as.character(femaleTransALL$femaledatatransPS))
femaleTransALL$femaledatatransLL <- as.numeric(as.character(femaleTransALL$femaledatatransLL))


##########################################
## running lsm code and then backtransforming the output for easier interpretation
## cannot use SE because of transforming so we will use the 95% CL 

##males 
aovMalemaledatatransDtG <- aov(maledatatransDtG~malePopulation, data=maleTransALL)
lsmMalemaledatatransDtG <- lsmeans(aovMalemaledatatransDtG, "malePopulation")
lsmMalemaledatatransDtGtable <- summary(lsmMalemaledatatransDtG)
backtransMaleDtGlsm <- (lsmMalemaledatatransDtGtable$lsmean)^(-1/2)
backtransMaleDtGlower <- (lsmMalemaledatatransDtGtable$lower.CL)^(-1/2)
backtransMaleDtGupper <- (lsmMalemaledatatransDtGtable$upper.CL)^(-1/2)

aovMalemaledatatransDt10G <- aov(maledatatransDt10G~malePopulation, data=maleTransALL)
lsmMalemaledatatransDt10G <- lsmeans(aovMalemaledatatransDt10G, "malePopulation")
lsmMalemaledatatransDt10Gtable<- summary(lsmMalemaledatatransDt10G)
backtransMaleDt10Glsm <- (lsmMalemaledatatransDt10Gtable$lsmean)^(-2)
backtransMaleDt10Glower <- (lsmMalemaledatatransDt10Gtable$lower.CL)^(-2)
backtransMaleDt10Gupper <- (lsmMalemaledatatransDt10Gtable$upper.CL)^(-2)

aovMalemaledatatransPS <- aov(maledatatransPS~malePopulation, data=maleTransALL)
lsmMalemaledatatransPS <- lsmeans(aovMalemaledatatransPS, "malePopulation")
lsmMalemaledatatransPStable<- summary(lsmMalemaledatatransPS)
backtransMalePSlsm <- (lsmMalemaledatatransPStable$lsmean)
backtransMalePSlower <- (lsmMalemaledatatransPStable$lower.CL)
backtransMalePSupper <- (lsmMalemaledatatransPStable$upper.CL)

aovMalemaledatatransLL <- aov(maledatatransLL~malePopulation, data=maleTransALL)
lsmMalemaledatatransLL <- lsmeans(aovMalemaledatatransLL, "malePopulation")  
lsmMalemaledatatransLLtable<- summary(lsmMalemaledatatransLL)
backtransMaleLLlsm <- exp(lsmMalemaledatatransLLtable$lsmean)
backtransMaleLLlower <- exp(lsmMalemaledatatransLLtable$lower.CL)
backtransMaleLLupper <- exp(lsmMalemaledatatransLLtable$upper.CL)

aovMalemaledatatransDtP <- aov(maledatatransDtP~malePopulation, data=maleTransALL)
lsmMalemaledatatransDtP <- lsmeans(aovMalemaledatatransDtP, "malePopulation")  
lsmMalemaledatatransDtPtable<- summary(lsmMalemaledatatransDtP)
backtransMaleDtPlsm <- (lsmMalemaledatatransDtPtable$lsmean)^(1/2)
backtransMaleDtPlower <- (lsmMalemaledatatransDtPtable$lower.CL)^(1/2)
backtransMaleDtPupper <- (lsmMalemaledatatransDtPtable$upper.CL)^(1/2)

aovMalemaledatatransNP <- aov(maledatatransNP~malePopulation, data=maleTransALL)
lsmMalemaledatatransNP <- lsmeans(aovMalemaledatatransNP, "malePopulation")  
lsmMalemaledatatransNPtable<- summary(lsmMalemaledatatransNP)
backtransMaleNPlsm <- (lsmMalemaledatatransNPtable$lsmean)^(2)
backtransMaleNPlower <- (lsmMalemaledatatransNPtable$lower.CL)^(2)
backtransMaleNPupper <- (lsmMalemaledatatransNPtable$upper.CL)^(2)

##########################################

####females

aovFemalefemaledatatransDtG <- aov(femaledatatransDtG~femalePopulation, data=femaleTransALL)
lsmFemalefemaledatatransDtG <- lsmeans(aovFemalefemaledatatransDtG, "femalePopulation")
lsmFemalefemaledatatransDtGtable <- summary(lsmFemalefemaledatatransDtG)
backtransFemaleDtGlsm <- (lsmFemalefemaledatatransDtGtable$lsmean)^(-1/1.25)
backtransFemaleDtGlower <- (lsmFemalefemaledatatransDtGtable$lower.CL)^(-1/1.25)
backtransFemaleDtGupper <- (lsmFemalefemaledatatransDtGtable$upper.CL)^(-1/1.25)

aovFemalefemaledatatransDt10G <- aov(femaledatatransDt10G~femalePopulation, data=femaleTransALL)
lsmFemalefemaledatatransDt10G <- lsmeans(aovFemalefemaledatatransDt10G, "femalePopulation") 
lsmFemalefemaledatatransDt10Gtable<- summary(lsmFemalefemaledatatransDt10G)
backtransFemaleDt10Glsm <- (lsmFemalefemaledatatransDt10Gtable$lsmean)^(-4)
backtransFemaleDt10Glower <- (lsmFemalefemaledatatransDt10Gtable$lower.CL)^(-4)
backtransFemaleDt10Gupper <- (lsmFemalefemaledatatransDt10Gtable$upper.CL)^(-4)

aovFemalefemaledatatransPS <- aov(femaledatatransPS~femalePopulation, data=femaleTransALL)
lsmFemalefemaledatatransPS <- lsmeans(aovFemalefemaledatatransPS, "femalePopulation") 
lsmFemalefemaledatatransPStable<- summary(lsmFemalefemaledatatransPS)
backtransFemalePSlsm <- (lsmFemalefemaledatatransPStable$lsmean)^(1/0.75)
backtransFemalePSlower <- (lsmFemalefemaledatatransPStable$lower.CL)^(1/0.75)
backtransFemalePSupper <- (lsmFemalefemaledatatransPStable$upper.CL)^(1/0.75)

aovFemalefemaledatatransLL <- aov(femaledatatransLL~femalePopulation, data=femaleTransALL)
lsmFemalefemaledatatransLL <- lsmeans(aovFemalefemaledatatransLL, "femalePopulation")  
lsmFemalefemaledatatransLLtable<- summary(lsmFemalefemaledatatransLL)
backtransFemaleLLlsm <- (lsmFemalefemaledatatransLLtable$lsmean)^(4)
backtransFemaleLLlower <- (lsmFemalefemaledatatransLLtable$lower.CL)^(4)
backtransFemaleLLupper <- (lsmFemalefemaledatatransLLtable$upper.CL)^(4)


##########################################
##write tables

DtGlsm <- as.data.frame(cbind(backtransMaleDtGlsm,backtransMaleDtGupper,backtransMaleDtGlower,backtransFemaleDtGlsm,backtransFemaleDtGupper,backtransFemaleDtGlower))
colnames(DtGlsm) <- c("Males", "", "", "Females", "", "")
rownames(DtGlsm) <- c("A","D","E","N","R","S","W")

Dt10Glsm <- as.data.frame(cbind(backtransMaleDt10Glsm,backtransMaleDt10Gupper,backtransMaleDt10Glower,backtransFemaleDt10Glsm,backtransFemaleDt10Gupper,backtransFemaleDt10Glower))
colnames(Dt10Glsm) <- c("Males", "", "", "Females", "", "")
rownames(Dt10Glsm) <- c("A","D","E","N","R","S","W")

PSlsm <- as.data.frame(cbind(backtransMalePSlsm,backtransMalePSupper,backtransMalePSlower,backtransFemalePSlsm,backtransFemalePSupper,backtransFemalePSlower))
colnames(PSlsm) <- c("Males", "", "", "Females", "", "")
rownames(PSlsm) <- c("A","D","E","N","R","S","W")

LLlsm <- as.data.frame(cbind(backtransMaleLLlsm,backtransMaleLLupper,backtransMaleLLlower,backtransFemaleLLlsm,backtransFemaleLLupper,backtransFemaleLLlower))
colnames(LLlsm) <- c("Males", "", "", "Females", "", "")
rownames(LLlsm) <- c("A","D","E","N","R","S","W")

Dtplsm <- as.data.frame(cbind(backtransMaleDtPlsm,backtransMaleDtPupper,backtransMaleDtPlower,"","",""))
colnames(Dtplsm) <- c("Males", "", "", "Females", "", "")
rownames(Dtplsm) <- c("A","D","E","N","R","S","W")

NPlsm <- as.data.frame(cbind(backtransMaleNPlsm,backtransMaleNPupper,backtransMaleNPlower,"","",""))
colnames(NPlsm) <- c("Males", "", "", "Females", "", "")
rownames(NPlsm) <- c("A","D","E","N","R","S","W")

traitCol <- c("DtG","","","","","","",
                 "Dt10G","","","","","","",
                 "PS","","","","","","",
                 "LL","","","","","","",
              "DtP","","","","","","",
              "NP","","","","","","")

lsmALL <- rbind(DtGlsm,Dt10Glsm,PSlsm,LLlsm,Dtplsm,NPlsm)
lsmALLwTrait <- cbind(traitCol,lsmALL)

write.csv(lsmALLwTrait,"tables/lsmALLwTrait.csv")

