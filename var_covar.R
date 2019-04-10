
## Variance and covariance matrices

## analyses for McDaniel et al done with R version 3.5.3

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


## data needs to be scaled
femaleTransALLscaled <- as.data.frame(scale(femaleTransALL[1:4]))
maleTransALLscaled <- as.data.frame(scale(maleTransALL[1:6]))


################################################
#### variance covariance matrices

# global male
CoVarMaleGlobal <- cov(maleTransALLscaled[1:6], use = "complete.obs")
colnames(CoVarMaleGlobal) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleGlobal) <- c("DtG","Dt10G","PS","LL","DtP","NP")

# global female
CoVarFemaleGlobal <- cov(femaleTransALLscaled[1:4], use = "complete.obs")
colnames(CoVarFemaleGlobal) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleGlobal) <- c("DtG","Dt10G","PS","LL")


# male within pops
maleAtrans <- subset(maleTransALLscaled , malePopulation %in% "A")
CoVarMaleA <- cov(maleAtrans[1:6], use = "complete.obs")
colnames(CoVarMaleA) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleA) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleDtrans <- subset(maleTransALLscaled , malePopulation %in% "D")
CoVarMaleD <- cov(maleDtrans[1:6], use = "complete.obs")
colnames(CoVarMaleD) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleD) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleEtrans <- subset(maleTransALLscaled , malePopulation %in% "E")
CoVarMaleE <- cov(maleEtrans[1:6], use = "complete.obs")
colnames(CoVarMaleE) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleE) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleNtrans <- subset(maleTransALLscaled , malePopulation %in% "N")
CoVarMaleN <- cov(maleNtrans[1:6], use = "complete.obs")
colnames(CoVarMaleN) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleN) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleRtrans <- subset(maleTransALLscaled , malePopulation %in% "R")
CoVarMaleR <- cov(maleRtrans[1:6], use = "complete.obs")
colnames(CoVarMaleR) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleR) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleStrans <- subset(maleTransALLscaled , malePopulation %in% "S")
CoVarMaleS <- cov(maleStrans[1:6], use = "complete.obs")
colnames(CoVarMaleS) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleS) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleWtrans <- subset(maleTransALLscaled , malePopulation %in% "W")
CoVarMaleW <- cov(maleWtrans[1:6], use = "complete.obs")
colnames(CoVarMaleW) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(CoVarMaleW) <- c("DtG","Dt10G","PS","LL","DtP","NP")

# female within pops
femaleAtrans <- subset(femaleTransALLscaled , femalePopulation %in% "A")
CoVarFemaleA <- cov(femaleAtrans[1:4], use = "complete.obs")
colnames(CoVarFemaleA) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleA) <- c("DtG","Dt10G","PS","LL")

femaleDtrans <- subset(femaleTransALLscaled , femalePopulation %in% "D")
CoVarFemaleD <- cov(femaleDtrans[1:4], use = "complete.obs")
colnames(CoVarFemaleD) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleD) <- c("DtG","Dt10G","PS","LL")

femaleEtrans <- subset(femaleTransALLscaled , femalePopulation %in% "E")
CoVarFemaleE <- cov(femaleEtrans[1:4], use = "complete.obs")
colnames(CoVarFemaleE) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleE) <- c("DtG","Dt10G","PS","LL")

femaleNtrans <- subset(femaleTransALLscaled , femalePopulation %in% "N")
CoVarFemaleN <- cov(femaleNtrans[1:4], use = "complete.obs")
colnames(CoVarFemaleN) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleN) <- c("DtG","Dt10G","PS","LL")

femaleRtrans <- subset(femaleTransALLscaled , femalePopulation %in% "R")
CoVarFemaleR <- cov(femaleRtrans[1:4], use = "complete.obs")
colnames(CoVarFemaleR) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleR) <- c("DtG","Dt10G","PS","LL")

femaleStrans <- subset(femaleTransALLscaled , femalePopulation %in% "S")
CoVarFemaleS <- cov(femaleStrans[1:4], use = "complete.obs")
colnames(CoVarFemaleS) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleS) <- c("DtG","Dt10G","PS","LL")

femaleWtrans <- subset(femaleTransALLscaled , femalePopulation %in% "W")
CoVarFemaleW <- cov(femaleWtrans[1:4], use = "complete.obs")
colnames(CoVarFemaleW) <- c("DtG","Dt10G","PS","LL")
rownames(CoVarFemaleW) <- c("DtG","Dt10G","PS","LL")


########################
## writing outputs

popnamesF <- c("A","","","",
              "D","","","",
              "E","","","",
              "N","","","",
              "R","","","",
              "S","","","",
              "W","","","")

popnamesM <- c("A","","","","","",
               "D","","","","","",
               "E","","","","","",
               "N","","","","","",
               "R","","","","","",
               "S","","","","","",
               "W","","","","","")

corrRowNamesF <- c("DtG","Dt10G","PS","LL",
                  "DtG","Dt10G","PS","LL",
                  "DtG","Dt10G","PS","LL",
                  "DtG","Dt10G","PS","LL",
                  "DtG","Dt10G","PS","LL",
                  "DtG","Dt10G","PS","LL",
                  "DtG","Dt10G","PS","LL")

corrRowNamesM <- c("DtG","Dt10G","PS","LL","DtP","NP",
                   "DtG","Dt10G","PS","LL","DtP","NP",
                   "DtG","Dt10G","PS","LL","DtP","NP",
                   "DtG","Dt10G","PS","LL","DtP","NP",
                   "DtG","Dt10G","PS","LL","DtP","NP",
                   "DtG","Dt10G","PS","LL","DtP","NP",
                   "DtG","Dt10G","PS","LL","DtP","NP")

## female population variance/covariances

femaleTraitVarCovar <- as.data.frame(rbind(CoVarFemaleA,CoVarFemaleD,CoVarFemaleE,CoVarFemaleN,
                                          CoVarFemaleR,CoVarFemaleS,CoVarFemaleW))

femaleTraitVarCovarwPop <- as.data.frame(cbind(popnamesF,corrRowNamesF,femaleTraitVarCovar))

write.csv(femaleTraitVarCovarwPop,"tables/femaleTraitVarCovarwPop.csv")


## male population variance/covariances
maleTraitVarCovar  <- as.data.frame(rbind(CoVarMaleA,CoVarMaleD,CoVarMaleE,CoVarMaleN,
                                          CoVarMaleR,CoVarMaleS,CoVarMaleW))

maleTraitVarCovarwPop <- as.data.frame(cbind(popnamesM,corrRowNamesM,maleTraitVarCovar))

write.csv(maleTraitVarCovarwPop,"tables/maleTraitVarCovarwPop.csv")


## global male
write.csv(CoVarMaleGlobal,"tables/CoVarMaleGlobal.csv")

#global female
write.csv(CoVarFemaleGlobal,"tables/CoVarFemaleGlobal.csv")



