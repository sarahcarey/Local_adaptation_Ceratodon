
## Within population trait correlations

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


#####################
## running correlations

# males

## global
maletransCor <- cor(maleTransALL[1:6], use="everything")
colnames(maletransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maletransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")

## by population
maleAtrans <- subset(maleTransALL, malePopulation %in% "A")
maleAtransCor <- cor(maleAtrans[1:6], use="everything")
colnames(maleAtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maleAtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleDtrans <- subset(maleTransALL, malePopulation %in% "D")
maleDtransCor <- cor(maleDtrans[1:6], use="everything")
colnames(maleDtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maleDtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleEtrans <- subset(maleTransALL, malePopulation %in% "E")
maleEtransCor <- cor(maleEtrans[1:6], use="everything")
colnames(maleEtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maleEtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleNtrans <- subset(maleTransALL, malePopulation %in% "N")
maleNtransCor <- cor(maleNtrans[1:6], use="everything")
colnames(maleNtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maleNtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleRtrans <- subset(maleTransALL, malePopulation %in% "R")
maleRtransCor <- cor(maleRtrans[1:6], use="everything")
colnames(maleRtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maleRtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleStrans <- subset(maleTransALL, malePopulation %in% "S")
maleStransCor <- cor(maleStrans[1:6], use="everything")
colnames(maleStransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maleStransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")

maleWtrans <- subset(maleTransALL, malePopulation %in% "W")
maleWtransCor <- cor(maleAtrans[1:6], use="everything")
colnames(maleWtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")
rownames(maleWtransCor) <- c("DtG","Dt10G","PS","LL","DtP","NP")


######################################

# females

## global

femaletransCor <- cor(femaleTransALL[1:4], use="everything")
colnames(femaletransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaletransCor) <- c("DtG","Dt10G","PS","LL")

## by population
femaleAtrans <- subset(femaleTransALL, femalePopulation %in% "A")
femaleAtransCor <- cor(femaleAtrans[1:4], use="everything")
colnames(femaleAtransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaleAtransCor) <- c("DtG","Dt10G","PS","LL")

femaleDtrans <- subset(femaleTransALL, femalePopulation %in% "D")
femaleDtransCor <- cor(femaleDtrans[1:4], use="everything")
colnames(femaleDtransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaleDtransCor) <- c("DtG","Dt10G","PS","LL")

femaleEtrans <- subset(femaleTransALL, femalePopulation %in% "E")
femaleEtransCor <- cor(femaleEtrans[1:4], use="everything")
colnames(femaleEtransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaleEtransCor) <- c("DtG","Dt10G","PS","LL")

femaleNtrans <- subset(femaleTransALL, femalePopulation %in% "N")
femaleNtransCor <- cor(femaleNtrans[1:4], use="everything")
colnames(femaleNtransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaleNtransCor) <- c("DtG","Dt10G","PS","LL")

femaleRtrans <- subset(femaleTransALL, femalePopulation %in% "R")
femaleRtransCor <- cor(femaleRtrans[1:4], use="everything")
colnames(femaleRtransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaleRtransCor) <- c("DtG","Dt10G","PS","LL")

femaleStrans <- subset(femaleTransALL, femalePopulation %in% "S")
femaleStransCor <- cor(femaleStrans[1:4], use="everything")
colnames(femaleStransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaleStransCor) <- c("DtG","Dt10G","PS","LL")

femaleWtrans <- subset(femaleTransALL, femalePopulation %in% "W")
femaleWtransCor <- cor(femaleWtrans[1:4], use="everything")
colnames(femaleWtransCor) <- c("DtG","Dt10G","PS","LL")
rownames(femaleWtransCor) <- c("DtG","Dt10G","PS","LL")


########################

## writing outputs

## global
write.csv(maletransCor,"tables/maleTraitCorrGlobal.csv")
write.csv(femaletransCor,"tables/femaleTraitCorrGlobal.csv")


## by population
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

femaleTraitCorrs <- as.data.frame(rbind(femaleAtransCor,femaleDtransCor,femaleEtransCor,femaleNtransCor,
                          femaleRtransCor,femaleStransCor,femaleWtransCor))

femaleTraitCorrswPop <- as.data.frame(cbind(popnamesF,corrRowNamesF,femaleTraitCorrs))

write.csv(femaleTraitCorrswPop,"tables/femaleTraitCorrswPop.csv")

maleTraitCorrs <- as.data.frame(rbind(maleAtransCor,maleDtransCor,maleEtransCor,maleNtransCor,
                                        maleRtransCor,maleStransCor,maleWtransCor))

maleTraitCorrswPop <- as.data.frame(cbind(popnamesM,corrRowNamesM,maleTraitCorrs))

write.csv(maleTraitCorrswPop,"tables/maleTraitCorrswPop.csv")

