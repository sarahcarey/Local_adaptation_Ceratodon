
## Coefficients of variation

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


#maleA
maleAtrans <- subset(maleTransALL, malePopulation %in% "A")

CV <- function(m, std){(std/m)}
m <- mean(maleAtrans$maledatatransDtG, na.rm = TRUE)
std <- sd(maleAtrans$maledatatransDtG, na.rm = TRUE)
CVmaleAtransDtG <- CV(m, std)

m <- mean(maleAtrans$maledatatransDt10G, na.rm = TRUE)
std <- sd(maleAtrans$maledatatransDt10G, na.rm = TRUE)
CVmaleAtransDt10G <- CV(m, std)

m <- mean(maleAtrans$maledatatransPS, na.rm = TRUE)
std <- sd(maleAtrans$maledatatransPS, na.rm = TRUE)
CVmaleAtransPS <- CV(m, std)

m <- mean(maleAtrans$maledatatransLL, na.rm = TRUE)
std <- sd(maleAtrans$maledatatransLL, na.rm = TRUE)
CVmaleAtransLL <- CV(m, std)

m <- mean(maleAtrans$maledatatransDtP, na.rm = TRUE)
std <- sd(maleAtrans$maledatatransDtP, na.rm = TRUE)
CVmaleAtransDtP <- CV(m, std)

m <- mean(maleAtrans$maledatatransNP, na.rm = TRUE)
std <- sd(maleAtrans$maledatatransNP, na.rm = TRUE)
CVmaleAtransNP <- CV(m, std)



#maleDtrans
maleDtrans <- subset(maleTransALL, malePopulation %in% "D")

CV <- function(m, std){(std/m)}
m <- mean(maleDtrans$maledatatransDtG, na.rm = TRUE)
std <- sd(maleDtrans$maledatatransDtG, na.rm = TRUE)
CVmaleDtransDtG <- CV(m, std)

m <- mean(maleDtrans$maledatatransDt10G, na.rm = TRUE)
std <- sd(maleDtrans$maledatatransDt10G, na.rm = TRUE)
CVmaleDtransDt10G <- CV(m, std)

m <- mean(maleDtrans$maledatatransPS, na.rm = TRUE)
std <- sd(maleDtrans$maledatatransPS, na.rm = TRUE)
CVmaleDtransPS <- CV(m, std)

m <- mean(maleDtrans$maledatatransLL, na.rm = TRUE)
std <- sd(maleDtrans$maledatatransLL, na.rm = TRUE)
CVmaleDtransLL <- CV(m, std)

m <- mean(maleDtrans$maledatatransDtP, na.rm = TRUE)
std <- sd(maleDtrans$maledatatransDtP, na.rm = TRUE)
CVmaleDtransDtP <- CV(m, std)

m <- mean(maleDtrans$maledatatransNP, na.rm = TRUE)
std <- sd(maleDtrans$maledatatransNP, na.rm = TRUE)
CVmaleDtransNP <- CV(m, std)


#maleEtrans
maleEtrans <- subset(maleTransALL, malePopulation %in% "E")

CV <- function(m, std){(std/m)}
m <- mean(maleEtrans$maledatatransDtG, na.rm = TRUE)
std <- sd(maleEtrans$maledatatransDtG, na.rm = TRUE)
CVmaleEtransDtG <- CV(m, std)

m <- mean(maleEtrans$maledatatransDt10G, na.rm = TRUE)
std <- sd(maleEtrans$maledatatransDt10G, na.rm = TRUE)
CVmaleEtransDt10G <- CV(m, std)

m <- mean(maleEtrans$maledatatransPS, na.rm = TRUE)
std <- sd(maleEtrans$maledatatransPS, na.rm = TRUE)
CVmaleEtransPS <- CV(m, std)

m <- mean(maleEtrans$maledatatransLL, na.rm = TRUE)
std <- sd(maleEtrans$maledatatransLL, na.rm = TRUE)
CVmaleEtransLL <- CV(m, std)

m <- mean(maleEtrans$maledatatransDtP, na.rm = TRUE)
std <- sd(maleEtrans$maledatatransDtP, na.rm = TRUE)
CVmaleEtransDtP <- CV(m, std)

m <- mean(maleEtrans$maledatatransNP, na.rm = TRUE)
std <- sd(maleEtrans$maledatatransNP, na.rm = TRUE)
CVmaleEtransNP <- CV(m, std)


#maleNtrans
maleNtrans <- subset(maleTransALL, malePopulation %in% "N")

CV <- function(m, std){(std/m)}
m <- mean(maleNtrans$maledatatransDtG, na.rm = TRUE)
std <- sd(maleNtrans$maledatatransDtG, na.rm = TRUE)
CVmaleNtransDtG <- CV(m, std)

m <- mean(maleNtrans$maledatatransDt10G, na.rm = TRUE)
std <- sd(maleNtrans$maledatatransDt10G, na.rm = TRUE)
CVmaleNtransDt10G <- CV(m, std)

m <- mean(maleNtrans$maledatatransPS, na.rm = TRUE)
std <- sd(maleNtrans$maledatatransPS, na.rm = TRUE)
CVmaleNtransPS <- CV(m, std)

m <- mean(maleNtrans$maledatatransLL, na.rm = TRUE)
std <- sd(maleNtrans$maledatatransLL, na.rm = TRUE)
CVmaleNtransLL <- CV(m, std)

m <- mean(maleNtrans$maledatatransDtP, na.rm = TRUE)
std <- sd(maleNtrans$maledatatransDtP, na.rm = TRUE)
CVmaleNtransDtP <- CV(m, std)

m <- mean(maleNtrans$maledatatransNP, na.rm = TRUE)
std <- sd(maleNtrans$maledatatransNP, na.rm = TRUE)
CVmaleNtransNP <- CV(m, std)


#maleRtrans
maleRtrans <- subset(maleTransALL, malePopulation %in% "R")

CV <- function(m, std){(std/m)}
m <- mean(maleRtrans$maledatatransDtG, na.rm = TRUE)
std <- sd(maleRtrans$maledatatransDtG, na.rm = TRUE)
CVmaleRtransDtG <- CV(m, std)

m <- mean(maleRtrans$maledatatransDt10G, na.rm = TRUE)
std <- sd(maleRtrans$maledatatransDt10G, na.rm = TRUE)
CVmaleRtransDt10G <- CV(m, std)

m <- mean(maleRtrans$maledatatransPS, na.rm = TRUE)
std <- sd(maleRtrans$maledatatransPS, na.rm = TRUE)
CVmaleRtransPS <- CV(m, std)

m <- mean(maleRtrans$maledatatransLL, na.rm = TRUE)
std <- sd(maleRtrans$maledatatransLL, na.rm = TRUE)
CVmaleRtransLL <- CV(m, std)

m <- mean(maleRtrans$maledatatransDtP, na.rm = TRUE)
std <- sd(maleRtrans$maledatatransDtP, na.rm = TRUE)
CVmaleRtransDtP <- CV(m, std)

m <- mean(maleRtrans$maledatatransNP, na.rm = TRUE)
std <- sd(maleRtrans$maledatatransNP, na.rm = TRUE)
CVmaleRtransNP <- CV(m, std)


#maleStrans
maleStrans <- subset(maleTransALL, malePopulation %in% "S")

CV <- function(m, std){(std/m)}
m <- mean(maleStrans$maledatatransDtG, na.rm = TRUE)
std <- sd(maleStrans$maledatatransDtG, na.rm = TRUE)
CVmaleStransDtG <- CV(m, std)

m <- mean(maleStrans$maledatatransDt10G, na.rm = TRUE)
std <- sd(maleStrans$maledatatransDt10G, na.rm = TRUE)
CVmaleStransDt10G <- CV(m, std)

m <- mean(maleStrans$maledatatransPS, na.rm = TRUE)
std <- sd(maleStrans$maledatatransPS, na.rm = TRUE)
CVmaleStransPS <- CV(m, std)

m <- mean(maleStrans$maledatatransLL, na.rm = TRUE)
std <- sd(maleStrans$maledatatransLL, na.rm = TRUE)
CVmaleStransLL <- CV(m, std)

m <- mean(maleStrans$maledatatransDtP, na.rm = TRUE)
std <- sd(maleStrans$maledatatransDtP, na.rm = TRUE)
CVmaleStransDtP <- CV(m, std)

m <- mean(maleStrans$maledatatransNP, na.rm = TRUE)
std <- sd(maleStrans$maledatatransNP, na.rm = TRUE)
CVmaleStransNP <- CV(m, std)


#maleWtrans
maleWtrans <- subset(maleTransALL, malePopulation %in% "W")

CV <- function(m, std){(std/m)}
m <- mean(maleWtrans$maledatatransDtG, na.rm = TRUE)
std <- sd(maleWtrans$maledatatransDtG, na.rm = TRUE)
CVmaleWtransDtG <- CV(m, std)

m <- mean(maleWtrans$maledatatransDt10G, na.rm = TRUE)
std <- sd(maleWtrans$maledatatransDt10G, na.rm = TRUE)
CVmaleWtransDt10G <- CV(m, std)

m <- mean(maleWtrans$maledatatransPS, na.rm = TRUE)
std <- sd(maleWtrans$maledatatransPS, na.rm = TRUE)
CVmaleWtransPS <- CV(m, std)

m <- mean(maleWtrans$maledatatransLL, na.rm = TRUE)
std <- sd(maleWtrans$maledatatransLL, na.rm = TRUE)
CVmaleWtransLL <- CV(m, std)

m <- mean(maleWtrans$maledatatransDtP, na.rm = TRUE)
std <- sd(maleWtrans$maledatatransDtP, na.rm = TRUE)
CVmaleWtransDtP <- CV(m, std)

m <- mean(maleWtrans$maledatatransNP, na.rm = TRUE)
std <- sd(maleWtrans$maledatatransNP, na.rm = TRUE)
CVmaleWtransNP <- CV(m, std)


#femaleAtrans
femaleAtrans <- subset(femaleTransALL, femalePopulation %in% "A")

CV <- function(m, std){(std/m)}
m <- mean(femaleAtrans$femaledatatransDtG, na.rm = TRUE)
std <- sd(femaleAtrans$femaledatatransDtG, na.rm = TRUE)
CVfemaleAtransDtG <- CV(m, std)

m <- mean(femaleAtrans$femaledatatransDt10G, na.rm = TRUE)
std <- sd(femaleAtrans$femaledatatransDt10G, na.rm = TRUE)
CVfemaleAtransDt10G <- CV(m, std)

m <- mean(femaleAtrans$femaledatatransPS, na.rm = TRUE)
std <- sd(femaleAtrans$femaledatatransPS, na.rm = TRUE)
CVfemaleAtransPS <- CV(m, std)

m <- mean(femaleAtrans$femaledatatransLL, na.rm = TRUE)
std <- sd(femaleAtrans$femaledatatransLL, na.rm = TRUE)
CVfemaleAtransLL <- CV(m, std)


#femaleDtrans
femaleDtrans <- subset(femaleTransALL, femalePopulation %in% "D")

CV <- function(m, std){(std/m)}
m <- mean(femaleDtrans$femaledatatransDtG, na.rm = TRUE)
std <- sd(femaleDtrans$femaledatatransDtG, na.rm = TRUE)
CVfemaleDtransDtG <- CV(m, std)

m <- mean(femaleDtrans$femaledatatransDt10G, na.rm = TRUE)
std <- sd(femaleDtrans$femaledatatransDt10G, na.rm = TRUE)
CVfemaleDtransDt10G <- CV(m, std)

m <- mean(femaleDtrans$femaledatatransPS, na.rm = TRUE)
std <- sd(femaleDtrans$femaledatatransPS, na.rm = TRUE)
CVfemaleDtransPS <- CV(m, std)

m <- mean(femaleDtrans$femaledatatransLL, na.rm = TRUE)
std <- sd(femaleDtrans$femaledatatransLL, na.rm = TRUE)
CVfemaleDtransLL <- CV(m, std)


#femaleEtrans
femaleEtrans <- subset(femaleTransALL, femalePopulation %in% "E")

CV <- function(m, std){(std/m)}
m <- mean(femaleEtrans$femaledatatransDtG, na.rm = TRUE)
std <- sd(femaleEtrans$femaledatatransDtG, na.rm = TRUE)
CVfemaleEtransDtG <- CV(m, std)

m <- mean(femaleEtrans$femaledatatransDt10G, na.rm = TRUE)
std <- sd(femaleEtrans$femaledatatransDt10G, na.rm = TRUE)
CVfemaleEtransDt10G <- CV(m, std)

m <- mean(femaleEtrans$femaledatatransPS, na.rm = TRUE)
std <- sd(femaleEtrans$femaledatatransPS, na.rm = TRUE)
CVfemaleEtransPS <- CV(m, std)

m <- mean(femaleEtrans$femaledatatransLL, na.rm = TRUE)
std <- sd(femaleEtrans$femaledatatransLL, na.rm = TRUE)
CVfemaleEtransLL <- CV(m, std)


#femaleNtrans
femaleNtrans <- subset(femaleTransALL, femalePopulation %in% "N")

CV <- function(m, std){(std/m)}
m <- mean(femaleNtrans$femaledatatransDtG, na.rm = TRUE)
std <- sd(femaleNtrans$femaledatatransDtG, na.rm = TRUE)
CVfemaleNtransDtG <- CV(m, std)

m <- mean(femaleNtrans$femaledatatransDt10G, na.rm = TRUE)
std <- sd(femaleNtrans$femaledatatransDt10G, na.rm = TRUE)
CVfemaleNtransDt10G <- CV(m, std)

m <- mean(femaleNtrans$femaledatatransPS, na.rm = TRUE)
std <- sd(femaleNtrans$femaledatatransPS, na.rm = TRUE)
CVfemaleNtransPS <- CV(m, std)

m <- mean(femaleNtrans$femaledatatransLL, na.rm = TRUE)
std <- sd(femaleNtrans$femaledatatransLL, na.rm = TRUE)
CVfemaleNtransLL <- CV(m, std)


#femaleRtrans
femaleRtrans <- subset(femaleTransALL, femalePopulation %in% "R")

CV <- function(m, std){(std/m)}
m <- mean(femaleRtrans$femaledatatransDtG, na.rm = TRUE)
std <- sd(femaleRtrans$femaledatatransDtG, na.rm = TRUE)
CVfemaleRtransDtG <- CV(m, std)

m <- mean(femaleRtrans$femaledatatransDt10G, na.rm = TRUE)
std <- sd(femaleRtrans$femaledatatransDt10G, na.rm = TRUE)
CVfemaleRtransDt10G <- CV(m, std)

m <- mean(femaleRtrans$femaledatatransPS, na.rm = TRUE)
std <- sd(femaleRtrans$femaledatatransPS, na.rm = TRUE)
CVfemaleRtransPS <- CV(m, std)

m <- mean(femaleRtrans$femaledatatransLL, na.rm = TRUE)
std <- sd(femaleRtrans$femaledatatransLL, na.rm = TRUE)
CVfemaleRtransLL <- CV(m, std)


#femaleStrans
femaleStrans <- subset(femaleTransALL, femalePopulation %in% "S")

CV <- function(m, std){(std/m)}
m <- mean(femaleStrans$femaledatatransDtG, na.rm = TRUE)
std <- sd(femaleStrans$femaledatatransDtG, na.rm = TRUE)
CVfemaleStransDtG <- CV(m, std)

m <- mean(femaleStrans$femaledatatransDt10G, na.rm = TRUE)
std <- sd(femaleStrans$femaledatatransDt10G, na.rm = TRUE)
CVfemaleStransDt10G <- CV(m, std)

m <- mean(femaleStrans$femaledatatransPS, na.rm = TRUE)
std <- sd(femaleStrans$femaledatatransPS, na.rm = TRUE)
CVfemaleStransPS <- CV(m, std)

m <- mean(femaleStrans$femaledatatransLL, na.rm = TRUE)
std <- sd(femaleStrans$femaledatatransLL, na.rm = TRUE)
CVfemaleStransLL <- CV(m, std)


#femaleWtrans
femaleWtrans <- subset(femaleTransALL, femalePopulation %in% "W")

CV <- function(m, std){(std/m)}
m <- mean(femaleWtrans$femaledatatransDtG, na.rm = TRUE)
std <- sd(femaleWtrans$femaledatatransDtG, na.rm = TRUE)
CVfemaleWtransDtG <- CV(m, std)

m <- mean(femaleWtrans$femaledatatransDt10G, na.rm = TRUE)
std <- sd(femaleWtrans$femaledatatransDt10G, na.rm = TRUE)
CVfemaleWtransDt10G <- CV(m, std)

m <- mean(femaleWtrans$femaledatatransPS, na.rm = TRUE)
std <- sd(femaleWtrans$femaledatatransPS, na.rm = TRUE)
CVfemaleWtransPS <- CV(m, std)

m <- mean(femaleWtrans$femaledatatransLL, na.rm = TRUE)
std <- sd(femaleWtrans$femaledatatransLL, na.rm = TRUE)
CVfemaleWtransLL <- CV(m, std)

######### writing outputs

CVDtGtransM <- as.data.frame(c(CVmaleAtransDtG,CVmaleDtransDtG,CVmaleEtransDtG,
                              CVmaleNtransDtG,CVmaleRtransDtG,
                              CVmaleStransDtG,CVmaleWtransDtG))

CVDt10GtransM <- as.data.frame(c(CVmaleAtransDt10G,CVmaleDtransDt10G,CVmaleEtransDt10G,
                                CVmaleNtransDt10G,CVmaleRtransDt10G,
                                CVmaleStransDt10G,CVmaleWtransDt10G))

CVPStransM <- as.data.frame(c(CVmaleAtransPS,CVmaleDtransPS,CVmaleEtransPS,
                             CVmaleNtransPS,CVmaleRtransPS,
                             CVmaleStransPS,CVmaleWtransPS))

CVLLtransM <- as.data.frame(c(CVmaleAtransLL,CVmaleDtransLL,CVmaleEtransLL,
                             CVmaleNtransLL,CVmaleRtransLL,
                             CVmaleStransLL,CVmaleWtransLL))

CVDtPtransM <- as.data.frame(c(CVmaleAtransDtP,CVmaleDtransDtP,CVmaleEtransDtP,
                              CVmaleNtransDtP,CVmaleRtransDtP,
                              CVmaleStransDtP,CVmaleWtransDtP))

CVNPtransM <- as.data.frame(c(CVmaleAtransNP,CVmaleDtransNP,CVmaleEtransNP,
                              CVmaleNtransNP,CVmaleRtransNP,
                              CVmaleStransNP,CVmaleWtransNP))

colnames(CVDtGtransM) <- c("DtG")
colnames(CVDt10GtransM) <- c("Dt10G")
colnames(CVPStransM) <- c("PS")
colnames(CVLLtransM) <- c("LL")
colnames(CVDtPtransM) <- c("DtP")
colnames(CVNPtransM) <- c("NP")

CVDtGtransF <- as.data.frame(c(CVfemaleAtransDtG,CVfemaleDtransDtG,CVfemaleEtransDtG,
                               CVfemaleNtransDtG,CVfemaleRtransDtG,
                               CVfemaleStransDtG,CVfemaleWtransDtG))

CVDt10GtransF <- as.data.frame(c(CVfemaleAtransDt10G,CVfemaleDtransDt10G,CVfemaleEtransDt10G,
                                 CVfemaleNtransDt10G,CVfemaleRtransDt10G,
                                 CVfemaleStransDt10G,CVfemaleWtransDt10G))

CVPStransF <- as.data.frame(c(CVfemaleAtransPS,CVfemaleDtransPS,CVfemaleEtransPS,
                              CVfemaleNtransPS,CVfemaleRtransPS,
                              CVfemaleStransPS,CVfemaleWtransPS))

CVLLtransF <- as.data.frame(c(CVfemaleAtransLL,CVfemaleDtransLL,CVfemaleEtransLL,
                              CVfemaleNtransLL,CVfemaleRtransLL,
                              CVfemaleStransLL,CVfemaleWtransLL))

colnames(CVDtGtransF) <- c("DtG")
colnames(CVDt10GtransF) <- c("Dt10G")
colnames(CVPStransF) <- c("PS")
colnames(CVLLtransF) <- c("LL")

colPops <- c("A","D","E","N","R","S","W")

coeffOfVarM <- as.data.frame(cbind(colPops,CVDtGtransM,CVDt10GtransM,CVPStransM,CVLLtransM,CVDtPtransM,CVNPtransM))
coeffOfVarF <- as.data.frame(cbind(colPops,CVDtGtransF,CVDt10GtransF,CVPStransF,CVLLtransF))


write.csv(coeffOfVarM, "tables/coeffOfVarM.csv")
write.csv(coeffOfVarF, "tables/coeffOfVarF.csv")





####################################################################################

#####this code is exclusively for making the output for use in the beta dispersion
## makings columns for output by trait, in alphabetical order by population, then by M then F 


#DtG
CVDtGtrans <- as.data.frame(c(CVmaleAtransDtG,CVfemaleAtransDtG,CVmaleDtransDtG,CVfemaleDtransDtG,CVmaleEtransDtG,
                              CVfemaleEtransDtG,CVmaleNtransDtG,CVfemaleNtransDtG,CVmaleRtransDtG,CVfemaleRtransDtG,
                              CVmaleStransDtG,CVfemaleStransDtG,CVmaleWtransDtG,CVfemaleWtransDtG))

CVDt10Gtrans <- as.data.frame(c(CVmaleAtransDt10G,CVfemaleAtransDt10G,CVmaleDtransDt10G,CVfemaleDtransDt10G,CVmaleEtransDt10G,
                                CVfemaleEtransDt10G,CVmaleNtransDt10G,CVfemaleNtransDt10G,CVmaleRtransDt10G,CVfemaleRtransDt10G,
                                CVmaleStransDt10G,CVfemaleStransDt10G,CVmaleWtransDt10G,CVfemaleWtransDt10G))

CVPStrans <- as.data.frame(c(CVmaleAtransPS,CVfemaleAtransPS,CVmaleDtransPS,CVfemaleDtransPS,CVmaleEtransPS,
                             CVfemaleEtransPS,CVmaleNtransPS,CVfemaleNtransPS,CVmaleRtransPS,CVfemaleRtransPS,
                             CVmaleStransPS,CVfemaleStransPS,CVmaleWtransPS,CVfemaleWtransPS))

CVLLtrans <- as.data.frame(c(CVmaleAtransLL,CVfemaleAtransLL,CVmaleDtransLL,CVfemaleDtransLL,CVmaleEtransLL,
                             CVfemaleEtransLL,CVmaleNtransLL,CVfemaleNtransLL,CVmaleRtransLL,CVfemaleRtransLL,
                             CVmaleStransLL,CVfemaleStransLL,CVmaleWtransLL,CVfemaleWtransLL))

colnames(CVDtGtrans) <- c("DtG")
colnames(CVDt10Gtrans) <- c("Dt10G")
colnames(CVPStrans) <- c("PS")
colnames(CVLLtrans) <- c("LL")

dataforbetatrans <- as.data.frame(cbind(CVDtGtrans,CVDt10Gtrans,CVPStrans,CVLLtrans))

write.csv(dataforbetatrans, "tables/dataforbetatrans.csv")






