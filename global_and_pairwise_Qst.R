
## Qst calculations

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


colnames(maleTransALL) <- c("DtG","Dt10G","PS","LL","DtP","NP","Population","Individual")
colnames(femaleTransALL) <- c("DtG","Dt10G","PS","LL","Population","Individual")


###################################################
### global Qst 

## females

## manova for Qst
FEMALEmanova <- manova(cbind(PS, DtG, Dt10G, LL) ~
                         Population + Population/Individual, 
                       data = femaleTransALL)

# pull out summary information containing sum of squares
sum_aov_female_global <- summary.aov(FEMALEmanova)

# extracting each trait's sum of squares
SSDtGFG <- as.data.frame(sum_aov_female_global$` Response DtG`$`Sum Sq`)
SSDt10GFG <- as.data.frame(sum_aov_female_global$` Response Dt10G`$`Sum Sq`)
SSPSFG <- as.data.frame(sum_aov_female_global$` Response PS`$`Sum Sq`)
SSLLFG <- as.data.frame(sum_aov_female_global$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGFG <- SSDtGFG[1,1]/(SSDtGFG[2,1]+SSDtGFG[1,1])
Dt10GFG <- SSDt10GFG[1,1]/(SSDt10GFG[2,1]+SSDt10GFG[1,1])
PSFG <- SSPSFG[1,1]/(SSPSFG[2,1]+SSPSFG[1,1])
LLFG <-  SSLLFG[1,1]/(SSLLFG[2,1]+SSLLFG[1,1])


QstFemaleGlobal <- as.data.frame(c(DtGFG, Dt10GFG, PSFG, LLFG), c("DtG",
                                                          "Dt10G", "PS", "LL"))
write.csv(QstFemaleGlobal, "tables/QstFemaleGlobal.csv")


FemaleGlobalBetweenSS <- as.data.frame(c(SSDtGFG[1,1], SSDt10GFG[1,1], SSPSFG[1,1], SSLLFG[1,1]), c("DtG",
                                                                  "Dt10G", "PS", "LL"))
write.csv(FemaleGlobalBetweenSS, "tables/FemaleGlobalBetweenSS.csv")


FemaleGlobalWithinSS <- as.data.frame(c(SSDtGFG[2,1], SSDt10GFG[2,1], SSPSFG[2,1], SSLLFG[2,1]), c("DtG",
                                                                                                    "Dt10G", "PS", "LL"))
write.csv(FemaleGlobalWithinSS, "tables/FemaleGlobalWithinSS.csv")


########################################################

##### males

MALEmanova <- manova(cbind(PS, DtG, Dt10G, LL, DtP, NP) ~
                       Population + Population/Individual, 
                     data = maleTransALL)

#pull out summary information containing sum of squares
sum_aov_male_global <- summary.aov(MALEmanova)

#extracting each trait's sum of squares
SSDtGMG <- as.data.frame(sum_aov_male_global$` Response DtG`$`Sum Sq`)
SSDt10GMG <- as.data.frame(sum_aov_male_global$` Response Dt10G`$`Sum Sq`)
SSPSMG <- as.data.frame(sum_aov_male_global$` Response PS`$`Sum Sq`)
SSLLMG <- as.data.frame(sum_aov_male_global$` Response LL`$`Sum Sq`)
SSDtPMG <- as.data.frame(sum_aov_male_global$` Response DtP`$`Sum Sq`)
SSNPMG <- as.data.frame(sum_aov_male_global$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGMG <- SSDtGMG[1,1]/(SSDtGMG[2,1]+SSDtGMG[1,1])
Dt10GMG <- SSDt10GMG[1,1]/(SSDt10GMG[2,1]+SSDt10GMG[1,1])
PSMG <- SSPSMG[1,1]/(SSPSMG[2,1]+SSPSMG[1,1])
LLMG <-  SSLLMG[1,1]/(SSLLMG[2,1]+SSLLMG[1,1])
DtPMG <-  SSDtPMG[1,1]/(SSDtPMG[2,1]+SSDtPMG[1,1])
NPMG <-  SSNPMG[1,1]/(SSNPMG[2,1]+SSNPMG[1,1])

QstMaleGlobal <- as.data.frame(c(DtGMG, Dt10GMG, PSMG, LLMG, DtPMG, NPMG), c("DtG",
                                                          "Dt10G", "PS", "LL", 
                                                          "DtP", "NP"))
write.csv(QstMaleGlobal, "tables/QstMaleGlobal.csv")


MaleGlobalBetweenSS <- as.data.frame(c(SSDtGMG[1,1], SSDt10GMG[1,1], SSPSMG[1,1], SSLLMG[1,1], SSDtPMG[1,1], SSNPMG[1,1]),
                                     c("DtG","Dt10G", "PS", "LL","DtP","NP"))

write.csv(MaleGlobalBetweenSS, "tables/MaleGlobalBetweenSS.csv")


MaleGlobalWithinSS <- as.data.frame(c(SSDtGMG[2,1], SSDt10GMG[2,1], SSPSMG[2,1], SSLLMG[2,1], SSDtPMG[2,1], SSNPMG[2,1]), 
                                    c("DtG","Dt10G", "PS", "LL","DtP","NP"))

write.csv(MaleGlobalWithinSS, "tables/MaleGlobalWithinSS.csv")


########################################################

### pairwise Qst

##females

FemalemanovaAD <- manova(cbind(PS, DtG, Dt10G, LL) ~
                         Population + Population/Individual, 
                       data = femaleTransALL,
                       subset = Population %in% c("A", "D"))


sum_aov_female_AD <- summary.aov(FemalemanovaAD)

#extracting each trait's sum of squares
SSDtGfemale_AD <- as.data.frame(sum_aov_female_AD$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_AD <- as.data.frame(sum_aov_female_AD$` Response Dt10G`$`Sum Sq`)
SSPSfemale_AD <- as.data.frame(sum_aov_female_AD$` Response PS`$`Sum Sq`)
SSLLfemale_AD<- as.data.frame(sum_aov_female_AD$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_AD <- SSDtGfemale_AD[1,1]/(SSDtGfemale_AD[2,1]+SSDtGfemale_AD[1,1])
Dt10Gfemale_AD <- SSDt10Gfemale_AD[1,1]/(SSDt10Gfemale_AD[2,1]+SSDt10Gfemale_AD[1,1])
PSfemale_AD <- SSPSfemale_AD[1,1]/(SSPSfemale_AD[2,1]+SSPSfemale_AD[1,1])
LLfemale_AD <-  SSLLfemale_AD[1,1]/(SSLLfemale_AD[2,1]+SSLLfemale_AD[1,1])


QstFemaleAD <- as.data.frame(c(DtGfemale_AD, Dt10Gfemale_AD, PSfemale_AD, LLfemale_AD), c("DtG",
                                                                  "Dt10G", "PS", "LL"))


#####
FemalemanovaAE <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("A", "E"))


sum_aov_female_AE <- summary.aov(FemalemanovaAE)

#extracting each trait's sum of squares
SSDtGfemale_AE <- as.data.frame(sum_aov_female_AE$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_AE <- as.data.frame(sum_aov_female_AE$` Response Dt10G`$`Sum Sq`)
SSPSfemale_AE <- as.data.frame(sum_aov_female_AE$` Response PS`$`Sum Sq`)
SSLLfemale_AE<- as.data.frame(sum_aov_female_AE$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_AE <- SSDtGfemale_AE[1,1]/(SSDtGfemale_AE[2,1]+SSDtGfemale_AE[1,1])
Dt10Gfemale_AE <- SSDt10Gfemale_AE[1,1]/(SSDt10Gfemale_AE[2,1]+SSDt10Gfemale_AE[1,1])
PSfemale_AE <- SSPSfemale_AE[1,1]/(SSPSfemale_AE[2,1]+SSPSfemale_AE[1,1])
LLfemale_AE <-  SSLLfemale_AE[1,1]/(SSLLfemale_AE[2,1]+SSLLfemale_AE[1,1])


QstFemaleAE <- as.data.frame(c(DtGfemale_AE, Dt10Gfemale_AE, PSfemale_AE, LLfemale_AE), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaAN <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("A", "N"))


sum_aov_female_AN <- summary.aov(FemalemanovaAN)

#extracting each trait's sum of squares
SSDtGfemale_AN <- as.data.frame(sum_aov_female_AN$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_AN <- as.data.frame(sum_aov_female_AN$` Response Dt10G`$`Sum Sq`)
SSPSfemale_AN <- as.data.frame(sum_aov_female_AN$` Response PS`$`Sum Sq`)
SSLLfemale_AN<- as.data.frame(sum_aov_female_AN$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_AN <- SSDtGfemale_AN[1,1]/(SSDtGfemale_AN[2,1]+SSDtGfemale_AN[1,1])
Dt10Gfemale_AN <- SSDt10Gfemale_AN[1,1]/(SSDt10Gfemale_AN[2,1]+SSDt10Gfemale_AN[1,1])
PSfemale_AN <- SSPSfemale_AN[1,1]/(SSPSfemale_AN[2,1]+SSPSfemale_AN[1,1])
LLfemale_AN <-  SSLLfemale_AN[1,1]/(SSLLfemale_AN[2,1]+SSLLfemale_AN[1,1])


QstFemaleAN <- as.data.frame(c(DtGfemale_AN, Dt10Gfemale_AN, PSfemale_AN, LLfemale_AN), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaAR <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("A", "R"))


sum_aov_female_AR <- summary.aov(FemalemanovaAR)

#extracting each trait's sum of squares
SSDtGfemale_AR <- as.data.frame(sum_aov_female_AR$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_AR <- as.data.frame(sum_aov_female_AR$` Response Dt10G`$`Sum Sq`)
SSPSfemale_AR <- as.data.frame(sum_aov_female_AR$` Response PS`$`Sum Sq`)
SSLLfemale_AR<- as.data.frame(sum_aov_female_AR$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_AR <- SSDtGfemale_AR[1,1]/(SSDtGfemale_AR[2,1]+SSDtGfemale_AR[1,1])
Dt10Gfemale_AR <- SSDt10Gfemale_AR[1,1]/(SSDt10Gfemale_AR[2,1]+SSDt10Gfemale_AR[1,1])
PSfemale_AR <- SSPSfemale_AR[1,1]/(SSPSfemale_AR[2,1]+SSPSfemale_AR[1,1])
LLfemale_AR <-  SSLLfemale_AR[1,1]/(SSLLfemale_AR[2,1]+SSLLfemale_AR[1,1])


QstFemaleAR <- as.data.frame(c(DtGfemale_AR, Dt10Gfemale_AR, PSfemale_AR, LLfemale_AR), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaAS <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("A", "S"))


sum_aov_female_AS <- summary.aov(FemalemanovaAS)

#extracting each trait's sum of squares
SSDtGfemale_AS <- as.data.frame(sum_aov_female_AS$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_AS <- as.data.frame(sum_aov_female_AS$` Response Dt10G`$`Sum Sq`)
SSPSfemale_AS <- as.data.frame(sum_aov_female_AS$` Response PS`$`Sum Sq`)
SSLLfemale_AS<- as.data.frame(sum_aov_female_AS$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_AS <- SSDtGfemale_AS[1,1]/(SSDtGfemale_AS[2,1]+SSDtGfemale_AS[1,1])
Dt10Gfemale_AS <- SSDt10Gfemale_AS[1,1]/(SSDt10Gfemale_AS[2,1]+SSDt10Gfemale_AS[1,1])
PSfemale_AS <- SSPSfemale_AS[1,1]/(SSPSfemale_AS[2,1]+SSPSfemale_AS[1,1])
LLfemale_AS <-  SSLLfemale_AS[1,1]/(SSLLfemale_AS[2,1]+SSLLfemale_AS[1,1])


QstFemaleAS <- as.data.frame(c(DtGfemale_AS, Dt10Gfemale_AS, PSfemale_AS, LLfemale_AS), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaAW <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("A", "W"))


sum_aov_female_AW <- summary.aov(FemalemanovaAW)

#extracting each trait's sum of squares
SSDtGfemale_AW <- as.data.frame(sum_aov_female_AW$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_AW <- as.data.frame(sum_aov_female_AW$` Response Dt10G`$`Sum Sq`)
SSPSfemale_AW <- as.data.frame(sum_aov_female_AW$` Response PS`$`Sum Sq`)
SSLLfemale_AW<- as.data.frame(sum_aov_female_AW$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_AW <- SSDtGfemale_AW[1,1]/(SSDtGfemale_AW[2,1]+SSDtGfemale_AW[1,1])
Dt10Gfemale_AW <- SSDt10Gfemale_AW[1,1]/(SSDt10Gfemale_AW[2,1]+SSDt10Gfemale_AW[1,1])
PSfemale_AW <- SSPSfemale_AW[1,1]/(SSPSfemale_AW[2,1]+SSPSfemale_AW[1,1])
LLfemale_AW <-  SSLLfemale_AW[1,1]/(SSLLfemale_AW[2,1]+SSLLfemale_AW[1,1])


QstFemaleAW <- as.data.frame(c(DtGfemale_AW, Dt10Gfemale_AW, PSfemale_AW, LLfemale_AW), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaDE <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("D", "E"))


sum_aov_female_DE <- summary.aov(FemalemanovaDE)

#extracting each trait's sum of squares
SSDtGfemale_DE <- as.data.frame(sum_aov_female_DE$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_DE <- as.data.frame(sum_aov_female_DE$` Response Dt10G`$`Sum Sq`)
SSPSfemale_DE <- as.data.frame(sum_aov_female_DE$` Response PS`$`Sum Sq`)
SSLLfemale_DE<- as.data.frame(sum_aov_female_DE$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_DE <- SSDtGfemale_DE[1,1]/(SSDtGfemale_DE[2,1]+SSDtGfemale_DE[1,1])
Dt10Gfemale_DE <- SSDt10Gfemale_DE[1,1]/(SSDt10Gfemale_DE[2,1]+SSDt10Gfemale_DE[1,1])
PSfemale_DE <- SSPSfemale_DE[1,1]/(SSPSfemale_DE[2,1]+SSPSfemale_DE[1,1])
LLfemale_DE <-  SSLLfemale_DE[1,1]/(SSLLfemale_DE[2,1]+SSLLfemale_DE[1,1])


QstFemaleDE <- as.data.frame(c(DtGfemale_DE, Dt10Gfemale_DE, PSfemale_DE, LLfemale_DE), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaDN <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("D", "N"))


sum_aov_female_DN <- summary.aov(FemalemanovaDN)

#extracting each trait's sum of squares
SSDtGfemale_DN <- as.data.frame(sum_aov_female_DN$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_DN <- as.data.frame(sum_aov_female_DN$` Response Dt10G`$`Sum Sq`)
SSPSfemale_DN <- as.data.frame(sum_aov_female_DN$` Response PS`$`Sum Sq`)
SSLLfemale_DN<- as.data.frame(sum_aov_female_DN$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_DN <- SSDtGfemale_DN[1,1]/(SSDtGfemale_DN[2,1]+SSDtGfemale_DN[1,1])
Dt10Gfemale_DN <- SSDt10Gfemale_DN[1,1]/(SSDt10Gfemale_DN[2,1]+SSDt10Gfemale_DN[1,1])
PSfemale_DN <- SSPSfemale_DN[1,1]/(SSPSfemale_DN[2,1]+SSPSfemale_DN[1,1])
LLfemale_DN <-  SSLLfemale_DN[1,1]/(SSLLfemale_DN[2,1]+SSLLfemale_DN[1,1])


QstFemaleDN <- as.data.frame(c(DtGfemale_DN, Dt10Gfemale_DN, PSfemale_DN, LLfemale_DN), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaDR <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("D", "R"))


sum_aov_female_DR <- summary.aov(FemalemanovaDR)

#extracting each trait's sum of squares
SSDtGfemale_DR <- as.data.frame(sum_aov_female_DR$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_DR <- as.data.frame(sum_aov_female_DR$` Response Dt10G`$`Sum Sq`)
SSPSfemale_DR <- as.data.frame(sum_aov_female_DR$` Response PS`$`Sum Sq`)
SSLLfemale_DR<- as.data.frame(sum_aov_female_DR$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_DR <- SSDtGfemale_DR[1,1]/(SSDtGfemale_DR[2,1]+SSDtGfemale_DR[1,1])
Dt10Gfemale_DR <- SSDt10Gfemale_DR[1,1]/(SSDt10Gfemale_DR[2,1]+SSDt10Gfemale_DR[1,1])
PSfemale_DR <- SSPSfemale_DR[1,1]/(SSPSfemale_DR[2,1]+SSPSfemale_DR[1,1])
LLfemale_DR <-  SSLLfemale_DR[1,1]/(SSLLfemale_DR[2,1]+SSLLfemale_DR[1,1])


QstFemaleDR <- as.data.frame(c(DtGfemale_DR, Dt10Gfemale_DR, PSfemale_DR, LLfemale_DR), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaDS <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("D", "S"))


sum_aov_female_DS <- summary.aov(FemalemanovaDS)

#extracting each trait's sum of squares
SSDtGfemale_DS <- as.data.frame(sum_aov_female_DS$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_DS <- as.data.frame(sum_aov_female_DS$` Response Dt10G`$`Sum Sq`)
SSPSfemale_DS <- as.data.frame(sum_aov_female_DS$` Response PS`$`Sum Sq`)
SSLLfemale_DS<- as.data.frame(sum_aov_female_DS$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_DS <- SSDtGfemale_DS[1,1]/(SSDtGfemale_DS[2,1]+SSDtGfemale_DS[1,1])
Dt10Gfemale_DS <- SSDt10Gfemale_DS[1,1]/(SSDt10Gfemale_DS[2,1]+SSDt10Gfemale_DS[1,1])
PSfemale_DS <- SSPSfemale_DS[1,1]/(SSPSfemale_DS[2,1]+SSPSfemale_DS[1,1])
LLfemale_DS <-  SSLLfemale_DS[1,1]/(SSLLfemale_DS[2,1]+SSLLfemale_DS[1,1])


QstFemaleDS <- as.data.frame(c(DtGfemale_DS, Dt10Gfemale_DS, PSfemale_DS, LLfemale_DS), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaDW <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("D", "W"))


sum_aov_female_DW <- summary.aov(FemalemanovaDW)

#extracting each trait's sum of squares
SSDtGfemale_DW <- as.data.frame(sum_aov_female_DW$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_DW <- as.data.frame(sum_aov_female_DW$` Response Dt10G`$`Sum Sq`)
SSPSfemale_DW <- as.data.frame(sum_aov_female_DW$` Response PS`$`Sum Sq`)
SSLLfemale_DW<- as.data.frame(sum_aov_female_DW$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_DW <- SSDtGfemale_DW[1,1]/(SSDtGfemale_DW[2,1]+SSDtGfemale_DW[1,1])
Dt10Gfemale_DW <- SSDt10Gfemale_DW[1,1]/(SSDt10Gfemale_DW[2,1]+SSDt10Gfemale_DW[1,1])
PSfemale_DW <- SSPSfemale_DW[1,1]/(SSPSfemale_DW[2,1]+SSPSfemale_DW[1,1])
LLfemale_DW <-  SSLLfemale_DW[1,1]/(SSLLfemale_DW[2,1]+SSLLfemale_DW[1,1])


QstFemaleDW <- as.data.frame(c(DtGfemale_DW, Dt10Gfemale_DW, PSfemale_DW, LLfemale_DW), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))



#####
FemalemanovaEN <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("E", "N"))


sum_aov_female_EN <- summary.aov(FemalemanovaEN)

#extracting each trait's sum of squares
SSDtGfemale_EN <- as.data.frame(sum_aov_female_EN$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_EN <- as.data.frame(sum_aov_female_EN$` Response Dt10G`$`Sum Sq`)
SSPSfemale_EN <- as.data.frame(sum_aov_female_EN$` Response PS`$`Sum Sq`)
SSLLfemale_EN<- as.data.frame(sum_aov_female_EN$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_EN <- SSDtGfemale_EN[1,1]/(SSDtGfemale_EN[2,1]+SSDtGfemale_EN[1,1])
Dt10Gfemale_EN <- SSDt10Gfemale_EN[1,1]/(SSDt10Gfemale_EN[2,1]+SSDt10Gfemale_EN[1,1])
PSfemale_EN <- SSPSfemale_EN[1,1]/(SSPSfemale_EN[2,1]+SSPSfemale_EN[1,1])
LLfemale_EN <-  SSLLfemale_EN[1,1]/(SSLLfemale_EN[2,1]+SSLLfemale_EN[1,1])


QstFemaleEN <- as.data.frame(c(DtGfemale_EN, Dt10Gfemale_EN, PSfemale_EN, LLfemale_EN), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaER <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("E", "R"))


sum_aov_female_ER <- summary.aov(FemalemanovaER)

#extracting each trait's sum of squares
SSDtGfemale_ER <- as.data.frame(sum_aov_female_ER$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_ER <- as.data.frame(sum_aov_female_ER$` Response Dt10G`$`Sum Sq`)
SSPSfemale_ER <- as.data.frame(sum_aov_female_ER$` Response PS`$`Sum Sq`)
SSLLfemale_ER<- as.data.frame(sum_aov_female_ER$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_ER <- SSDtGfemale_ER[1,1]/(SSDtGfemale_ER[2,1]+SSDtGfemale_ER[1,1])
Dt10Gfemale_ER <- SSDt10Gfemale_ER[1,1]/(SSDt10Gfemale_ER[2,1]+SSDt10Gfemale_ER[1,1])
PSfemale_ER <- SSPSfemale_ER[1,1]/(SSPSfemale_ER[2,1]+SSPSfemale_ER[1,1])
LLfemale_ER <-  SSLLfemale_ER[1,1]/(SSLLfemale_ER[2,1]+SSLLfemale_ER[1,1])


QstFemaleER <- as.data.frame(c(DtGfemale_ER, Dt10Gfemale_ER, PSfemale_ER, LLfemale_ER), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaES <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("E", "S"))


sum_aov_female_ES <- summary.aov(FemalemanovaES)

#extracting each trait's sum of squares
SSDtGfemale_ES <- as.data.frame(sum_aov_female_ES$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_ES <- as.data.frame(sum_aov_female_ES$` Response Dt10G`$`Sum Sq`)
SSPSfemale_ES <- as.data.frame(sum_aov_female_ES$` Response PS`$`Sum Sq`)
SSLLfemale_ES<- as.data.frame(sum_aov_female_ES$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_ES <- SSDtGfemale_ES[1,1]/(SSDtGfemale_ES[2,1]+SSDtGfemale_ES[1,1])
Dt10Gfemale_ES <- SSDt10Gfemale_ES[1,1]/(SSDt10Gfemale_ES[2,1]+SSDt10Gfemale_ES[1,1])
PSfemale_ES <- SSPSfemale_ES[1,1]/(SSPSfemale_ES[2,1]+SSPSfemale_ES[1,1])
LLfemale_ES <-  SSLLfemale_ES[1,1]/(SSLLfemale_ES[2,1]+SSLLfemale_ES[1,1])


QstFemaleES <- as.data.frame(c(DtGfemale_ES, Dt10Gfemale_ES, PSfemale_ES, LLfemale_ES), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaEW <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("E", "W"))


sum_aov_female_EW <- summary.aov(FemalemanovaEW)

#extracting each trait's sum of squares
SSDtGfemale_EW <- as.data.frame(sum_aov_female_EW$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_EW <- as.data.frame(sum_aov_female_EW$` Response Dt10G`$`Sum Sq`)
SSPSfemale_EW <- as.data.frame(sum_aov_female_EW$` Response PS`$`Sum Sq`)
SSLLfemale_EW<- as.data.frame(sum_aov_female_EW$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_EW <- SSDtGfemale_EW[1,1]/(SSDtGfemale_EW[2,1]+SSDtGfemale_EW[1,1])
Dt10Gfemale_EW <- SSDt10Gfemale_EW[1,1]/(SSDt10Gfemale_EW[2,1]+SSDt10Gfemale_EW[1,1])
PSfemale_EW <- SSPSfemale_EW[1,1]/(SSPSfemale_EW[2,1]+SSPSfemale_EW[1,1])
LLfemale_EW <-  SSLLfemale_EW[1,1]/(SSLLfemale_EW[2,1]+SSLLfemale_EW[1,1])


QstFemaleEW <- as.data.frame(c(DtGfemale_EW, Dt10Gfemale_EW, PSfemale_EW, LLfemale_EW), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))



#####
FemalemanovaNR <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("N", "R"))


sum_aov_female_NR <- summary.aov(FemalemanovaNR)

#extracting each trait's sum of squares
SSDtGfemale_NR <- as.data.frame(sum_aov_female_NR$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_NR <- as.data.frame(sum_aov_female_NR$` Response Dt10G`$`Sum Sq`)
SSPSfemale_NR <- as.data.frame(sum_aov_female_NR$` Response PS`$`Sum Sq`)
SSLLfemale_NR<- as.data.frame(sum_aov_female_NR$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_NR <- SSDtGfemale_NR[1,1]/(SSDtGfemale_NR[2,1]+SSDtGfemale_NR[1,1])
Dt10Gfemale_NR <- SSDt10Gfemale_NR[1,1]/(SSDt10Gfemale_NR[2,1]+SSDt10Gfemale_NR[1,1])
PSfemale_NR <- SSPSfemale_NR[1,1]/(SSPSfemale_NR[2,1]+SSPSfemale_NR[1,1])
LLfemale_NR <-  SSLLfemale_NR[1,1]/(SSLLfemale_NR[2,1]+SSLLfemale_NR[1,1])


QstFemaleNR <- as.data.frame(c(DtGfemale_NR, Dt10Gfemale_NR, PSfemale_NR, LLfemale_NR), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaNS <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("N", "S"))


sum_aov_female_NS <- summary.aov(FemalemanovaNS)

#extracting each trait's sum of squares
SSDtGfemale_NS <- as.data.frame(sum_aov_female_NS$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_NS <- as.data.frame(sum_aov_female_NS$` Response Dt10G`$`Sum Sq`)
SSPSfemale_NS <- as.data.frame(sum_aov_female_NS$` Response PS`$`Sum Sq`)
SSLLfemale_NS<- as.data.frame(sum_aov_female_NS$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_NS <- SSDtGfemale_NS[1,1]/(SSDtGfemale_NS[2,1]+SSDtGfemale_NS[1,1])
Dt10Gfemale_NS <- SSDt10Gfemale_NS[1,1]/(SSDt10Gfemale_NS[2,1]+SSDt10Gfemale_NS[1,1])
PSfemale_NS <- SSPSfemale_NS[1,1]/(SSPSfemale_NS[2,1]+SSPSfemale_NS[1,1])
LLfemale_NS <-  SSLLfemale_NS[1,1]/(SSLLfemale_NS[2,1]+SSLLfemale_NS[1,1])


QstFemaleNS <- as.data.frame(c(DtGfemale_NS, Dt10Gfemale_NS, PSfemale_NS, LLfemale_NS), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaNW <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("N", "W"))


sum_aov_female_NW <- summary.aov(FemalemanovaNW)

#extracting each trait's sum of squares
SSDtGfemale_NW <- as.data.frame(sum_aov_female_NW$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_NW <- as.data.frame(sum_aov_female_NW$` Response Dt10G`$`Sum Sq`)
SSPSfemale_NW <- as.data.frame(sum_aov_female_NW$` Response PS`$`Sum Sq`)
SSLLfemale_NW<- as.data.frame(sum_aov_female_NW$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_NW <- SSDtGfemale_NW[1,1]/(SSDtGfemale_NW[2,1]+SSDtGfemale_NW[1,1])
Dt10Gfemale_NW <- SSDt10Gfemale_NW[1,1]/(SSDt10Gfemale_NW[2,1]+SSDt10Gfemale_NW[1,1])
PSfemale_NW <- SSPSfemale_NW[1,1]/(SSPSfemale_NW[2,1]+SSPSfemale_NW[1,1])
LLfemale_NW <-  SSLLfemale_NW[1,1]/(SSLLfemale_NW[2,1]+SSLLfemale_NW[1,1])


QstFemaleNW <- as.data.frame(c(DtGfemale_NW, Dt10Gfemale_NW, PSfemale_NW, LLfemale_NW), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))


#####
FemalemanovaRS <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("R", "S"))


sum_aov_female_RS <- summary.aov(FemalemanovaNR)

#extracting each trait's sum of squares
SSDtGfemale_RS <- as.data.frame(sum_aov_female_RS$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_RS <- as.data.frame(sum_aov_female_RS$` Response Dt10G`$`Sum Sq`)
SSPSfemale_RS <- as.data.frame(sum_aov_female_RS$` Response PS`$`Sum Sq`)
SSLLfemale_RS<- as.data.frame(sum_aov_female_RS$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_RS <- SSDtGfemale_RS[1,1]/(SSDtGfemale_RS[2,1]+SSDtGfemale_RS[1,1])
Dt10Gfemale_RS <- SSDt10Gfemale_RS[1,1]/(SSDt10Gfemale_RS[2,1]+SSDt10Gfemale_RS[1,1])
PSfemale_RS <- SSPSfemale_RS[1,1]/(SSPSfemale_RS[2,1]+SSPSfemale_RS[1,1])
LLfemale_RS <-  SSLLfemale_RS[1,1]/(SSLLfemale_RS[2,1]+SSLLfemale_RS[1,1])


QstFemaleRS <- as.data.frame(c(DtGfemale_RS, Dt10Gfemale_RS, PSfemale_RS, LLfemale_RS), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))



#####
FemalemanovaRW <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("R", "W"))


sum_aov_female_RW <- summary.aov(FemalemanovaRW)

#extracting each trait's sum of squares
SSDtGfemale_RW <- as.data.frame(sum_aov_female_RW$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_RW <- as.data.frame(sum_aov_female_RW$` Response Dt10G`$`Sum Sq`)
SSPSfemale_RW <- as.data.frame(sum_aov_female_RW$` Response PS`$`Sum Sq`)
SSLLfemale_RW<- as.data.frame(sum_aov_female_RW$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_RW <- SSDtGfemale_RW[1,1]/(SSDtGfemale_RW[2,1]+SSDtGfemale_RW[1,1])
Dt10Gfemale_RW <- SSDt10Gfemale_RW[1,1]/(SSDt10Gfemale_RW[2,1]+SSDt10Gfemale_RW[1,1])
PSfemale_RW <- SSPSfemale_RW[1,1]/(SSPSfemale_RW[2,1]+SSPSfemale_RW[1,1])
LLfemale_RW <-  SSLLfemale_RW[1,1]/(SSLLfemale_RW[2,1]+SSLLfemale_RW[1,1])


QstFemaleRW <- as.data.frame(c(DtGfemale_RW, Dt10Gfemale_RW, PSfemale_RW, LLfemale_RW), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

#####
FemalemanovaSW <- manova(cbind(PS, DtG, Dt10G, LL) ~
                           Population + Population/Individual, 
                         data = femaleTransALL,
                         subset = Population %in% c("S", "W"))


sum_aov_female_SW <- summary.aov(FemalemanovaSW)

#extracting each trait's sum of squares
SSDtGfemale_SW <- as.data.frame(sum_aov_female_SW$` Response DtG`$`Sum Sq`)
SSDt10Gfemale_SW <- as.data.frame(sum_aov_female_SW$` Response Dt10G`$`Sum Sq`)
SSPSfemale_SW <- as.data.frame(sum_aov_female_SW$` Response PS`$`Sum Sq`)
SSLLfemale_SW<- as.data.frame(sum_aov_female_SW$` Response LL`$`Sum Sq`)


#haploid Qst calculation
DtGfemale_SW <- SSDtGfemale_SW[1,1]/(SSDtGfemale_SW[2,1]+SSDtGfemale_SW[1,1])
Dt10Gfemale_SW <- SSDt10Gfemale_SW[1,1]/(SSDt10Gfemale_SW[2,1]+SSDt10Gfemale_SW[1,1])
PSfemale_SW <- SSPSfemale_SW[1,1]/(SSPSfemale_SW[2,1]+SSPSfemale_SW[1,1])
LLfemale_SW <-  SSLLfemale_SW[1,1]/(SSLLfemale_SW[2,1]+SSLLfemale_SW[1,1])


QstFemaleSW <- as.data.frame(c(DtGfemale_SW, Dt10Gfemale_SW, PSfemale_SW, LLfemale_SW), c("DtG",
                                                                                          "Dt10G", "PS", "LL"))

###########################################################################
#####writing female outputs


### DtG

femalePairwiseQstsDtGA <- as.data.frame(c(0,QstFemaleAD[1,],QstFemaleAE[1,],QstFemaleAN[1,],QstFemaleAR[1,],QstFemaleAS[1,],QstFemaleAW[1,]))
femalePairwiseQstsDtGD <- as.data.frame(c(QstFemaleAD[1,],0,QstFemaleDE[1,],QstFemaleDN[1,],QstFemaleDR[1,],QstFemaleDS[1,],QstFemaleDW[1,]))
femalePairwiseQstsDtGE <- as.data.frame(c(QstFemaleAE[1,],QstFemaleDE[1,],0,QstFemaleEN[1,],QstFemaleER[1,],QstFemaleES[1,],QstFemaleEW[1,]))
femalePairwiseQstsDtGN <- as.data.frame(c(QstFemaleAN[1,],QstFemaleDN[1,],QstFemaleEN[1,],0,QstFemaleNR[1,],QstFemaleNS[1,],QstFemaleNW[1,]))
femalePairwiseQstsDtGR <- as.data.frame(c(QstFemaleAR[1,],QstFemaleDR[1,],QstFemaleER[1,],QstFemaleNR[1,],0,QstFemaleRS[1,],QstFemaleRW[1,]))
femalePairwiseQstsDtGS <- as.data.frame(c(QstFemaleAS[1,],QstFemaleDS[1,],QstFemaleES[1,],QstFemaleNS[1,],QstFemaleRS[1,],0,QstFemaleSW[1,]))
femalePairwiseQstsDtGW <- as.data.frame(c(QstFemaleAW[1,],QstFemaleDW[1,],QstFemaleEW[1,],QstFemaleNW[1,],QstFemaleRW[1,],QstFemaleSW[1,],0))


femalePairwiseALLDtG <- as.data.frame(cbind(femalePairwiseQstsDtGA,femalePairwiseQstsDtGD,femalePairwiseQstsDtGE,
                                         femalePairwiseQstsDtGN,femalePairwiseQstsDtGR,femalePairwiseQstsDtGS,
                                         femalePairwiseQstsDtGW))
colnames(femalePairwiseALLDtG) <- c("A","D","E","R","N","S","W")
rownames(femalePairwiseALLDtG) <- c("A","D","E","R","N","S","W")


### Dt10G

femalePairwiseQstsDt10GA <- as.data.frame(c(0,QstFemaleAD[2,],QstFemaleAE[2,],QstFemaleAN[2,],QstFemaleAR[2,],QstFemaleAS[2,],QstFemaleAW[2,]))
femalePairwiseQstsDt10GD <- as.data.frame(c(QstFemaleAD[2,],0,QstFemaleDE[2,],QstFemaleDN[2,],QstFemaleDR[2,],QstFemaleDS[2,],QstFemaleDW[2,]))
femalePairwiseQstsDt10GE <- as.data.frame(c(QstFemaleAE[2,],QstFemaleDE[2,],0,QstFemaleEN[2,],QstFemaleER[2,],QstFemaleES[2,],QstFemaleEW[2,]))
femalePairwiseQstsDt10GN <- as.data.frame(c(QstFemaleAN[2,],QstFemaleDN[2,],QstFemaleEN[2,],0,QstFemaleNR[2,],QstFemaleNS[2,],QstFemaleNW[2,]))
femalePairwiseQstsDt10GR <- as.data.frame(c(QstFemaleAR[2,],QstFemaleDR[2,],QstFemaleER[2,],QstFemaleNR[2,],0,QstFemaleRS[2,],QstFemaleRW[2,]))
femalePairwiseQstsDt10GS <- as.data.frame(c(QstFemaleAS[2,],QstFemaleDS[2,],QstFemaleES[2,],QstFemaleNS[2,],QstFemaleRS[2,],0,QstFemaleSW[2,]))
femalePairwiseQstsDt10GW <- as.data.frame(c(QstFemaleAW[2,],QstFemaleDW[2,],QstFemaleEW[2,],QstFemaleNW[2,],QstFemaleRW[2,],QstFemaleSW[2,],0))


femalePairwiseALLDt10G <- as.data.frame(cbind(femalePairwiseQstsDt10GA,femalePairwiseQstsDt10GD,femalePairwiseQstsDt10GE,
                                            femalePairwiseQstsDt10GN,femalePairwiseQstsDt10GR,femalePairwiseQstsDt10GS,
                                            femalePairwiseQstsDt10GW))
colnames(femalePairwiseALLDt10G) <- c("A","D","E","R","N","S","W")
rownames(femalePairwiseALLDt10G) <- c("A","D","E","R","N","S","W")


### PS

femalePairwiseQstsPSA <- as.data.frame(c(0,QstFemaleAD[3,],QstFemaleAE[3,],QstFemaleAN[3,],QstFemaleAR[3,],QstFemaleAS[3,],QstFemaleAW[3,]))
femalePairwiseQstsPSD <- as.data.frame(c(QstFemaleAD[3,],0,QstFemaleDE[3,],QstFemaleDN[3,],QstFemaleDR[3,],QstFemaleDS[3,],QstFemaleDW[3,]))
femalePairwiseQstsPSE <- as.data.frame(c(QstFemaleAE[3,],QstFemaleDE[3,],0,QstFemaleEN[3,],QstFemaleER[3,],QstFemaleES[3,],QstFemaleEW[3,]))
femalePairwiseQstsPSN <- as.data.frame(c(QstFemaleAN[3,],QstFemaleDN[3,],QstFemaleEN[3,],0,QstFemaleNR[3,],QstFemaleNS[3,],QstFemaleNW[3,]))
femalePairwiseQstsPSR <- as.data.frame(c(QstFemaleAR[3,],QstFemaleDR[3,],QstFemaleER[3,],QstFemaleNR[3,],0,QstFemaleRS[3,],QstFemaleRW[3,]))
femalePairwiseQstsPSS <- as.data.frame(c(QstFemaleAS[3,],QstFemaleDS[3,],QstFemaleES[3,],QstFemaleNS[3,],QstFemaleRS[3,],0,QstFemaleSW[3,]))
femalePairwiseQstsPSW <- as.data.frame(c(QstFemaleAW[3,],QstFemaleDW[3,],QstFemaleEW[3,],QstFemaleNW[3,],QstFemaleRW[3,],QstFemaleSW[3,],0))


femalePairwiseALLPS <- as.data.frame(cbind(femalePairwiseQstsPSA,femalePairwiseQstsPSD,femalePairwiseQstsPSE,
                                            femalePairwiseQstsPSN,femalePairwiseQstsPSR,femalePairwiseQstsPSS,
                                            femalePairwiseQstsPSW))
colnames(femalePairwiseALLPS) <- c("A","D","E","R","N","S","W")
rownames(femalePairwiseALLPS) <- c("A","D","E","R","N","S","W")

### LL

femalePairwiseQstsLLA <- as.data.frame(c(0,QstFemaleAD[4,],QstFemaleAE[4,],QstFemaleAN[4,],QstFemaleAR[4,],QstFemaleAS[4,],QstFemaleAW[4,]))
femalePairwiseQstsLLD <- as.data.frame(c(QstFemaleAD[4,],0,QstFemaleDE[4,],QstFemaleDN[4,],QstFemaleDR[4,],QstFemaleDS[4,],QstFemaleDW[4,]))
femalePairwiseQstsLLE <- as.data.frame(c(QstFemaleAE[4,],QstFemaleDE[4,],0,QstFemaleEN[4,],QstFemaleER[4,],QstFemaleES[4,],QstFemaleEW[4,]))
femalePairwiseQstsLLN <- as.data.frame(c(QstFemaleAN[4,],QstFemaleDN[4,],QstFemaleEN[4,],0,QstFemaleNR[4,],QstFemaleNS[4,],QstFemaleNW[4,]))
femalePairwiseQstsLLR <- as.data.frame(c(QstFemaleAR[4,],QstFemaleDR[4,],QstFemaleER[4,],QstFemaleNR[4,],0,QstFemaleRS[4,],QstFemaleRW[4,]))
femalePairwiseQstsLLS <- as.data.frame(c(QstFemaleAS[4,],QstFemaleDS[4,],QstFemaleES[4,],QstFemaleNS[4,],QstFemaleRS[4,],0,QstFemaleSW[4,]))
femalePairwiseQstsLLW <- as.data.frame(c(QstFemaleAW[4,],QstFemaleDW[4,],QstFemaleEW[4,],QstFemaleNW[4,],QstFemaleRW[4,],QstFemaleSW[4,],0))


femalePairwiseALLLL <- as.data.frame(cbind(femalePairwiseQstsLLA,femalePairwiseQstsLLD,femalePairwiseQstsLLE,
                                            femalePairwiseQstsLLN,femalePairwiseQstsLLR,femalePairwiseQstsLLS,
                                            femalePairwiseQstsLLW))
colnames(femalePairwiseALLLL) <- c("A","D","E","R","N","S","W")
rownames(femalePairwiseALLLL) <- c("A","D","E","R","N","S","W")



traitNames <- c("DtG","","","","","","",
                "Dt10G","","","","","","",
                "PS","","","","","","",
                "LL","","","","","","")

femalePairwiseQst <- as.data.frame(rbind(femalePairwiseALLDtG,femalePairwiseALLDt10G,femalePairwiseALLPS,femalePairwiseALLLL))
femalePairwiseQstwTrait <- as.data.frame(cbind(traitNames,femalePairwiseQst))


write.csv(femalePairwiseQstwTrait,"tables/femalePairwiseQstwTrait.csv")


##############################################################


########## Between population variances (sum of squares)

### DtG

femaleBetweenPopSSDtGA <- as.data.frame(c(0,SSDtGfemale_AD[1,1],SSDtGfemale_AE[1,1],SSDtGfemale_AN[1,1],SSDtGfemale_AR[1,1],SSDtGfemale_AS[1,1],SSDtGfemale_AW[1,1]))
femaleBetweenPopSSDtGD <- as.data.frame(c(SSDtGfemale_AD[1,1],0,SSDtGfemale_DE[1,1],SSDtGfemale_DN[1,1],SSDtGfemale_DR[1,1],SSDtGfemale_DS[1,1],SSDtGfemale_DW[1,1]))
femaleBetweenPopSSDtGE <- as.data.frame(c(SSDtGfemale_AE[1,1],SSDtGfemale_DE[1,1],0,SSDtGfemale_EN[1,1],SSDtGfemale_ER[1,1],SSDtGfemale_ES[1,1],SSDtGfemale_EW[1,1]))
femaleBetweenPopSSDtGN <- as.data.frame(c(SSDtGfemale_AN[1,1],SSDtGfemale_DN[1,1],SSDtGfemale_EN[1,1],0,SSDtGfemale_NR[1,1],SSDtGfemale_NS[1,1],SSDtGfemale_NW[1,1]))
femaleBetweenPopSSDtGR <- as.data.frame(c(SSDtGfemale_AR[1,1],SSDtGfemale_DR[1,1],SSDtGfemale_ER[1,1],SSDtGfemale_NR[1,1],0,SSDtGfemale_RS[1,1],SSDtGfemale_RW[1,1]))
femaleBetweenPopSSDtGS <- as.data.frame(c(SSDtGfemale_AS[1,1],SSDtGfemale_DS[1,1],SSDtGfemale_ES[1,1],SSDtGfemale_NS[1,1],SSDtGfemale_RS[1,1],0,SSDtGfemale_SW[1,1]))
femaleBetweenPopSSDtGW <- as.data.frame(c(SSDtGfemale_AW[1,1],SSDtGfemale_DW[1,1],SSDtGfemale_EW[1,1],SSDtGfemale_NW[1,1],SSDtGfemale_RW[1,1],SSDtGfemale_SW[1,1],0))

femaleBetweenPopSSALLDtG <- as.data.frame(cbind(femaleBetweenPopSSDtGA,femaleBetweenPopSSDtGD,femaleBetweenPopSSDtGE,
                                                femaleBetweenPopSSDtGN,femaleBetweenPopSSDtGR,femaleBetweenPopSSDtGS,
                                                femaleBetweenPopSSDtGW))
colnames(femaleBetweenPopSSALLDtG) <- c("A","D","E","R","N","S","W")
rownames(femaleBetweenPopSSALLDtG) <- c("A","D","E","R","N","S","W")



### Dt10G

femaleBetweenPopSSDt10GA <- as.data.frame(c(0,SSDt10Gfemale_AD[1,1],SSDt10Gfemale_AE[1,1],SSDt10Gfemale_AN[1,1],SSDt10Gfemale_AR[1,1],SSDt10Gfemale_AS[1,1],SSDt10Gfemale_AW[1,1]))
femaleBetweenPopSSDt10GD <- as.data.frame(c(SSDt10Gfemale_AD[1,1],0,SSDt10Gfemale_DE[1,1],SSDt10Gfemale_DN[1,1],SSDt10Gfemale_DR[1,1],SSDt10Gfemale_DS[1,1],SSDt10Gfemale_DW[1,1]))
femaleBetweenPopSSDt10GE <- as.data.frame(c(SSDt10Gfemale_AE[1,1],SSDt10Gfemale_DE[1,1],0,SSDt10Gfemale_EN[1,1],SSDt10Gfemale_ER[1,1],SSDt10Gfemale_ES[1,1],SSDt10Gfemale_EW[1,1]))
femaleBetweenPopSSDt10GN <- as.data.frame(c(SSDt10Gfemale_AN[1,1],SSDt10Gfemale_DN[1,1],SSDt10Gfemale_EN[1,1],0,SSDt10Gfemale_NR[1,1],SSDt10Gfemale_NS[1,1],SSDt10Gfemale_NW[1,1]))
femaleBetweenPopSSDt10GR <- as.data.frame(c(SSDt10Gfemale_AR[1,1],SSDt10Gfemale_DR[1,1],SSDt10Gfemale_ER[1,1],SSDt10Gfemale_NR[1,1],0,SSDt10Gfemale_RS[1,1],SSDt10Gfemale_RW[1,1]))
femaleBetweenPopSSDt10GS <- as.data.frame(c(SSDt10Gfemale_AS[1,1],SSDt10Gfemale_DS[1,1],SSDt10Gfemale_ES[1,1],SSDt10Gfemale_NS[1,1],SSDt10Gfemale_RS[1,1],0,SSDt10Gfemale_SW[1,1]))
femaleBetweenPopSSDt10GW <- as.data.frame(c(SSDt10Gfemale_AW[1,1],SSDt10Gfemale_DW[1,1],SSDt10Gfemale_EW[1,1],SSDt10Gfemale_NW[1,1],SSDt10Gfemale_RW[1,1],SSDt10Gfemale_SW[1,1],0))

femaleBetweenPopSSALLDt10G <- as.data.frame(cbind(femaleBetweenPopSSDt10GA,femaleBetweenPopSSDt10GD,femaleBetweenPopSSDt10GE,
                                                  femaleBetweenPopSSDt10GN,femaleBetweenPopSSDt10GR,femaleBetweenPopSSDt10GS,
                                                  femaleBetweenPopSSDt10GW))
colnames(femaleBetweenPopSSALLDt10G) <- c("A","D","E","R","N","S","W")
rownames(femaleBetweenPopSSALLDt10G) <- c("A","D","E","R","N","S","W")



### PS

femaleBetweenPopSSPSA <- as.data.frame(c(0,SSPSfemale_AD[1,1],SSPSfemale_AE[1,1],SSPSfemale_AN[1,1],SSPSfemale_AR[1,1],SSPSfemale_AS[1,1],SSPSfemale_AW[1,1]))
femaleBetweenPopSSPSD <- as.data.frame(c(SSPSfemale_AD[1,1],0,SSPSfemale_DE[1,1],SSPSfemale_DN[1,1],SSPSfemale_DR[1,1],SSPSfemale_DS[1,1],SSPSfemale_DW[1,1]))
femaleBetweenPopSSPSE <- as.data.frame(c(SSPSfemale_AE[1,1],SSPSfemale_DE[1,1],0,SSPSfemale_EN[1,1],SSPSfemale_ER[1,1],SSPSfemale_ES[1,1],SSPSfemale_EW[1,1]))
femaleBetweenPopSSPSN <- as.data.frame(c(SSPSfemale_AN[1,1],SSPSfemale_DN[1,1],SSPSfemale_EN[1,1],0,SSPSfemale_NR[1,1],SSPSfemale_NS[1,1],SSPSfemale_NW[1,1]))
femaleBetweenPopSSPSR <- as.data.frame(c(SSPSfemale_AR[1,1],SSPSfemale_DR[1,1],SSPSfemale_ER[1,1],SSPSfemale_NR[1,1],0,SSPSfemale_RS[1,1],SSPSfemale_RW[1,1]))
femaleBetweenPopSSPSS <- as.data.frame(c(SSPSfemale_AS[1,1],SSPSfemale_DS[1,1],SSPSfemale_ES[1,1],SSPSfemale_NS[1,1],SSPSfemale_RS[1,1],0,SSPSfemale_SW[1,1]))
femaleBetweenPopSSPSW <- as.data.frame(c(SSPSfemale_AW[1,1],SSPSfemale_DW[1,1],SSPSfemale_EW[1,1],SSPSfemale_NW[1,1],SSPSfemale_RW[1,1],SSPSfemale_SW[1,1],0))

femaleBetweenPopSSALLPS <- as.data.frame(cbind(femaleBetweenPopSSPSA,femaleBetweenPopSSPSD,femaleBetweenPopSSPSE,
                                               femaleBetweenPopSSPSN,femaleBetweenPopSSPSR,femaleBetweenPopSSPSS,
                                               femaleBetweenPopSSPSW))
colnames(femaleBetweenPopSSALLPS) <- c("A","D","E","R","N","S","W")
rownames(femaleBetweenPopSSALLPS) <- c("A","D","E","R","N","S","W")



### LL

femaleBetweenPopSSLLA <- as.data.frame(c(0,SSLLfemale_AD[1,1],SSLLfemale_AE[1,1],SSLLfemale_AN[1,1],SSLLfemale_AR[1,1],SSLLfemale_AS[1,1],SSLLfemale_AW[1,1]))
femaleBetweenPopSSLLD <- as.data.frame(c(SSLLfemale_AD[1,1],0,SSLLfemale_DE[1,1],SSLLfemale_DN[1,1],SSLLfemale_DR[1,1],SSLLfemale_DS[1,1],SSLLfemale_DW[1,1]))
femaleBetweenPopSSLLE <- as.data.frame(c(SSLLfemale_AE[1,1],SSLLfemale_DE[1,1],0,SSLLfemale_EN[1,1],SSLLfemale_ER[1,1],SSLLfemale_ES[1,1],SSLLfemale_EW[1,1]))
femaleBetweenPopSSLLN <- as.data.frame(c(SSLLfemale_AN[1,1],SSLLfemale_DN[1,1],SSLLfemale_EN[1,1],0,SSLLfemale_NR[1,1],SSLLfemale_NS[1,1],SSLLfemale_NW[1,1]))
femaleBetweenPopSSLLR <- as.data.frame(c(SSLLfemale_AR[1,1],SSLLfemale_DR[1,1],SSLLfemale_ER[1,1],SSLLfemale_NR[1,1],0,SSLLfemale_RS[1,1],SSLLfemale_RW[1,1]))
femaleBetweenPopSSLLS <- as.data.frame(c(SSLLfemale_AS[1,1],SSLLfemale_DS[1,1],SSLLfemale_ES[1,1],SSLLfemale_NS[1,1],SSLLfemale_RS[1,1],0,SSLLfemale_SW[1,1]))
femaleBetweenPopSSLLW <- as.data.frame(c(SSLLfemale_AW[1,1],SSLLfemale_DW[1,1],SSLLfemale_EW[1,1],SSLLfemale_NW[1,1],SSLLfemale_RW[1,1],SSLLfemale_SW[1,1],0))



femaleBetweenPopSSALLLL <- as.data.frame(cbind(femaleBetweenPopSSLLA,femaleBetweenPopSSLLD,femaleBetweenPopSSLLE,
                                               femaleBetweenPopSSLLN,femaleBetweenPopSSLLR,femaleBetweenPopSSLLS,
                                               femaleBetweenPopSSLLW))
colnames(femaleBetweenPopSSALLLL) <- c("A","D","E","R","N","S","W")
rownames(femaleBetweenPopSSALLLL) <- c("A","D","E","R","N","S","W")



femaleBetweenPopSumofSquares <- as.data.frame(rbind(femaleBetweenPopSSALLDtG,femaleBetweenPopSSALLDt10G,femaleBetweenPopSSALLPS,femaleBetweenPopSSALLLL))


write.csv(femaleBetweenPopSumofSquares,"tables/femaleBetweenPopSumofSquares.csv")



##############################################################


########## within population variances (sum of squares)

### DtG

femaleWithinPopSSDtGA <- as.data.frame(c(0,SSDtGfemale_AD[2,1],SSDtGfemale_AE[2,1],SSDtGfemale_AN[2,1],SSDtGfemale_AR[2,1],SSDtGfemale_AS[2,1],SSDtGfemale_AW[2,1]))
femaleWithinPopSSDtGD <- as.data.frame(c(SSDtGfemale_AD[2,1],0,SSDtGfemale_DE[2,1],SSDtGfemale_DN[2,1],SSDtGfemale_DR[2,1],SSDtGfemale_DS[2,1],SSDtGfemale_DW[2,1]))
femaleWithinPopSSDtGE <- as.data.frame(c(SSDtGfemale_AE[2,1],SSDtGfemale_DE[2,1],0,SSDtGfemale_EN[2,1],SSDtGfemale_ER[2,1],SSDtGfemale_ES[2,1],SSDtGfemale_EW[2,1]))
femaleWithinPopSSDtGN <- as.data.frame(c(SSDtGfemale_AN[2,1],SSDtGfemale_DN[2,1],SSDtGfemale_EN[2,1],0,SSDtGfemale_NR[2,1],SSDtGfemale_NS[2,1],SSDtGfemale_NW[2,1]))
femaleWithinPopSSDtGR <- as.data.frame(c(SSDtGfemale_AR[2,1],SSDtGfemale_DR[2,1],SSDtGfemale_ER[2,1],SSDtGfemale_NR[2,1],0,SSDtGfemale_RS[2,1],SSDtGfemale_RW[2,1]))
femaleWithinPopSSDtGS <- as.data.frame(c(SSDtGfemale_AS[2,1],SSDtGfemale_DS[2,1],SSDtGfemale_ES[2,1],SSDtGfemale_NS[2,1],SSDtGfemale_RS[2,1],0,SSDtGfemale_SW[2,1]))
femaleWithinPopSSDtGW <- as.data.frame(c(SSDtGfemale_AW[2,1],SSDtGfemale_DW[2,1],SSDtGfemale_EW[2,1],SSDtGfemale_NW[2,1],SSDtGfemale_RW[2,1],SSDtGfemale_SW[2,1],0))

femaleWithinPopSSALLDtG <- as.data.frame(cbind(femaleWithinPopSSDtGA,femaleWithinPopSSDtGD,femaleWithinPopSSDtGE,
                                               femaleWithinPopSSDtGN,femaleWithinPopSSDtGR,femaleWithinPopSSDtGS,
                                               femaleWithinPopSSDtGW))
colnames(femaleWithinPopSSALLDtG) <- c("A","D","E","R","N","S","W")
rownames(femaleWithinPopSSALLDtG) <- c("A","D","E","R","N","S","W")




### Dt10G

femaleWithinPopSSDt10GA <- as.data.frame(c(0,SSDt10Gfemale_AD[2,1],SSDt10Gfemale_AE[2,1],SSDt10Gfemale_AN[2,1],SSDt10Gfemale_AR[2,1],SSDt10Gfemale_AS[2,1],SSDt10Gfemale_AW[2,1]))
femaleWithinPopSSDt10GD <- as.data.frame(c(SSDt10Gfemale_AD[2,1],0,SSDt10Gfemale_DE[2,1],SSDt10Gfemale_DN[2,1],SSDt10Gfemale_DR[2,1],SSDt10Gfemale_DS[2,1],SSDt10Gfemale_DW[2,1]))
femaleWithinPopSSDt10GE <- as.data.frame(c(SSDt10Gfemale_AE[2,1],SSDt10Gfemale_DE[2,1],0,SSDt10Gfemale_EN[2,1],SSDt10Gfemale_ER[2,1],SSDt10Gfemale_ES[2,1],SSDt10Gfemale_EW[2,1]))
femaleWithinPopSSDt10GN <- as.data.frame(c(SSDt10Gfemale_AN[2,1],SSDt10Gfemale_DN[2,1],SSDt10Gfemale_EN[2,1],0,SSDt10Gfemale_NR[2,1],SSDt10Gfemale_NS[2,1],SSDt10Gfemale_NW[2,1]))
femaleWithinPopSSDt10GR <- as.data.frame(c(SSDt10Gfemale_AR[2,1],SSDt10Gfemale_DR[2,1],SSDt10Gfemale_ER[2,1],SSDt10Gfemale_NR[2,1],0,SSDt10Gfemale_RS[2,1],SSDt10Gfemale_RW[2,1]))
femaleWithinPopSSDt10GS <- as.data.frame(c(SSDt10Gfemale_AS[2,1],SSDt10Gfemale_DS[2,1],SSDt10Gfemale_ES[2,1],SSDt10Gfemale_NS[2,1],SSDt10Gfemale_RS[2,1],0,SSDt10Gfemale_SW[2,1]))
femaleWithinPopSSDt10GW <- as.data.frame(c(SSDt10Gfemale_AW[2,1],SSDt10Gfemale_DW[2,1],SSDt10Gfemale_EW[2,1],SSDt10Gfemale_NW[2,1],SSDt10Gfemale_RW[2,1],SSDt10Gfemale_SW[2,1],0))

femaleWithinPopSSALLDt10G <- as.data.frame(cbind(femaleWithinPopSSDt10GA,femaleWithinPopSSDt10GD,femaleWithinPopSSDt10GE,
                                                 femaleWithinPopSSDt10GN,femaleWithinPopSSDt10GR,femaleWithinPopSSDt10GS,
                                                 femaleWithinPopSSDt10GW))
colnames(femaleWithinPopSSALLDt10G) <- c("A","D","E","R","N","S","W")
rownames(femaleWithinPopSSALLDt10G) <- c("A","D","E","R","N","S","W")


### PS

femaleWithinPopSSPSA <- as.data.frame(c(0,SSPSfemale_AD[2,1],SSPSfemale_AE[2,1],SSPSfemale_AN[2,1],SSPSfemale_AR[2,1],SSPSfemale_AS[2,1],SSPSfemale_AW[2,1]))
femaleWithinPopSSPSD <- as.data.frame(c(SSPSfemale_AD[2,1],0,SSPSfemale_DE[2,1],SSPSfemale_DN[2,1],SSPSfemale_DR[2,1],SSPSfemale_DS[2,1],SSPSfemale_DW[2,1]))
femaleWithinPopSSPSE <- as.data.frame(c(SSPSfemale_AE[2,1],SSPSfemale_DE[2,1],0,SSPSfemale_EN[2,1],SSPSfemale_ER[2,1],SSPSfemale_ES[2,1],SSPSfemale_EW[2,1]))
femaleWithinPopSSPSN <- as.data.frame(c(SSPSfemale_AN[2,1],SSPSfemale_DN[2,1],SSPSfemale_EN[2,1],0,SSPSfemale_NR[2,1],SSPSfemale_NS[2,1],SSPSfemale_NW[2,1]))
femaleWithinPopSSPSR <- as.data.frame(c(SSPSfemale_AR[2,1],SSPSfemale_DR[2,1],SSPSfemale_ER[2,1],SSPSfemale_NR[2,1],0,SSPSfemale_RS[2,1],SSPSfemale_RW[2,1]))
femaleWithinPopSSPSS <- as.data.frame(c(SSPSfemale_AS[2,1],SSPSfemale_DS[2,1],SSPSfemale_ES[2,1],SSPSfemale_NS[2,1],SSPSfemale_RS[2,1],0,SSPSfemale_SW[2,1]))
femaleWithinPopSSPSW <- as.data.frame(c(SSPSfemale_AW[2,1],SSPSfemale_DW[2,1],SSPSfemale_EW[2,1],SSPSfemale_NW[2,1],SSPSfemale_RW[2,1],SSPSfemale_SW[2,1],0))

femaleWithinPopSSALLPS <- as.data.frame(cbind(femaleWithinPopSSPSA,femaleWithinPopSSPSD,femaleWithinPopSSPSE,
                                              femaleWithinPopSSPSN,femaleWithinPopSSPSR,femaleWithinPopSSPSS,
                                              femaleWithinPopSSPSW))
colnames(femaleWithinPopSSALLPS) <- c("A","D","E","R","N","S","W")
rownames(femaleWithinPopSSALLPS) <- c("A","D","E","R","N","S","W")


### LL

femaleWithinPopSSLLA <- as.data.frame(c(0,SSLLfemale_AD[2,1],SSLLfemale_AE[2,1],SSLLfemale_AN[2,1],SSLLfemale_AR[2,1],SSLLfemale_AS[2,1],SSLLfemale_AW[2,1]))
femaleWithinPopSSLLD <- as.data.frame(c(SSLLfemale_AD[2,1],0,SSLLfemale_DE[2,1],SSLLfemale_DN[2,1],SSLLfemale_DR[2,1],SSLLfemale_DS[2,1],SSLLfemale_DW[2,1]))
femaleWithinPopSSLLE <- as.data.frame(c(SSLLfemale_AE[2,1],SSLLfemale_DE[2,1],0,SSLLfemale_EN[2,1],SSLLfemale_ER[2,1],SSLLfemale_ES[2,1],SSLLfemale_EW[2,1]))
femaleWithinPopSSLLN <- as.data.frame(c(SSLLfemale_AN[2,1],SSLLfemale_DN[2,1],SSLLfemale_EN[2,1],0,SSLLfemale_NR[2,1],SSLLfemale_NS[2,1],SSLLfemale_NW[2,1]))
femaleWithinPopSSLLR <- as.data.frame(c(SSLLfemale_AR[2,1],SSLLfemale_DR[2,1],SSLLfemale_ER[2,1],SSLLfemale_NR[2,1],0,SSLLfemale_RS[2,1],SSLLfemale_RW[2,1]))
femaleWithinPopSSLLS <- as.data.frame(c(SSLLfemale_AS[2,1],SSLLfemale_DS[2,1],SSLLfemale_ES[2,1],SSLLfemale_NS[2,1],SSLLfemale_RS[2,1],0,SSLLfemale_SW[2,1]))
femaleWithinPopSSLLW <- as.data.frame(c(SSLLfemale_AW[2,1],SSLLfemale_DW[2,1],SSLLfemale_EW[2,1],SSLLfemale_NW[2,1],SSLLfemale_RW[2,1],SSLLfemale_SW[2,1],0))



femaleWithinPopSSALLLL <- as.data.frame(cbind(femaleWithinPopSSLLA,femaleWithinPopSSLLD,femaleWithinPopSSLLE,
                                              femaleWithinPopSSLLN,femaleWithinPopSSLLR,femaleWithinPopSSLLS,
                                              femaleWithinPopSSLLW))
colnames(femaleWithinPopSSALLLL) <- c("A","D","E","R","N","S","W")
rownames(femaleWithinPopSSALLLL) <- c("A","D","E","R","N","S","W")



traitNames <- c("DtG","","","","","","",
              "Dt10G","","","","","","",
              "PS","","","","","","",
              "LL","","","","","","")


femaleWithinPopSumofSquares <- as.data.frame(rbind(femaleWithinPopSSALLDtG,femaleWithinPopSSALLDt10G,femaleWithinPopSSALLPS,femaleWithinPopSSALLLL))
femaleWithinPopSumofSquares <- as.data.frame(cbind(traitNames,femaleWithinPopSumofSquares))


write.csv(femaleWithinPopSumofSquares,"tables/femaleWithinPopSumofSquares.csv")




########################################################

### pairwise Qst

##males

malemanovaAD <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("A", "D"))


sum_aov_male_AD <- summary.aov(malemanovaAD)

#extracting each trait's sum of squares
SSDtGmale_AD <- as.data.frame(sum_aov_male_AD$` Response DtG`$`Sum Sq`)
SSDt10Gmale_AD <- as.data.frame(sum_aov_male_AD$` Response Dt10G`$`Sum Sq`)
SSPSmale_AD <- as.data.frame(sum_aov_male_AD$` Response PS`$`Sum Sq`)
SSLLmale_AD<- as.data.frame(sum_aov_male_AD$` Response LL`$`Sum Sq`)
SSDtPmale_AD<- as.data.frame(sum_aov_male_AD$` Response DtP`$`Sum Sq`)
SSNPmale_AD<- as.data.frame(sum_aov_male_AD$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_AD <- SSDtGmale_AD[1,1]/(SSDtGmale_AD[2,1]+SSDtGmale_AD[1,1])
Dt10Gmale_AD <- SSDt10Gmale_AD[1,1]/(SSDt10Gmale_AD[2,1]+SSDt10Gmale_AD[1,1])
PSmale_AD <- SSPSmale_AD[1,1]/(SSPSmale_AD[2,1]+SSPSmale_AD[1,1])
LLmale_AD <-  SSLLmale_AD[1,1]/(SSLLmale_AD[2,1]+SSLLmale_AD[1,1])
DtPmale_AD <-  SSDtPmale_AD[1,1]/(SSDtPmale_AD[2,1]+SSDtPmale_AD[1,1])
NPmale_AD <-  SSNPmale_AD[1,1]/(SSNPmale_AD[2,1]+SSNPmale_AD[1,1])

QstmaleAD <- as.data.frame(c(DtGmale_AD, Dt10Gmale_AD, PSmale_AD, LLmale_AD, DtPmale_AD, NPmale_AD ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaAE <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("A", "E"))


sum_aov_male_AE <- summary.aov(malemanovaAE)

#extracting each trait's sum of squares
SSDtGmale_AE <- as.data.frame(sum_aov_male_AE$` Response DtG`$`Sum Sq`)
SSDt10Gmale_AE <- as.data.frame(sum_aov_male_AE$` Response Dt10G`$`Sum Sq`)
SSPSmale_AE <- as.data.frame(sum_aov_male_AE$` Response PS`$`Sum Sq`)
SSLLmale_AE<- as.data.frame(sum_aov_male_AE$` Response LL`$`Sum Sq`)
SSDtPmale_AE<- as.data.frame(sum_aov_male_AE$` Response DtP`$`Sum Sq`)
SSNPmale_AE<- as.data.frame(sum_aov_male_AE$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_AE <- SSDtGmale_AE[1,1]/(SSDtGmale_AE[2,1]+SSDtGmale_AE[1,1])
Dt10Gmale_AE <- SSDt10Gmale_AE[1,1]/(SSDt10Gmale_AE[2,1]+SSDt10Gmale_AE[1,1])
PSmale_AE <- SSPSmale_AE[1,1]/(SSPSmale_AE[2,1]+SSPSmale_AE[1,1])
LLmale_AE <-  SSLLmale_AE[1,1]/(SSLLmale_AE[2,1]+SSLLmale_AE[1,1])
DtPmale_AE <-  SSDtPmale_AE[1,1]/(SSDtPmale_AE[2,1]+SSDtPmale_AE[1,1])
NPmale_AE <-  SSNPmale_AE[1,1]/(SSNPmale_AE[2,1]+SSNPmale_AE[1,1])

QstmaleAE <- as.data.frame(c(DtGmale_AE, Dt10Gmale_AE, PSmale_AE, LLmale_AE, DtPmale_AE, NPmale_AE ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaAN <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("A", "N"))


sum_aov_male_AN <- summary.aov(malemanovaAN)

#extracting each trait's sum of squares
SSDtGmale_AN <- as.data.frame(sum_aov_male_AN$` Response DtG`$`Sum Sq`)
SSDt10Gmale_AN <- as.data.frame(sum_aov_male_AN$` Response Dt10G`$`Sum Sq`)
SSPSmale_AN <- as.data.frame(sum_aov_male_AN$` Response PS`$`Sum Sq`)
SSLLmale_AN<- as.data.frame(sum_aov_male_AN$` Response LL`$`Sum Sq`)
SSDtPmale_AN<- as.data.frame(sum_aov_male_AN$` Response DtP`$`Sum Sq`)
SSNPmale_AN<- as.data.frame(sum_aov_male_AN$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_AN <- SSDtGmale_AN[1,1]/(SSDtGmale_AN[2,1]+SSDtGmale_AN[1,1])
Dt10Gmale_AN <- SSDt10Gmale_AN[1,1]/(SSDt10Gmale_AN[2,1]+SSDt10Gmale_AN[1,1])
PSmale_AN <- SSPSmale_AN[1,1]/(SSPSmale_AN[2,1]+SSPSmale_AN[1,1])
LLmale_AN <-  SSLLmale_AN[1,1]/(SSLLmale_AN[2,1]+SSLLmale_AN[1,1])
DtPmale_AN <-  SSDtPmale_AN[1,1]/(SSDtPmale_AN[2,1]+SSDtPmale_AN[1,1])
NPmale_AN <-  SSNPmale_AN[1,1]/(SSNPmale_AN[2,1]+SSNPmale_AN[1,1])

QstmaleAN <- as.data.frame(c(DtGmale_AN, Dt10Gmale_AN, PSmale_AN, LLmale_AN, DtPmale_AN, NPmale_AN ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaAR <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("A", "R"))


sum_aov_male_AR <- summary.aov(malemanovaAR)

#extracting each trait's sum of squares
SSDtGmale_AR <- as.data.frame(sum_aov_male_AR$` Response DtG`$`Sum Sq`)
SSDt10Gmale_AR <- as.data.frame(sum_aov_male_AR$` Response Dt10G`$`Sum Sq`)
SSPSmale_AR <- as.data.frame(sum_aov_male_AR$` Response PS`$`Sum Sq`)
SSLLmale_AR<- as.data.frame(sum_aov_male_AR$` Response LL`$`Sum Sq`)
SSDtPmale_AR<- as.data.frame(sum_aov_male_AR$` Response DtP`$`Sum Sq`)
SSNPmale_AR<- as.data.frame(sum_aov_male_AR$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_AR <- SSDtGmale_AR[1,1]/(SSDtGmale_AR[2,1]+SSDtGmale_AR[1,1])
Dt10Gmale_AR <- SSDt10Gmale_AR[1,1]/(SSDt10Gmale_AR[2,1]+SSDt10Gmale_AR[1,1])
PSmale_AR <- SSPSmale_AR[1,1]/(SSPSmale_AR[2,1]+SSPSmale_AR[1,1])
LLmale_AR <-  SSLLmale_AR[1,1]/(SSLLmale_AR[2,1]+SSLLmale_AR[1,1])
DtPmale_AR <-  SSDtPmale_AR[1,1]/(SSDtPmale_AR[2,1]+SSDtPmale_AR[1,1])
NPmale_AR <-  SSNPmale_AR[1,1]/(SSNPmale_AR[2,1]+SSNPmale_AR[1,1])

QstmaleAR <- as.data.frame(c(DtGmale_AR, Dt10Gmale_AR, PSmale_AR, LLmale_AR, DtPmale_AR, NPmale_AR ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaAS <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("A", "S"))


sum_aov_male_AS <- summary.aov(malemanovaAS)

#extracting each trait's sum of squares
SSDtGmale_AS <- as.data.frame(sum_aov_male_AS$` Response DtG`$`Sum Sq`)
SSDt10Gmale_AS <- as.data.frame(sum_aov_male_AS$` Response Dt10G`$`Sum Sq`)
SSPSmale_AS <- as.data.frame(sum_aov_male_AS$` Response PS`$`Sum Sq`)
SSLLmale_AS<- as.data.frame(sum_aov_male_AS$` Response LL`$`Sum Sq`)
SSDtPmale_AS<- as.data.frame(sum_aov_male_AS$` Response DtP`$`Sum Sq`)
SSNPmale_AS<- as.data.frame(sum_aov_male_AS$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_AS <- SSDtGmale_AS[1,1]/(SSDtGmale_AS[2,1]+SSDtGmale_AS[1,1])
Dt10Gmale_AS <- SSDt10Gmale_AS[1,1]/(SSDt10Gmale_AS[2,1]+SSDt10Gmale_AS[1,1])
PSmale_AS <- SSPSmale_AS[1,1]/(SSPSmale_AS[2,1]+SSPSmale_AS[1,1])
LLmale_AS <-  SSLLmale_AS[1,1]/(SSLLmale_AS[2,1]+SSLLmale_AS[1,1])
DtPmale_AS <-  SSDtPmale_AS[1,1]/(SSDtPmale_AS[2,1]+SSDtPmale_AS[1,1])
NPmale_AS <-  SSNPmale_AS[1,1]/(SSNPmale_AS[2,1]+SSNPmale_AS[1,1])

QstmaleAS <- as.data.frame(c(DtGmale_AS, Dt10Gmale_AS, PSmale_AS, LLmale_AS, DtPmale_AS, NPmale_AS ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaAW <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("A", "W"))


sum_aov_male_AW <- summary.aov(malemanovaAW)

#extracting each trait's sum of squares
SSDtGmale_AW <- as.data.frame(sum_aov_male_AW$` Response DtG`$`Sum Sq`)
SSDt10Gmale_AW <- as.data.frame(sum_aov_male_AW$` Response Dt10G`$`Sum Sq`)
SSPSmale_AW <- as.data.frame(sum_aov_male_AW$` Response PS`$`Sum Sq`)
SSLLmale_AW<- as.data.frame(sum_aov_male_AW$` Response LL`$`Sum Sq`)
SSDtPmale_AW<- as.data.frame(sum_aov_male_AW$` Response DtP`$`Sum Sq`)
SSNPmale_AW<- as.data.frame(sum_aov_male_AW$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_AW <- SSDtGmale_AW[1,1]/(SSDtGmale_AW[2,1]+SSDtGmale_AW[1,1])
Dt10Gmale_AW <- SSDt10Gmale_AW[1,1]/(SSDt10Gmale_AW[2,1]+SSDt10Gmale_AW[1,1])
PSmale_AW <- SSPSmale_AW[1,1]/(SSPSmale_AW[2,1]+SSPSmale_AW[1,1])
LLmale_AW <-  SSLLmale_AW[1,1]/(SSLLmale_AW[2,1]+SSLLmale_AW[1,1])
DtPmale_AW <-  SSDtPmale_AW[1,1]/(SSDtPmale_AW[2,1]+SSDtPmale_AW[1,1])
NPmale_AW <-  SSNPmale_AW[1,1]/(SSNPmale_AW[2,1]+SSNPmale_AW[1,1])

QstmaleAW <- as.data.frame(c(DtGmale_AW, Dt10Gmale_AW, PSmale_AW, LLmale_AW, DtPmale_AW, NPmale_AW ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaDE <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("D", "E"))


sum_aov_male_DE <- summary.aov(malemanovaDE)

#extracting each trait's sum of squares
SSDtGmale_DE <- as.data.frame(sum_aov_male_DE$` Response DtG`$`Sum Sq`)
SSDt10Gmale_DE <- as.data.frame(sum_aov_male_DE$` Response Dt10G`$`Sum Sq`)
SSPSmale_DE <- as.data.frame(sum_aov_male_DE$` Response PS`$`Sum Sq`)
SSLLmale_DE<- as.data.frame(sum_aov_male_DE$` Response LL`$`Sum Sq`)
SSDtPmale_DE<- as.data.frame(sum_aov_male_DE$` Response DtP`$`Sum Sq`)
SSNPmale_DE<- as.data.frame(sum_aov_male_DE$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_DE <- SSDtGmale_DE[1,1]/(SSDtGmale_DE[2,1]+SSDtGmale_DE[1,1])
Dt10Gmale_DE <- SSDt10Gmale_DE[1,1]/(SSDt10Gmale_DE[2,1]+SSDt10Gmale_DE[1,1])
PSmale_DE <- SSPSmale_DE[1,1]/(SSPSmale_DE[2,1]+SSPSmale_DE[1,1])
LLmale_DE <-  SSLLmale_DE[1,1]/(SSLLmale_DE[2,1]+SSLLmale_DE[1,1])
DtPmale_DE <-  SSDtPmale_DE[1,1]/(SSDtPmale_DE[2,1]+SSDtPmale_DE[1,1])
NPmale_DE <-  SSNPmale_DE[1,1]/(SSNPmale_DE[2,1]+SSNPmale_DE[1,1])

QstmaleDE <- as.data.frame(c(DtGmale_DE, Dt10Gmale_DE, PSmale_DE, LLmale_DE, DtPmale_DE, NPmale_DE ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaDN <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("D", "N"))


sum_aov_male_DN <- summary.aov(malemanovaDN)

#extracting each trait's sum of squares
SSDtGmale_DN <- as.data.frame(sum_aov_male_DN$` Response DtG`$`Sum Sq`)
SSDt10Gmale_DN <- as.data.frame(sum_aov_male_DN$` Response Dt10G`$`Sum Sq`)
SSPSmale_DN <- as.data.frame(sum_aov_male_DN$` Response PS`$`Sum Sq`)
SSLLmale_DN<- as.data.frame(sum_aov_male_DN$` Response LL`$`Sum Sq`)
SSDtPmale_DN<- as.data.frame(sum_aov_male_DN$` Response DtP`$`Sum Sq`)
SSNPmale_DN<- as.data.frame(sum_aov_male_DN$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_DN <- SSDtGmale_DN[1,1]/(SSDtGmale_DN[2,1]+SSDtGmale_DN[1,1])
Dt10Gmale_DN <- SSDt10Gmale_DN[1,1]/(SSDt10Gmale_DN[2,1]+SSDt10Gmale_DN[1,1])
PSmale_DN <- SSPSmale_DN[1,1]/(SSPSmale_DN[2,1]+SSPSmale_DN[1,1])
LLmale_DN <-  SSLLmale_DN[1,1]/(SSLLmale_DN[2,1]+SSLLmale_DN[1,1])
DtPmale_DN <-  SSDtPmale_DN[1,1]/(SSDtPmale_DN[2,1]+SSDtPmale_DN[1,1])
NPmale_DN <-  SSNPmale_DN[1,1]/(SSNPmale_DN[2,1]+SSNPmale_DN[1,1])

QstmaleDN <- as.data.frame(c(DtGmale_DN, Dt10Gmale_DN, PSmale_DN, LLmale_DN, DtPmale_DN, NPmale_DN ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaDR <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("D", "R"))


sum_aov_male_DR <- summary.aov(malemanovaDR)

#extracting each trait's sum of squares
SSDtGmale_DR <- as.data.frame(sum_aov_male_DR$` Response DtG`$`Sum Sq`)
SSDt10Gmale_DR <- as.data.frame(sum_aov_male_DR$` Response Dt10G`$`Sum Sq`)
SSPSmale_DR <- as.data.frame(sum_aov_male_DR$` Response PS`$`Sum Sq`)
SSLLmale_DR<- as.data.frame(sum_aov_male_DR$` Response LL`$`Sum Sq`)
SSDtPmale_DR<- as.data.frame(sum_aov_male_DR$` Response DtP`$`Sum Sq`)
SSNPmale_DR<- as.data.frame(sum_aov_male_DR$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_DR <- SSDtGmale_DR[1,1]/(SSDtGmale_DR[2,1]+SSDtGmale_DR[1,1])
Dt10Gmale_DR <- SSDt10Gmale_DR[1,1]/(SSDt10Gmale_DR[2,1]+SSDt10Gmale_DR[1,1])
PSmale_DR <- SSPSmale_DR[1,1]/(SSPSmale_DR[2,1]+SSPSmale_DR[1,1])
LLmale_DR <-  SSLLmale_DR[1,1]/(SSLLmale_DR[2,1]+SSLLmale_DR[1,1])
DtPmale_DR <-  SSDtPmale_DR[1,1]/(SSDtPmale_DR[2,1]+SSDtPmale_DR[1,1])
NPmale_DR <-  SSNPmale_DR[1,1]/(SSNPmale_DR[2,1]+SSNPmale_DR[1,1])

QstmaleDR <- as.data.frame(c(DtGmale_DR, Dt10Gmale_DR, PSmale_DR, LLmale_DR, DtPmale_DR, NPmale_DR ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaDS <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("D", "S"))


sum_aov_male_DS <- summary.aov(malemanovaDS)

#extracting each trait's sum of squares
SSDtGmale_DS <- as.data.frame(sum_aov_male_DS$` Response DtG`$`Sum Sq`)
SSDt10Gmale_DS <- as.data.frame(sum_aov_male_DS$` Response Dt10G`$`Sum Sq`)
SSPSmale_DS <- as.data.frame(sum_aov_male_DS$` Response PS`$`Sum Sq`)
SSLLmale_DS<- as.data.frame(sum_aov_male_DS$` Response LL`$`Sum Sq`)
SSDtPmale_DS<- as.data.frame(sum_aov_male_DS$` Response DtP`$`Sum Sq`)
SSNPmale_DS<- as.data.frame(sum_aov_male_DS$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_DS <- SSDtGmale_DS[1,1]/(SSDtGmale_DS[2,1]+SSDtGmale_DS[1,1])
Dt10Gmale_DS <- SSDt10Gmale_DS[1,1]/(SSDt10Gmale_DS[2,1]+SSDt10Gmale_DS[1,1])
PSmale_DS <- SSPSmale_DS[1,1]/(SSPSmale_DS[2,1]+SSPSmale_DS[1,1])
LLmale_DS <-  SSLLmale_DS[1,1]/(SSLLmale_DS[2,1]+SSLLmale_DS[1,1])
DtPmale_DS <-  SSDtPmale_DS[1,1]/(SSDtPmale_DS[2,1]+SSDtPmale_DS[1,1])
NPmale_DS <-  SSNPmale_DS[1,1]/(SSNPmale_DS[2,1]+SSNPmale_DS[1,1])

QstmaleDS <- as.data.frame(c(DtGmale_DS, Dt10Gmale_DS, PSmale_DS, LLmale_DS, DtPmale_DS, NPmale_DS ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaDW <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("D", "W"))


sum_aov_male_DW <- summary.aov(malemanovaDW)

#extracting each trait's sum of squares
SSDtGmale_DW <- as.data.frame(sum_aov_male_DW$` Response DtG`$`Sum Sq`)
SSDt10Gmale_DW <- as.data.frame(sum_aov_male_DW$` Response Dt10G`$`Sum Sq`)
SSPSmale_DW <- as.data.frame(sum_aov_male_DW$` Response PS`$`Sum Sq`)
SSLLmale_DW<- as.data.frame(sum_aov_male_DW$` Response LL`$`Sum Sq`)
SSDtPmale_DW<- as.data.frame(sum_aov_male_DW$` Response DtP`$`Sum Sq`)
SSNPmale_DW<- as.data.frame(sum_aov_male_DW$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_DW <- SSDtGmale_DW[1,1]/(SSDtGmale_DW[2,1]+SSDtGmale_DW[1,1])
Dt10Gmale_DW <- SSDt10Gmale_DW[1,1]/(SSDt10Gmale_DW[2,1]+SSDt10Gmale_DW[1,1])
PSmale_DW <- SSPSmale_DW[1,1]/(SSPSmale_DW[2,1]+SSPSmale_DW[1,1])
LLmale_DW <-  SSLLmale_DW[1,1]/(SSLLmale_DW[2,1]+SSLLmale_DW[1,1])
DtPmale_DW <-  SSDtPmale_DW[1,1]/(SSDtPmale_DW[2,1]+SSDtPmale_DW[1,1])
NPmale_DW <-  SSNPmale_DW[1,1]/(SSNPmale_DW[2,1]+SSNPmale_DW[1,1])

QstmaleDW <- as.data.frame(c(DtGmale_DW, Dt10Gmale_DW, PSmale_DW, LLmale_DW, DtPmale_DW, NPmale_DW ), c("DtG","Dt10G","PS","LL","DtP","NP"))



#####
malemanovaEN <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("E", "N"))


sum_aov_male_EN <- summary.aov(malemanovaEN)

#extracting each trait's sum of squares
SSDtGmale_EN <- as.data.frame(sum_aov_male_EN$` Response DtG`$`Sum Sq`)
SSDt10Gmale_EN <- as.data.frame(sum_aov_male_EN$` Response Dt10G`$`Sum Sq`)
SSPSmale_EN <- as.data.frame(sum_aov_male_EN$` Response PS`$`Sum Sq`)
SSLLmale_EN<- as.data.frame(sum_aov_male_EN$` Response LL`$`Sum Sq`)
SSDtPmale_EN<- as.data.frame(sum_aov_male_EN$` Response DtP`$`Sum Sq`)
SSNPmale_EN<- as.data.frame(sum_aov_male_EN$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_EN <- SSDtGmale_EN[1,1]/(SSDtGmale_EN[2,1]+SSDtGmale_EN[1,1])
Dt10Gmale_EN <- SSDt10Gmale_EN[1,1]/(SSDt10Gmale_EN[2,1]+SSDt10Gmale_EN[1,1])
PSmale_EN <- SSPSmale_EN[1,1]/(SSPSmale_EN[2,1]+SSPSmale_EN[1,1])
LLmale_EN <-  SSLLmale_EN[1,1]/(SSLLmale_EN[2,1]+SSLLmale_EN[1,1])
DtPmale_EN <-  SSDtPmale_EN[1,1]/(SSDtPmale_EN[2,1]+SSDtPmale_EN[1,1])
NPmale_EN <-  SSNPmale_EN[1,1]/(SSNPmale_EN[2,1]+SSNPmale_EN[1,1])

QstmaleEN <- as.data.frame(c(DtGmale_EN, Dt10Gmale_EN, PSmale_EN, LLmale_EN, DtPmale_EN, NPmale_EN ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaER <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("E", "R"))


sum_aov_male_ER <- summary.aov(malemanovaER)

#extracting each trait's sum of squares
SSDtGmale_ER <- as.data.frame(sum_aov_male_ER$` Response DtG`$`Sum Sq`)
SSDt10Gmale_ER <- as.data.frame(sum_aov_male_ER$` Response Dt10G`$`Sum Sq`)
SSPSmale_ER <- as.data.frame(sum_aov_male_ER$` Response PS`$`Sum Sq`)
SSLLmale_ER<- as.data.frame(sum_aov_male_ER$` Response LL`$`Sum Sq`)
SSDtPmale_ER<- as.data.frame(sum_aov_male_ER$` Response DtP`$`Sum Sq`)
SSNPmale_ER<- as.data.frame(sum_aov_male_ER$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_ER <- SSDtGmale_ER[1,1]/(SSDtGmale_ER[2,1]+SSDtGmale_ER[1,1])
Dt10Gmale_ER <- SSDt10Gmale_ER[1,1]/(SSDt10Gmale_ER[2,1]+SSDt10Gmale_ER[1,1])
PSmale_ER <- SSPSmale_ER[1,1]/(SSPSmale_ER[2,1]+SSPSmale_ER[1,1])
LLmale_ER <-  SSLLmale_ER[1,1]/(SSLLmale_ER[2,1]+SSLLmale_ER[1,1])
DtPmale_ER <-  SSDtPmale_ER[1,1]/(SSDtPmale_ER[2,1]+SSDtPmale_ER[1,1])
NPmale_ER <-  SSNPmale_ER[1,1]/(SSNPmale_ER[2,1]+SSNPmale_ER[1,1])

QstmaleER <- as.data.frame(c(DtGmale_ER, Dt10Gmale_ER, PSmale_ER, LLmale_ER, DtPmale_ER, NPmale_ER ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaES <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("E", "S"))


sum_aov_male_ES <- summary.aov(malemanovaES)

#extracting each trait's sum of squares
SSDtGmale_ES <- as.data.frame(sum_aov_male_ES$` Response DtG`$`Sum Sq`)
SSDt10Gmale_ES <- as.data.frame(sum_aov_male_ES$` Response Dt10G`$`Sum Sq`)
SSPSmale_ES <- as.data.frame(sum_aov_male_ES$` Response PS`$`Sum Sq`)
SSLLmale_ES<- as.data.frame(sum_aov_male_ES$` Response LL`$`Sum Sq`)
SSDtPmale_ES<- as.data.frame(sum_aov_male_ES$` Response DtP`$`Sum Sq`)
SSNPmale_ES<- as.data.frame(sum_aov_male_ES$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_ES <- SSDtGmale_ES[1,1]/(SSDtGmale_ES[2,1]+SSDtGmale_ES[1,1])
Dt10Gmale_ES <- SSDt10Gmale_ES[1,1]/(SSDt10Gmale_ES[2,1]+SSDt10Gmale_ES[1,1])
PSmale_ES <- SSPSmale_ES[1,1]/(SSPSmale_ES[2,1]+SSPSmale_ES[1,1])
LLmale_ES <-  SSLLmale_ES[1,1]/(SSLLmale_ES[2,1]+SSLLmale_ES[1,1])
DtPmale_ES <-  SSDtPmale_ES[1,1]/(SSDtPmale_ES[2,1]+SSDtPmale_ES[1,1])
NPmale_ES <-  SSNPmale_ES[1,1]/(SSNPmale_ES[2,1]+SSNPmale_ES[1,1])

QstmaleES <- as.data.frame(c(DtGmale_ES, Dt10Gmale_ES, PSmale_ES, LLmale_ES, DtPmale_ES, NPmale_ES ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaEW <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("E", "W"))


sum_aov_male_EW <- summary.aov(malemanovaEW)

#extracting each trait's sum of squares
SSDtGmale_EW <- as.data.frame(sum_aov_male_EW$` Response DtG`$`Sum Sq`)
SSDt10Gmale_EW <- as.data.frame(sum_aov_male_EW$` Response Dt10G`$`Sum Sq`)
SSPSmale_EW <- as.data.frame(sum_aov_male_EW$` Response PS`$`Sum Sq`)
SSLLmale_EW<- as.data.frame(sum_aov_male_EW$` Response LL`$`Sum Sq`)
SSDtPmale_EW<- as.data.frame(sum_aov_male_EW$` Response DtP`$`Sum Sq`)
SSNPmale_EW<- as.data.frame(sum_aov_male_EW$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_EW <- SSDtGmale_EW[1,1]/(SSDtGmale_EW[2,1]+SSDtGmale_EW[1,1])
Dt10Gmale_EW <- SSDt10Gmale_EW[1,1]/(SSDt10Gmale_EW[2,1]+SSDt10Gmale_EW[1,1])
PSmale_EW <- SSPSmale_EW[1,1]/(SSPSmale_EW[2,1]+SSPSmale_EW[1,1])
LLmale_EW <-  SSLLmale_EW[1,1]/(SSLLmale_EW[2,1]+SSLLmale_EW[1,1])
DtPmale_EW <-  SSDtPmale_EW[1,1]/(SSDtPmale_EW[2,1]+SSDtPmale_EW[1,1])
NPmale_EW <-  SSNPmale_EW[1,1]/(SSNPmale_EW[2,1]+SSNPmale_EW[1,1])

QstmaleEW <- as.data.frame(c(DtGmale_EW, Dt10Gmale_EW, PSmale_EW, LLmale_EW, DtPmale_EW, NPmale_EW ), c("DtG","Dt10G","PS","LL","DtP","NP"))



#####
malemanovaNR <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("N", "R"))


sum_aov_male_NR <- summary.aov(malemanovaNR)

#extracting each trait's sum of squares
SSDtGmale_NR <- as.data.frame(sum_aov_male_NR$` Response DtG`$`Sum Sq`)
SSDt10Gmale_NR <- as.data.frame(sum_aov_male_NR$` Response Dt10G`$`Sum Sq`)
SSPSmale_NR <- as.data.frame(sum_aov_male_NR$` Response PS`$`Sum Sq`)
SSLLmale_NR<- as.data.frame(sum_aov_male_NR$` Response LL`$`Sum Sq`)
SSDtPmale_NR<- as.data.frame(sum_aov_male_NR$` Response DtP`$`Sum Sq`)
SSNPmale_NR<- as.data.frame(sum_aov_male_NR$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_NR <- SSDtGmale_NR[1,1]/(SSDtGmale_NR[2,1]+SSDtGmale_NR[1,1])
Dt10Gmale_NR <- SSDt10Gmale_NR[1,1]/(SSDt10Gmale_NR[2,1]+SSDt10Gmale_NR[1,1])
PSmale_NR <- SSPSmale_NR[1,1]/(SSPSmale_NR[2,1]+SSPSmale_NR[1,1])
LLmale_NR <-  SSLLmale_NR[1,1]/(SSLLmale_NR[2,1]+SSLLmale_NR[1,1])
DtPmale_NR <-  SSDtPmale_NR[1,1]/(SSDtPmale_NR[2,1]+SSDtPmale_NR[1,1])
NPmale_NR <-  SSNPmale_NR[1,1]/(SSNPmale_NR[2,1]+SSNPmale_NR[1,1])

QstmaleNR <- as.data.frame(c(DtGmale_NR, Dt10Gmale_NR, PSmale_NR, LLmale_NR, DtPmale_NR, NPmale_NR ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaNS <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("N", "S"))


sum_aov_male_NS <- summary.aov(malemanovaNS)

#extracting each trait's sum of squares
SSDtGmale_NS <- as.data.frame(sum_aov_male_NS$` Response DtG`$`Sum Sq`)
SSDt10Gmale_NS <- as.data.frame(sum_aov_male_NS$` Response Dt10G`$`Sum Sq`)
SSPSmale_NS <- as.data.frame(sum_aov_male_NS$` Response PS`$`Sum Sq`)
SSLLmale_NS<- as.data.frame(sum_aov_male_NS$` Response LL`$`Sum Sq`)
SSDtPmale_NS<- as.data.frame(sum_aov_male_NS$` Response DtP`$`Sum Sq`)
SSNPmale_NS<- as.data.frame(sum_aov_male_NS$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_NS <- SSDtGmale_NS[1,1]/(SSDtGmale_NS[2,1]+SSDtGmale_NS[1,1])
Dt10Gmale_NS <- SSDt10Gmale_NS[1,1]/(SSDt10Gmale_NS[2,1]+SSDt10Gmale_NS[1,1])
PSmale_NS <- SSPSmale_NS[1,1]/(SSPSmale_NS[2,1]+SSPSmale_NS[1,1])
LLmale_NS <-  SSLLmale_NS[1,1]/(SSLLmale_NS[2,1]+SSLLmale_NS[1,1])
DtPmale_NS <-  SSDtPmale_NS[1,1]/(SSDtPmale_NS[2,1]+SSDtPmale_NS[1,1])
NPmale_NS <-  SSNPmale_NS[1,1]/(SSNPmale_NS[2,1]+SSNPmale_NS[1,1])

QstmaleNS <- as.data.frame(c(DtGmale_NS, Dt10Gmale_NS, PSmale_NS, LLmale_NS, DtPmale_NS, NPmale_NS ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaNW <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("N", "W"))


sum_aov_male_NW <- summary.aov(malemanovaNW)

#extracting each trait's sum of squares
SSDtGmale_NW <- as.data.frame(sum_aov_male_NW$` Response DtG`$`Sum Sq`)
SSDt10Gmale_NW <- as.data.frame(sum_aov_male_NW$` Response Dt10G`$`Sum Sq`)
SSPSmale_NW <- as.data.frame(sum_aov_male_NW$` Response PS`$`Sum Sq`)
SSLLmale_NW<- as.data.frame(sum_aov_male_NW$` Response LL`$`Sum Sq`)
SSDtPmale_NW<- as.data.frame(sum_aov_male_NW$` Response DtP`$`Sum Sq`)
SSNPmale_NW<- as.data.frame(sum_aov_male_NW$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_NW <- SSDtGmale_NW[1,1]/(SSDtGmale_NW[2,1]+SSDtGmale_NW[1,1])
Dt10Gmale_NW <- SSDt10Gmale_NW[1,1]/(SSDt10Gmale_NW[2,1]+SSDt10Gmale_NW[1,1])
PSmale_NW <- SSPSmale_NW[1,1]/(SSPSmale_NW[2,1]+SSPSmale_NW[1,1])
LLmale_NW <-  SSLLmale_NW[1,1]/(SSLLmale_NW[2,1]+SSLLmale_NW[1,1])
DtPmale_NW <-  SSDtPmale_NW[1,1]/(SSDtPmale_NW[2,1]+SSDtPmale_NW[1,1])
NPmale_NW <-  SSNPmale_NW[1,1]/(SSNPmale_NW[2,1]+SSNPmale_NW[1,1])

QstmaleNW <- as.data.frame(c(DtGmale_NW, Dt10Gmale_NW, PSmale_NW, LLmale_NW, DtPmale_NW, NPmale_NW ), c("DtG","Dt10G","PS","LL","DtP","NP"))


#####
malemanovaRS <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("R", "S"))


sum_aov_male_RS <- summary.aov(malemanovaNR)

#extracting each trait's sum of squares
SSDtGmale_RS <- as.data.frame(sum_aov_male_RS$` Response DtG`$`Sum Sq`)
SSDt10Gmale_RS <- as.data.frame(sum_aov_male_RS$` Response Dt10G`$`Sum Sq`)
SSPSmale_RS <- as.data.frame(sum_aov_male_RS$` Response PS`$`Sum Sq`)
SSLLmale_RS<- as.data.frame(sum_aov_male_RS$` Response LL`$`Sum Sq`)
SSDtPmale_RS<- as.data.frame(sum_aov_male_RS$` Response DtP`$`Sum Sq`)
SSNPmale_RS<- as.data.frame(sum_aov_male_RS$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_RS <- SSDtGmale_RS[1,1]/(SSDtGmale_RS[2,1]+SSDtGmale_RS[1,1])
Dt10Gmale_RS <- SSDt10Gmale_RS[1,1]/(SSDt10Gmale_RS[2,1]+SSDt10Gmale_RS[1,1])
PSmale_RS <- SSPSmale_RS[1,1]/(SSPSmale_RS[2,1]+SSPSmale_RS[1,1])
LLmale_RS <-  SSLLmale_RS[1,1]/(SSLLmale_RS[2,1]+SSLLmale_RS[1,1])
DtPmale_RS <-  SSDtPmale_RS[1,1]/(SSDtPmale_RS[2,1]+SSDtPmale_RS[1,1])
NPmale_RS <-  SSNPmale_RS[1,1]/(SSNPmale_RS[2,1]+SSNPmale_RS[1,1])

QstmaleRS <- as.data.frame(c(DtGmale_RS, Dt10Gmale_RS, PSmale_RS, LLmale_RS, DtPmale_RS, NPmale_RS ), c("DtG","Dt10G","PS","LL","DtP","NP"))



#####
malemanovaRW <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("R", "W"))


sum_aov_male_RW <- summary.aov(malemanovaRW)

#extracting each trait's sum of squares
SSDtGmale_RW <- as.data.frame(sum_aov_male_RW$` Response DtG`$`Sum Sq`)
SSDt10Gmale_RW <- as.data.frame(sum_aov_male_RW$` Response Dt10G`$`Sum Sq`)
SSPSmale_RW <- as.data.frame(sum_aov_male_RW$` Response PS`$`Sum Sq`)
SSLLmale_RW<- as.data.frame(sum_aov_male_RW$` Response LL`$`Sum Sq`)
SSDtPmale_RW<- as.data.frame(sum_aov_male_RW$` Response DtP`$`Sum Sq`)
SSNPmale_RW<- as.data.frame(sum_aov_male_RW$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_RW <- SSDtGmale_RW[1,1]/(SSDtGmale_RW[2,1]+SSDtGmale_RW[1,1])
Dt10Gmale_RW <- SSDt10Gmale_RW[1,1]/(SSDt10Gmale_RW[2,1]+SSDt10Gmale_RW[1,1])
PSmale_RW <- SSPSmale_RW[1,1]/(SSPSmale_RW[2,1]+SSPSmale_RW[1,1])
LLmale_RW <-  SSLLmale_RW[1,1]/(SSLLmale_RW[2,1]+SSLLmale_RW[1,1])
DtPmale_RW <-  SSDtPmale_RW[1,1]/(SSDtPmale_RW[2,1]+SSDtPmale_RW[1,1])
NPmale_RW <-  SSNPmale_RW[1,1]/(SSNPmale_RW[2,1]+SSNPmale_RW[1,1])

QstmaleRW <- as.data.frame(c(DtGmale_RW, Dt10Gmale_RW, PSmale_RW, LLmale_RW, DtPmale_RW, NPmale_RW ), c("DtG","Dt10G","PS","LL","DtP","NP"))

#####
malemanovaSW <- manova(cbind(DtG, Dt10G, PS, LL, DtP, NP) ~
                           Population + Population/Individual, 
                         data = maleTransALL,
                         subset = Population %in% c("S", "W"))


sum_aov_male_SW <- summary.aov(malemanovaSW)

#extracting each trait's sum of squares
SSDtGmale_SW <- as.data.frame(sum_aov_male_SW$` Response DtG`$`Sum Sq`)
SSDt10Gmale_SW <- as.data.frame(sum_aov_male_SW$` Response Dt10G`$`Sum Sq`)
SSPSmale_SW <- as.data.frame(sum_aov_male_SW$` Response PS`$`Sum Sq`)
SSLLmale_SW<- as.data.frame(sum_aov_male_SW$` Response LL`$`Sum Sq`)
SSDtPmale_SW<- as.data.frame(sum_aov_male_SW$` Response DtP`$`Sum Sq`)
SSNPmale_SW<- as.data.frame(sum_aov_male_SW$` Response NP`$`Sum Sq`)

#haploid Qst calculation
DtGmale_SW <- SSDtGmale_SW[1,1]/(SSDtGmale_SW[2,1]+SSDtGmale_SW[1,1])
Dt10Gmale_SW <- SSDt10Gmale_SW[1,1]/(SSDt10Gmale_SW[2,1]+SSDt10Gmale_SW[1,1])
PSmale_SW <- SSPSmale_SW[1,1]/(SSPSmale_SW[2,1]+SSPSmale_SW[1,1])
LLmale_SW <-  SSLLmale_SW[1,1]/(SSLLmale_SW[2,1]+SSLLmale_SW[1,1])
DtPmale_SW <-  SSDtPmale_SW[1,1]/(SSDtPmale_SW[2,1]+SSDtPmale_SW[1,1])
NPmale_SW <-  SSNPmale_SW[1,1]/(SSNPmale_SW[2,1]+SSNPmale_SW[1,1])

QstmaleSW <- as.data.frame(c(DtGmale_SW, Dt10Gmale_SW, PSmale_SW, LLmale_SW, DtPmale_SW, NPmale_SW ), c("DtG","Dt10G","PS","LL","DtP","NP"))

###########################################################################
#####writing male outputs


### DtG

malePairwiseQstsDtGA <- as.data.frame(c(0,QstmaleAD[1,],QstmaleAE[1,],QstmaleAN[1,],QstmaleAR[1,],QstmaleAS[1,],QstmaleAW[1,]))
malePairwiseQstsDtGD <- as.data.frame(c(QstmaleAD[1,],0,QstmaleDE[1,],QstmaleDN[1,],QstmaleDR[1,],QstmaleDS[1,],QstmaleDW[1,]))
malePairwiseQstsDtGE <- as.data.frame(c(QstmaleAE[1,],QstmaleDE[1,],0,QstmaleEN[1,],QstmaleER[1,],QstmaleES[1,],QstmaleEW[1,]))
malePairwiseQstsDtGN <- as.data.frame(c(QstmaleAN[1,],QstmaleDN[1,],QstmaleEN[1,],0,QstmaleNR[1,],QstmaleNS[1,],QstmaleNW[1,]))
malePairwiseQstsDtGR <- as.data.frame(c(QstmaleAR[1,],QstmaleDR[1,],QstmaleER[1,],QstmaleNR[1,],0,QstmaleRS[1,],QstmaleRW[1,]))
malePairwiseQstsDtGS <- as.data.frame(c(QstmaleAS[1,],QstmaleDS[1,],QstmaleES[1,],QstmaleNS[1,],QstmaleRS[1,],0,QstmaleSW[1,]))
malePairwiseQstsDtGW <- as.data.frame(c(QstmaleAW[1,],QstmaleDW[1,],QstmaleEW[1,],QstmaleNW[1,],QstmaleRW[1,],QstmaleSW[1,],0))


malePairwiseALLDtG <- as.data.frame(cbind(malePairwiseQstsDtGA,malePairwiseQstsDtGD,malePairwiseQstsDtGE,
                                            malePairwiseQstsDtGN,malePairwiseQstsDtGR,malePairwiseQstsDtGS,
                                            malePairwiseQstsDtGW))
colnames(malePairwiseALLDtG) <- c("A","D","E","R","N","S","W")
rownames(malePairwiseALLDtG) <- c("A","D","E","R","N","S","W")



### Dt10G

malePairwiseQstsDt10GA <- as.data.frame(c(0,QstmaleAD[2,],QstmaleAE[2,],QstmaleAN[2,],QstmaleAR[2,],QstmaleAS[2,],QstmaleAW[2,]))
malePairwiseQstsDt10GD <- as.data.frame(c(QstmaleAD[2,],0,QstmaleDE[2,],QstmaleDN[2,],QstmaleDR[2,],QstmaleDS[2,],QstmaleDW[2,]))
malePairwiseQstsDt10GE <- as.data.frame(c(QstmaleAE[2,],QstmaleDE[2,],0,QstmaleEN[2,],QstmaleER[2,],QstmaleES[2,],QstmaleEW[2,]))
malePairwiseQstsDt10GN <- as.data.frame(c(QstmaleAN[2,],QstmaleDN[2,],QstmaleEN[2,],0,QstmaleNR[2,],QstmaleNS[2,],QstmaleNW[2,]))
malePairwiseQstsDt10GR <- as.data.frame(c(QstmaleAR[2,],QstmaleDR[2,],QstmaleER[2,],QstmaleNR[2,],0,QstmaleRS[2,],QstmaleRW[2,]))
malePairwiseQstsDt10GS <- as.data.frame(c(QstmaleAS[2,],QstmaleDS[2,],QstmaleES[2,],QstmaleNS[2,],QstmaleRS[2,],0,QstmaleSW[2,]))
malePairwiseQstsDt10GW <- as.data.frame(c(QstmaleAW[2,],QstmaleDW[2,],QstmaleEW[2,],QstmaleNW[2,],QstmaleRW[2,],QstmaleSW[2,],0))


malePairwiseALLDt10G <- as.data.frame(cbind(malePairwiseQstsDt10GA,malePairwiseQstsDt10GD,malePairwiseQstsDt10GE,
                                              malePairwiseQstsDt10GN,malePairwiseQstsDt10GR,malePairwiseQstsDt10GS,
                                              malePairwiseQstsDt10GW))
colnames(malePairwiseALLDt10G) <- c("A","D","E","R","N","S","W")
rownames(malePairwiseALLDt10G) <- c("A","D","E","R","N","S","W")



### PS

malePairwiseQstsPSA <- as.data.frame(c(0,QstmaleAD[3,],QstmaleAE[3,],QstmaleAN[3,],QstmaleAR[3,],QstmaleAS[3,],QstmaleAW[3,]))
malePairwiseQstsPSD <- as.data.frame(c(QstmaleAD[3,],0,QstmaleDE[3,],QstmaleDN[3,],QstmaleDR[3,],QstmaleDS[3,],QstmaleDW[3,]))
malePairwiseQstsPSE <- as.data.frame(c(QstmaleAE[3,],QstmaleDE[3,],0,QstmaleEN[3,],QstmaleER[3,],QstmaleES[3,],QstmaleEW[3,]))
malePairwiseQstsPSN <- as.data.frame(c(QstmaleAN[3,],QstmaleDN[3,],QstmaleEN[3,],0,QstmaleNR[3,],QstmaleNS[3,],QstmaleNW[3,]))
malePairwiseQstsPSR <- as.data.frame(c(QstmaleAR[3,],QstmaleDR[3,],QstmaleER[3,],QstmaleNR[3,],0,QstmaleRS[3,],QstmaleRW[3,]))
malePairwiseQstsPSS <- as.data.frame(c(QstmaleAS[3,],QstmaleDS[3,],QstmaleES[3,],QstmaleNS[3,],QstmaleRS[3,],0,QstmaleSW[3,]))
malePairwiseQstsPSW <- as.data.frame(c(QstmaleAW[3,],QstmaleDW[3,],QstmaleEW[3,],QstmaleNW[3,],QstmaleRW[3,],QstmaleSW[3,],0))


malePairwiseALLPS <- as.data.frame(cbind(malePairwiseQstsPSA,malePairwiseQstsPSD,malePairwiseQstsPSE,
                                           malePairwiseQstsPSN,malePairwiseQstsPSR,malePairwiseQstsPSS,
                                           malePairwiseQstsPSW))
colnames(malePairwiseALLPS) <- c("A","D","E","R","N","S","W")
rownames(malePairwiseALLPS) <- c("A","D","E","R","N","S","W")


### LL

malePairwiseQstsLLA <- as.data.frame(c(0,QstmaleAD[4,],QstmaleAE[4,],QstmaleAN[4,],QstmaleAR[4,],QstmaleAS[4,],QstmaleAW[4,]))
malePairwiseQstsLLD <- as.data.frame(c(QstmaleAD[4,],0,QstmaleDE[4,],QstmaleDN[4,],QstmaleDR[4,],QstmaleDS[4,],QstmaleDW[4,]))
malePairwiseQstsLLE <- as.data.frame(c(QstmaleAE[4,],QstmaleDE[4,],0,QstmaleEN[4,],QstmaleER[4,],QstmaleES[4,],QstmaleEW[4,]))
malePairwiseQstsLLN <- as.data.frame(c(QstmaleAN[4,],QstmaleDN[4,],QstmaleEN[4,],0,QstmaleNR[4,],QstmaleNS[4,],QstmaleNW[4,]))
malePairwiseQstsLLR <- as.data.frame(c(QstmaleAR[4,],QstmaleDR[4,],QstmaleER[4,],QstmaleNR[4,],0,QstmaleRS[4,],QstmaleRW[4,]))
malePairwiseQstsLLS <- as.data.frame(c(QstmaleAS[4,],QstmaleDS[4,],QstmaleES[4,],QstmaleNS[4,],QstmaleRS[4,],0,QstmaleSW[4,]))
malePairwiseQstsLLW <- as.data.frame(c(QstmaleAW[4,],QstmaleDW[4,],QstmaleEW[4,],QstmaleNW[4,],QstmaleRW[4,],QstmaleSW[4,],0))


malePairwiseALLLL <- as.data.frame(cbind(malePairwiseQstsLLA,malePairwiseQstsLLD,malePairwiseQstsLLE,
                                         malePairwiseQstsLLN,malePairwiseQstsLLR,malePairwiseQstsLLS,
                                         malePairwiseQstsLLW))
colnames(malePairwiseALLLL) <- c("A","D","E","R","N","S","W")
rownames(malePairwiseALLLL) <- c("A","D","E","R","N","S","W")

### DtP

malePairwiseQstsDtPA <- as.data.frame(c(0,QstmaleAD[4,],QstmaleAE[4,],QstmaleAN[4,],QstmaleAR[4,],QstmaleAS[4,],QstmaleAW[4,]))
malePairwiseQstsDtPD <- as.data.frame(c(QstmaleAD[4,],0,QstmaleDE[4,],QstmaleDN[4,],QstmaleDR[4,],QstmaleDS[4,],QstmaleDW[4,]))
malePairwiseQstsDtPE <- as.data.frame(c(QstmaleAE[4,],QstmaleDE[4,],0,QstmaleEN[4,],QstmaleER[4,],QstmaleES[4,],QstmaleEW[4,]))
malePairwiseQstsDtPN <- as.data.frame(c(QstmaleAN[4,],QstmaleDN[4,],QstmaleEN[4,],0,QstmaleNR[4,],QstmaleNS[4,],QstmaleNW[4,]))
malePairwiseQstsDtPR <- as.data.frame(c(QstmaleAR[4,],QstmaleDR[4,],QstmaleER[4,],QstmaleNR[4,],0,QstmaleRS[4,],QstmaleRW[4,]))
malePairwiseQstsDtPS <- as.data.frame(c(QstmaleAS[4,],QstmaleDS[4,],QstmaleES[4,],QstmaleNS[4,],QstmaleRS[4,],0,QstmaleSW[4,]))
malePairwiseQstsDtPW <- as.data.frame(c(QstmaleAW[4,],QstmaleDW[4,],QstmaleEW[4,],QstmaleNW[4,],QstmaleRW[4,],QstmaleSW[4,],0))

malePairwiseALLDtP <- as.data.frame(cbind(malePairwiseQstsDtPA,malePairwiseQstsDtPD,malePairwiseQstsDtPE,
                                         malePairwiseQstsDtPN,malePairwiseQstsDtPR,malePairwiseQstsDtPS,
                                         malePairwiseQstsDtPW))
colnames(malePairwiseALLDtP) <- c("A","D","E","R","N","S","W")
rownames(malePairwiseALLDtP) <- c("A","D","E","R","N","S","W")


### NP

malePairwiseQstsNPA <- as.data.frame(c(0,QstmaleAD[4,],QstmaleAE[4,],QstmaleAN[4,],QstmaleAR[4,],QstmaleAS[4,],QstmaleAW[4,]))
malePairwiseQstsNPD <- as.data.frame(c(QstmaleAD[4,],0,QstmaleDE[4,],QstmaleDN[4,],QstmaleDR[4,],QstmaleDS[4,],QstmaleDW[4,]))
malePairwiseQstsNPE <- as.data.frame(c(QstmaleAE[4,],QstmaleDE[4,],0,QstmaleEN[4,],QstmaleER[4,],QstmaleES[4,],QstmaleEW[4,]))
malePairwiseQstsNPN <- as.data.frame(c(QstmaleAN[4,],QstmaleDN[4,],QstmaleEN[4,],0,QstmaleNR[4,],QstmaleNS[4,],QstmaleNW[4,]))
malePairwiseQstsNPR <- as.data.frame(c(QstmaleAR[4,],QstmaleDR[4,],QstmaleER[4,],QstmaleNR[4,],0,QstmaleRS[4,],QstmaleRW[4,]))
malePairwiseQstsNPS <- as.data.frame(c(QstmaleAS[4,],QstmaleDS[4,],QstmaleES[4,],QstmaleNS[4,],QstmaleRS[4,],0,QstmaleSW[4,]))
malePairwiseQstsNPW <- as.data.frame(c(QstmaleAW[4,],QstmaleDW[4,],QstmaleEW[4,],QstmaleNW[4,],QstmaleRW[4,],QstmaleSW[4,],0))



malePairwiseALLNP <- as.data.frame(cbind(malePairwiseQstsNPA,malePairwiseQstsNPD,malePairwiseQstsNPE,
                                           malePairwiseQstsNPN,malePairwiseQstsNPR,malePairwiseQstsNPS,
                                           malePairwiseQstsNPW))
colnames(malePairwiseALLNP) <- c("A","D","E","R","N","S","W")
rownames(malePairwiseALLNP) <- c("A","D","E","R","N","S","W")



traitNames <- c("DtG","","","","","","",
                "Dt10G","","","","","","",
                "PS","","","","","","",
                "LL","","","","","","",
                "DtP","","","","","","",
                "NP","","","","","","")


malePairwiseQst <- as.data.frame(rbind(malePairwiseALLDtG,malePairwiseALLDt10G,malePairwiseALLPS,malePairwiseALLLL,malePairwiseALLDtP,malePairwiseALLNP))
malePairwiseQstwTrait <- as.data.frame(cbind(traitNames,malePairwiseQst))


write.csv(malePairwiseQstwTrait,"tables/malePairwiseQstwTrait.csv")


##############################################################


########## Between population variances (sum of squares)

### DtG

maleBetweenPopSSDtGA <- as.data.frame(c(0,SSDtGmale_AD[1,1],SSDtGmale_AE[1,1],SSDtGmale_AN[1,1],SSDtGmale_AR[1,1],SSDtGmale_AS[1,1],SSDtGmale_AW[1,1]))
maleBetweenPopSSDtGD <- as.data.frame(c(SSDtGmale_AD[1,1],0,SSDtGmale_DE[1,1],SSDtGmale_DN[1,1],SSDtGmale_DR[1,1],SSDtGmale_DS[1,1],SSDtGmale_DW[1,1]))
maleBetweenPopSSDtGE <- as.data.frame(c(SSDtGmale_AE[1,1],SSDtGmale_DE[1,1],0,SSDtGmale_EN[1,1],SSDtGmale_ER[1,1],SSDtGmale_ES[1,1],SSDtGmale_EW[1,1]))
maleBetweenPopSSDtGN <- as.data.frame(c(SSDtGmale_AN[1,1],SSDtGmale_DN[1,1],SSDtGmale_EN[1,1],0,SSDtGmale_NR[1,1],SSDtGmale_NS[1,1],SSDtGmale_NW[1,1]))
maleBetweenPopSSDtGR <- as.data.frame(c(SSDtGmale_AR[1,1],SSDtGmale_DR[1,1],SSDtGmale_ER[1,1],SSDtGmale_NR[1,1],0,SSDtGmale_RS[1,1],SSDtGmale_RW[1,1]))
maleBetweenPopSSDtGS <- as.data.frame(c(SSDtGmale_AS[1,1],SSDtGmale_DS[1,1],SSDtGmale_ES[1,1],SSDtGmale_NS[1,1],SSDtGmale_RS[1,1],0,SSDtGmale_SW[1,1]))
maleBetweenPopSSDtGW <- as.data.frame(c(SSDtGmale_AW[1,1],SSDtGmale_DW[1,1],SSDtGmale_EW[1,1],SSDtGmale_NW[1,1],SSDtGmale_RW[1,1],SSDtGmale_SW[1,1],0))

maleBetweenPopSSALLDtG <- as.data.frame(cbind(maleBetweenPopSSDtGA,maleBetweenPopSSDtGD,maleBetweenPopSSDtGE,
                                          maleBetweenPopSSDtGN,maleBetweenPopSSDtGR,maleBetweenPopSSDtGS,
                                          maleBetweenPopSSDtGW))
colnames(maleBetweenPopSSALLDtG) <- c("A","D","E","R","N","S","W")
rownames(maleBetweenPopSSALLDtG) <- c("A","D","E","R","N","S","W")



### Dt10G

maleBetweenPopSSDt10GA <- as.data.frame(c(0,SSDt10Gmale_AD[1,1],SSDt10Gmale_AE[1,1],SSDt10Gmale_AN[1,1],SSDt10Gmale_AR[1,1],SSDt10Gmale_AS[1,1],SSDt10Gmale_AW[1,1]))
maleBetweenPopSSDt10GD <- as.data.frame(c(SSDt10Gmale_AD[1,1],0,SSDt10Gmale_DE[1,1],SSDt10Gmale_DN[1,1],SSDt10Gmale_DR[1,1],SSDt10Gmale_DS[1,1],SSDt10Gmale_DW[1,1]))
maleBetweenPopSSDt10GE <- as.data.frame(c(SSDt10Gmale_AE[1,1],SSDt10Gmale_DE[1,1],0,SSDt10Gmale_EN[1,1],SSDt10Gmale_ER[1,1],SSDt10Gmale_ES[1,1],SSDt10Gmale_EW[1,1]))
maleBetweenPopSSDt10GN <- as.data.frame(c(SSDt10Gmale_AN[1,1],SSDt10Gmale_DN[1,1],SSDt10Gmale_EN[1,1],0,SSDt10Gmale_NR[1,1],SSDt10Gmale_NS[1,1],SSDt10Gmale_NW[1,1]))
maleBetweenPopSSDt10GR <- as.data.frame(c(SSDt10Gmale_AR[1,1],SSDt10Gmale_DR[1,1],SSDt10Gmale_ER[1,1],SSDt10Gmale_NR[1,1],0,SSDt10Gmale_RS[1,1],SSDt10Gmale_RW[1,1]))
maleBetweenPopSSDt10GS <- as.data.frame(c(SSDt10Gmale_AS[1,1],SSDt10Gmale_DS[1,1],SSDt10Gmale_ES[1,1],SSDt10Gmale_NS[1,1],SSDt10Gmale_RS[1,1],0,SSDt10Gmale_SW[1,1]))
maleBetweenPopSSDt10GW <- as.data.frame(c(SSDt10Gmale_AW[1,1],SSDt10Gmale_DW[1,1],SSDt10Gmale_EW[1,1],SSDt10Gmale_NW[1,1],SSDt10Gmale_RW[1,1],SSDt10Gmale_SW[1,1],0))

maleBetweenPopSSALLDt10G <- as.data.frame(cbind(maleBetweenPopSSDt10GA,maleBetweenPopSSDt10GD,maleBetweenPopSSDt10GE,
                                            maleBetweenPopSSDt10GN,maleBetweenPopSSDt10GR,maleBetweenPopSSDt10GS,
                                            maleBetweenPopSSDt10GW))
colnames(maleBetweenPopSSALLDt10G) <- c("A","D","E","R","N","S","W")
rownames(maleBetweenPopSSALLDt10G) <- c("A","D","E","R","N","S","W")



### PS

maleBetweenPopSSPSA <- as.data.frame(c(0,SSPSmale_AD[1,1],SSPSmale_AE[1,1],SSPSmale_AN[1,1],SSPSmale_AR[1,1],SSPSmale_AS[1,1],SSPSmale_AW[1,1]))
maleBetweenPopSSPSD <- as.data.frame(c(SSPSmale_AD[1,1],0,SSPSmale_DE[1,1],SSPSmale_DN[1,1],SSPSmale_DR[1,1],SSPSmale_DS[1,1],SSPSmale_DW[1,1]))
maleBetweenPopSSPSE <- as.data.frame(c(SSPSmale_AE[1,1],SSPSmale_DE[1,1],0,SSPSmale_EN[1,1],SSPSmale_ER[1,1],SSPSmale_ES[1,1],SSPSmale_EW[1,1]))
maleBetweenPopSSPSN <- as.data.frame(c(SSPSmale_AN[1,1],SSPSmale_DN[1,1],SSPSmale_EN[1,1],0,SSPSmale_NR[1,1],SSPSmale_NS[1,1],SSPSmale_NW[1,1]))
maleBetweenPopSSPSR <- as.data.frame(c(SSPSmale_AR[1,1],SSPSmale_DR[1,1],SSPSmale_ER[1,1],SSPSmale_NR[1,1],0,SSPSmale_RS[1,1],SSPSmale_RW[1,1]))
maleBetweenPopSSPSS <- as.data.frame(c(SSPSmale_AS[1,1],SSPSmale_DS[1,1],SSPSmale_ES[1,1],SSPSmale_NS[1,1],SSPSmale_RS[1,1],0,SSPSmale_SW[1,1]))
maleBetweenPopSSPSW <- as.data.frame(c(SSPSmale_AW[1,1],SSPSmale_DW[1,1],SSPSmale_EW[1,1],SSPSmale_NW[1,1],SSPSmale_RW[1,1],SSPSmale_SW[1,1],0))

maleBetweenPopSSALLPS <- as.data.frame(cbind(maleBetweenPopSSPSA,maleBetweenPopSSPSD,maleBetweenPopSSPSE,
                                         maleBetweenPopSSPSN,maleBetweenPopSSPSR,maleBetweenPopSSPSS,
                                         maleBetweenPopSSPSW))
colnames(maleBetweenPopSSALLPS) <- c("A","D","E","R","N","S","W")
rownames(maleBetweenPopSSALLPS) <- c("A","D","E","R","N","S","W")




### LL

maleBetweenPopSSLLA <- as.data.frame(c(0,SSLLmale_AD[1,1],SSLLmale_AE[1,1],SSLLmale_AN[1,1],SSLLmale_AR[1,1],SSLLmale_AS[1,1],SSLLmale_AW[1,1]))
maleBetweenPopSSLLD <- as.data.frame(c(SSLLmale_AD[1,1],0,SSLLmale_DE[1,1],SSLLmale_DN[1,1],SSLLmale_DR[1,1],SSLLmale_DS[1,1],SSLLmale_DW[1,1]))
maleBetweenPopSSLLE <- as.data.frame(c(SSLLmale_AE[1,1],SSLLmale_DE[1,1],0,SSLLmale_EN[1,1],SSLLmale_ER[1,1],SSLLmale_ES[1,1],SSLLmale_EW[1,1]))
maleBetweenPopSSLLN <- as.data.frame(c(SSLLmale_AN[1,1],SSLLmale_DN[1,1],SSLLmale_EN[1,1],0,SSLLmale_NR[1,1],SSLLmale_NS[1,1],SSLLmale_NW[1,1]))
maleBetweenPopSSLLR <- as.data.frame(c(SSLLmale_AR[1,1],SSLLmale_DR[1,1],SSLLmale_ER[1,1],SSLLmale_NR[1,1],0,SSLLmale_RS[1,1],SSLLmale_RW[1,1]))
maleBetweenPopSSLLS <- as.data.frame(c(SSLLmale_AS[1,1],SSLLmale_DS[1,1],SSLLmale_ES[1,1],SSLLmale_NS[1,1],SSLLmale_RS[1,1],0,SSLLmale_SW[1,1]))
maleBetweenPopSSLLW <- as.data.frame(c(SSLLmale_AW[1,1],SSLLmale_DW[1,1],SSLLmale_EW[1,1],SSLLmale_NW[1,1],SSLLmale_RW[1,1],SSLLmale_SW[1,1],0))



maleBetweenPopSSALLLL <- as.data.frame(cbind(maleBetweenPopSSLLA,maleBetweenPopSSLLD,maleBetweenPopSSLLE,
                                         maleBetweenPopSSLLN,maleBetweenPopSSLLR,maleBetweenPopSSLLS,
                                         maleBetweenPopSSLLW))
colnames(maleBetweenPopSSALLLL) <- c("A","D","E","R","N","S","W")
rownames(maleBetweenPopSSALLLL) <- c("A","D","E","R","N","S","W")



### DtP

maleBetweenPopDtPA <- as.data.frame(c(0,QstmaleAD[4,],QstmaleAE[4,],QstmaleAN[4,],QstmaleAR[4,],QstmaleAS[4,],QstmaleAW[4,]))
maleBetweenPopDtPD <- as.data.frame(c(QstmaleAD[4,],0,QstmaleDE[4,],QstmaleDN[4,],QstmaleDR[4,],QstmaleDS[4,],QstmaleDW[4,]))
maleBetweenPopDtPE <- as.data.frame(c(QstmaleAE[4,],QstmaleDE[4,],0,QstmaleEN[4,],QstmaleER[4,],QstmaleES[4,],QstmaleEW[4,]))
maleBetweenPopDtPN <- as.data.frame(c(QstmaleAN[4,],QstmaleDN[4,],QstmaleEN[4,],0,QstmaleNR[4,],QstmaleNS[4,],QstmaleNW[4,]))
maleBetweenPopDtPR <- as.data.frame(c(QstmaleAR[4,],QstmaleDR[4,],QstmaleER[4,],QstmaleNR[4,],0,QstmaleRS[4,],QstmaleRW[4,]))
maleBetweenPopDtPS <- as.data.frame(c(QstmaleAS[4,],QstmaleDS[4,],QstmaleES[4,],QstmaleNS[4,],QstmaleRS[4,],0,QstmaleSW[4,]))
maleBetweenPopDtPW <- as.data.frame(c(QstmaleAW[4,],QstmaleDW[4,],QstmaleEW[4,],QstmaleNW[4,],QstmaleRW[4,],QstmaleSW[4,],0))

maleBetweenPopALLDtP <- as.data.frame(cbind(maleBetweenPopDtPA,maleBetweenPopDtPD,maleBetweenPopDtPE,
                                          maleBetweenPopDtPN,maleBetweenPopDtPR,maleBetweenPopDtPS,
                                          maleBetweenPopDtPW))
colnames(maleBetweenPopALLDtP) <- c("A","D","E","R","N","S","W")
rownames(maleBetweenPopALLDtP) <- c("A","D","E","R","N","S","W")


### NP

maleBetweenPopNPA <- as.data.frame(c(0,QstmaleAD[4,],QstmaleAE[4,],QstmaleAN[4,],QstmaleAR[4,],QstmaleAS[4,],QstmaleAW[4,]))
maleBetweenPopNPD <- as.data.frame(c(QstmaleAD[4,],0,QstmaleDE[4,],QstmaleDN[4,],QstmaleDR[4,],QstmaleDS[4,],QstmaleDW[4,]))
maleBetweenPopNPE <- as.data.frame(c(QstmaleAE[4,],QstmaleDE[4,],0,QstmaleEN[4,],QstmaleER[4,],QstmaleES[4,],QstmaleEW[4,]))
maleBetweenPopNPN <- as.data.frame(c(QstmaleAN[4,],QstmaleDN[4,],QstmaleEN[4,],0,QstmaleNR[4,],QstmaleNS[4,],QstmaleNW[4,]))
maleBetweenPopNPR <- as.data.frame(c(QstmaleAR[4,],QstmaleDR[4,],QstmaleER[4,],QstmaleNR[4,],0,QstmaleRS[4,],QstmaleRW[4,]))
maleBetweenPopNPS <- as.data.frame(c(QstmaleAS[4,],QstmaleDS[4,],QstmaleES[4,],QstmaleNS[4,],QstmaleRS[4,],0,QstmaleSW[4,]))
maleBetweenPopNPW <- as.data.frame(c(QstmaleAW[4,],QstmaleDW[4,],QstmaleEW[4,],QstmaleNW[4,],QstmaleRW[4,],QstmaleSW[4,],0))



maleBetweenPopALLNP <- as.data.frame(cbind(maleBetweenPopNPA,maleBetweenPopNPD,maleBetweenPopNPE,
                                         maleBetweenPopNPN,maleBetweenPopNPR,maleBetweenPopNPS,
                                         maleBetweenPopNPW))
colnames(maleBetweenPopALLNP) <- c("A","D","E","R","N","S","W")
rownames(maleBetweenPopALLNP) <- c("A","D","E","R","N","S","W")



maleBetweenPopSumofSquares <- as.data.frame(rbind(maleBetweenPopSSALLDtG,maleBetweenPopSSALLDt10G,maleBetweenPopSSALLPS,maleBetweenPopSSALLLL,maleBetweenPopALLDtP,maleBetweenPopALLNP))
maleBetweenPopSumofSquares <- as.data.frame(cbind(traitNames,maleBetweenPopSumofSquares))


write.csv(maleBetweenPopSumofSquares,"tables/maleBetweenPopSumofSquares.csv")



##############################################################


########## within population variances (sum of squares)

### DtG

maleWithinPopSSDtGA <- as.data.frame(c(0,SSDtGmale_AD[2,1],SSDtGmale_AE[2,1],SSDtGmale_AN[2,1],SSDtGmale_AR[2,1],SSDtGmale_AS[2,1],SSDtGmale_AW[2,1]))
maleWithinPopSSDtGD <- as.data.frame(c(SSDtGmale_AD[2,1],0,SSDtGmale_DE[2,1],SSDtGmale_DN[2,1],SSDtGmale_DR[2,1],SSDtGmale_DS[2,1],SSDtGmale_DW[2,1]))
maleWithinPopSSDtGE <- as.data.frame(c(SSDtGmale_AE[2,1],SSDtGmale_DE[2,1],0,SSDtGmale_EN[2,1],SSDtGmale_ER[2,1],SSDtGmale_ES[2,1],SSDtGmale_EW[2,1]))
maleWithinPopSSDtGN <- as.data.frame(c(SSDtGmale_AN[2,1],SSDtGmale_DN[2,1],SSDtGmale_EN[2,1],0,SSDtGmale_NR[2,1],SSDtGmale_NS[2,1],SSDtGmale_NW[2,1]))
maleWithinPopSSDtGR <- as.data.frame(c(SSDtGmale_AR[2,1],SSDtGmale_DR[2,1],SSDtGmale_ER[2,1],SSDtGmale_NR[2,1],0,SSDtGmale_RS[2,1],SSDtGmale_RW[2,1]))
maleWithinPopSSDtGS <- as.data.frame(c(SSDtGmale_AS[2,1],SSDtGmale_DS[2,1],SSDtGmale_ES[2,1],SSDtGmale_NS[2,1],SSDtGmale_RS[2,1],0,SSDtGmale_SW[2,1]))
maleWithinPopSSDtGW <- as.data.frame(c(SSDtGmale_AW[2,1],SSDtGmale_DW[2,1],SSDtGmale_EW[2,1],SSDtGmale_NW[2,1],SSDtGmale_RW[2,1],SSDtGmale_SW[2,1],0))

maleWithinPopSSALLDtG <- as.data.frame(cbind(maleWithinPopSSDtGA,maleWithinPopSSDtGD,maleWithinPopSSDtGE,
                                             maleWithinPopSSDtGN,maleWithinPopSSDtGR,maleWithinPopSSDtGS,
                                             maleWithinPopSSDtGW))
colnames(maleWithinPopSSALLDtG) <- c("A","D","E","R","N","S","W")
rownames(maleWithinPopSSALLDtG) <- c("A","D","E","R","N","S","W")



### Dt10G

maleWithinPopSSDt10GA <- as.data.frame(c(0,SSDt10Gmale_AD[2,1],SSDt10Gmale_AE[2,1],SSDt10Gmale_AN[2,1],SSDt10Gmale_AR[2,1],SSDt10Gmale_AS[2,1],SSDt10Gmale_AW[2,1]))
maleWithinPopSSDt10GD <- as.data.frame(c(SSDt10Gmale_AD[2,1],0,SSDt10Gmale_DE[2,1],SSDt10Gmale_DN[2,1],SSDt10Gmale_DR[2,1],SSDt10Gmale_DS[2,1],SSDt10Gmale_DW[2,1]))
maleWithinPopSSDt10GE <- as.data.frame(c(SSDt10Gmale_AE[2,1],SSDt10Gmale_DE[2,1],0,SSDt10Gmale_EN[2,1],SSDt10Gmale_ER[2,1],SSDt10Gmale_ES[2,1],SSDt10Gmale_EW[2,1]))
maleWithinPopSSDt10GN <- as.data.frame(c(SSDt10Gmale_AN[2,1],SSDt10Gmale_DN[2,1],SSDt10Gmale_EN[2,1],0,SSDt10Gmale_NR[2,1],SSDt10Gmale_NS[2,1],SSDt10Gmale_NW[2,1]))
maleWithinPopSSDt10GR <- as.data.frame(c(SSDt10Gmale_AR[2,1],SSDt10Gmale_DR[2,1],SSDt10Gmale_ER[2,1],SSDt10Gmale_NR[2,1],0,SSDt10Gmale_RS[2,1],SSDt10Gmale_RW[2,1]))
maleWithinPopSSDt10GS <- as.data.frame(c(SSDt10Gmale_AS[2,1],SSDt10Gmale_DS[2,1],SSDt10Gmale_ES[2,1],SSDt10Gmale_NS[2,1],SSDt10Gmale_RS[2,1],0,SSDt10Gmale_SW[2,1]))
maleWithinPopSSDt10GW <- as.data.frame(c(SSDt10Gmale_AW[2,1],SSDt10Gmale_DW[2,1],SSDt10Gmale_EW[2,1],SSDt10Gmale_NW[2,1],SSDt10Gmale_RW[2,1],SSDt10Gmale_SW[2,1],0))

maleWithinPopSSALLDt10G <- as.data.frame(cbind(maleWithinPopSSDt10GA,maleWithinPopSSDt10GD,maleWithinPopSSDt10GE,
                                               maleWithinPopSSDt10GN,maleWithinPopSSDt10GR,maleWithinPopSSDt10GS,
                                               maleWithinPopSSDt10GW))
colnames(maleWithinPopSSALLDt10G) <- c("A","D","E","R","N","S","W")
rownames(maleWithinPopSSALLDt10G) <- c("A","D","E","R","N","S","W")



### PS

maleWithinPopSSPSA <- as.data.frame(c(0,SSPSmale_AD[2,1],SSPSmale_AE[2,1],SSPSmale_AN[2,1],SSPSmale_AR[2,1],SSPSmale_AS[2,1],SSPSmale_AW[2,1]))
maleWithinPopSSPSD <- as.data.frame(c(SSPSmale_AD[2,1],0,SSPSmale_DE[2,1],SSPSmale_DN[2,1],SSPSmale_DR[2,1],SSPSmale_DS[2,1],SSPSmale_DW[2,1]))
maleWithinPopSSPSE <- as.data.frame(c(SSPSmale_AE[2,1],SSPSmale_DE[2,1],0,SSPSmale_EN[2,1],SSPSmale_ER[2,1],SSPSmale_ES[2,1],SSPSmale_EW[2,1]))
maleWithinPopSSPSN <- as.data.frame(c(SSPSmale_AN[2,1],SSPSmale_DN[2,1],SSPSmale_EN[2,1],0,SSPSmale_NR[2,1],SSPSmale_NS[2,1],SSPSmale_NW[2,1]))
maleWithinPopSSPSR <- as.data.frame(c(SSPSmale_AR[2,1],SSPSmale_DR[2,1],SSPSmale_ER[2,1],SSPSmale_NR[2,1],0,SSPSmale_RS[2,1],SSPSmale_RW[2,1]))
maleWithinPopSSPSS <- as.data.frame(c(SSPSmale_AS[2,1],SSPSmale_DS[2,1],SSPSmale_ES[2,1],SSPSmale_NS[2,1],SSPSmale_RS[2,1],0,SSPSmale_SW[2,1]))
maleWithinPopSSPSW <- as.data.frame(c(SSPSmale_AW[2,1],SSPSmale_DW[2,1],SSPSmale_EW[2,1],SSPSmale_NW[2,1],SSPSmale_RW[2,1],SSPSmale_SW[2,1],0))

maleWithinPopSSALLPS <- as.data.frame(cbind(maleWithinPopSSPSA,maleWithinPopSSPSD,maleWithinPopSSPSE,
                                            maleWithinPopSSPSN,maleWithinPopSSPSR,maleWithinPopSSPSS,
                                            maleWithinPopSSPSW))
colnames(maleWithinPopSSALLPS) <- c("A","D","E","R","N","S","W")
rownames(maleWithinPopSSALLPS) <- c("A","D","E","R","N","S","W")




### LL

maleWithinPopSSLLA <- as.data.frame(c(0,SSLLmale_AD[2,1],SSLLmale_AE[2,1],SSLLmale_AN[2,1],SSLLmale_AR[2,1],SSLLmale_AS[2,1],SSLLmale_AW[2,1]))
maleWithinPopSSLLD <- as.data.frame(c(SSLLmale_AD[2,1],0,SSLLmale_DE[2,1],SSLLmale_DN[2,1],SSLLmale_DR[2,1],SSLLmale_DS[2,1],SSLLmale_DW[2,1]))
maleWithinPopSSLLE <- as.data.frame(c(SSLLmale_AE[2,1],SSLLmale_DE[2,1],0,SSLLmale_EN[2,1],SSLLmale_ER[2,1],SSLLmale_ES[2,1],SSLLmale_EW[2,1]))
maleWithinPopSSLLN <- as.data.frame(c(SSLLmale_AN[2,1],SSLLmale_DN[2,1],SSLLmale_EN[2,1],0,SSLLmale_NR[2,1],SSLLmale_NS[2,1],SSLLmale_NW[2,1]))
maleWithinPopSSLLR <- as.data.frame(c(SSLLmale_AR[2,1],SSLLmale_DR[2,1],SSLLmale_ER[2,1],SSLLmale_NR[2,1],0,SSLLmale_RS[2,1],SSLLmale_RW[2,1]))
maleWithinPopSSLLS <- as.data.frame(c(SSLLmale_AS[2,1],SSLLmale_DS[2,1],SSLLmale_ES[2,1],SSLLmale_NS[2,1],SSLLmale_RS[2,1],0,SSLLmale_SW[2,1]))
maleWithinPopSSLLW <- as.data.frame(c(SSLLmale_AW[2,1],SSLLmale_DW[2,1],SSLLmale_EW[2,1],SSLLmale_NW[2,1],SSLLmale_RW[2,1],SSLLmale_SW[2,1],0))



maleWithinPopSSALLLL <- as.data.frame(cbind(maleWithinPopSSLLA,maleWithinPopSSLLD,maleWithinPopSSLLE,
                                            maleWithinPopSSLLN,maleWithinPopSSLLR,maleWithinPopSSLLS,
                                            maleWithinPopSSLLW))
colnames(maleWithinPopSSALLLL) <- c("A","D","E","R","N","S","W")
rownames(maleWithinPopSSALLLL) <- c("A","D","E","R","N","S","W")



### DtP

maleWithinPopDtPA <- as.data.frame(c(0,QstmaleAD[4,],QstmaleAE[4,],QstmaleAN[4,],QstmaleAR[4,],QstmaleAS[4,],QstmaleAW[4,]))
maleWithinPopDtPD <- as.data.frame(c(QstmaleAD[4,],0,QstmaleDE[4,],QstmaleDN[4,],QstmaleDR[4,],QstmaleDS[4,],QstmaleDW[4,]))
maleWithinPopDtPE <- as.data.frame(c(QstmaleAE[4,],QstmaleDE[4,],0,QstmaleEN[4,],QstmaleER[4,],QstmaleES[4,],QstmaleEW[4,]))
maleWithinPopDtPN <- as.data.frame(c(QstmaleAN[4,],QstmaleDN[4,],QstmaleEN[4,],0,QstmaleNR[4,],QstmaleNS[4,],QstmaleNW[4,]))
maleWithinPopDtPR <- as.data.frame(c(QstmaleAR[4,],QstmaleDR[4,],QstmaleER[4,],QstmaleNR[4,],0,QstmaleRS[4,],QstmaleRW[4,]))
maleWithinPopDtPS <- as.data.frame(c(QstmaleAS[4,],QstmaleDS[4,],QstmaleES[4,],QstmaleNS[4,],QstmaleRS[4,],0,QstmaleSW[4,]))
maleWithinPopDtPW <- as.data.frame(c(QstmaleAW[4,],QstmaleDW[4,],QstmaleEW[4,],QstmaleNW[4,],QstmaleRW[4,],QstmaleSW[4,],0))

maleWithinPopALLDtP <- as.data.frame(cbind(maleWithinPopDtPA,maleWithinPopDtPD,maleWithinPopDtPE,
                                            maleWithinPopDtPN,maleWithinPopDtPR,maleWithinPopDtPS,
                                            maleWithinPopDtPW))
colnames(maleWithinPopALLDtP) <- c("A","D","E","R","N","S","W")
rownames(maleWithinPopALLDtP) <- c("A","D","E","R","N","S","W")


### NP

maleWithinPopNPA <- as.data.frame(c(0,QstmaleAD[4,],QstmaleAE[4,],QstmaleAN[4,],QstmaleAR[4,],QstmaleAS[4,],QstmaleAW[4,]))
maleWithinPopNPD <- as.data.frame(c(QstmaleAD[4,],0,QstmaleDE[4,],QstmaleDN[4,],QstmaleDR[4,],QstmaleDS[4,],QstmaleDW[4,]))
maleWithinPopNPE <- as.data.frame(c(QstmaleAE[4,],QstmaleDE[4,],0,QstmaleEN[4,],QstmaleER[4,],QstmaleES[4,],QstmaleEW[4,]))
maleWithinPopNPN <- as.data.frame(c(QstmaleAN[4,],QstmaleDN[4,],QstmaleEN[4,],0,QstmaleNR[4,],QstmaleNS[4,],QstmaleNW[4,]))
maleWithinPopNPR <- as.data.frame(c(QstmaleAR[4,],QstmaleDR[4,],QstmaleER[4,],QstmaleNR[4,],0,QstmaleRS[4,],QstmaleRW[4,]))
maleWithinPopNPS <- as.data.frame(c(QstmaleAS[4,],QstmaleDS[4,],QstmaleES[4,],QstmaleNS[4,],QstmaleRS[4,],0,QstmaleSW[4,]))
maleWithinPopNPW <- as.data.frame(c(QstmaleAW[4,],QstmaleDW[4,],QstmaleEW[4,],QstmaleNW[4,],QstmaleRW[4,],QstmaleSW[4,],0))



maleWithinPopALLNP <- as.data.frame(cbind(maleWithinPopNPA,maleWithinPopNPD,maleWithinPopNPE,
                                           maleWithinPopNPN,maleWithinPopNPR,maleWithinPopNPS,
                                           maleWithinPopNPW))
colnames(maleWithinPopALLNP) <- c("A","D","E","R","N","S","W")
rownames(maleWithinPopALLNP) <- c("A","D","E","R","N","S","W")




maleWithinPopSumofSquares <- as.data.frame(rbind(maleWithinPopSSALLDtG,maleWithinPopSSALLDt10G,maleWithinPopSSALLPS,maleWithinPopSSALLLL,maleWithinPopALLDtP,maleWithinPopALLNP))
maleWithinPopSumofSquares <- as.data.frame(cbind(traitNames,maleWithinPopSumofSquares))


write.csv(maleWithinPopSumofSquares,"tables/maleWithinPopSumofSquares.csv")










