
### ANOVA with sexes combined

## analyses for McDaniel et al done with R version 3.5.3

#### load in data
combined <- na.omit(read.csv("data/CombinedENA.csv"))
combined$Individual <- as.character(combined$Individual)
combined$Clone <- as.character(combined$Clone)

##### transform data
combinedtransDtG <- ((combined$DtG)^(-1.75))
combinedtransDt10G <- ((combined$Dt10G)^(-0.25))
combinedtransPS <- ((combined$PS)^(0.75))
combinedtransLL <- ((combined$LL)^(0.25))

### write new data frame for transformed data 
combinedPopulation <- as.character(combined$Population)
combinedIndividual <- as.character(combined$Individual)
combinedSex <- as.character(combined$Sex)
combinedTransALL <- as.data.frame(cbind(combinedtransDtG,combinedtransDt10G,combinedtransPS,combinedtransLL,
                                        combinedPopulation,combinedIndividual,combinedSex))
combinedTransALL$combinedtransDtG <- as.numeric(as.character(combinedTransALL$combinedtransDtG))
combinedTransALL$combinedtransDt10G <- as.numeric(as.character(combinedTransALL$combinedtransDt10G))
combinedTransALL$combinedtransPS <- as.numeric(as.character(combinedTransALL$combinedtransPS))
combinedTransALL$combinedtransLL <- as.numeric(as.character(combinedTransALL$combinedtransLL))


##########################################################################

### ANOVAs for all shared traits btwn males and females
aovALLDtGtrans <- aov(combinedtransDtG~combinedSex*combinedPopulation + 
                        combinedPopulation/combinedIndividual, data=combinedTransALL)
summary(aovALLDtGtrans)


aovALLDt10Gtrans <- aov(combinedtransDt10G~combinedSex*combinedPopulation +
                          combinedPopulation/combinedIndividual, data=combinedTransALL)
summary(aovALLDt10Gtrans)


aovALLPStrans <- aov(combinedtransPS~combinedSex*combinedPopulation + 
                       combinedPopulation/combinedIndividual, data=combinedTransALL)
summary(aovALLPStrans)


aovALLLLtrans <- aov(combinedtransLL~combinedSex*combinedPopulation +
                       combinedPopulation/combinedIndividual, data=combinedTransALL)
summary(aovALLLLtrans)


## MANOVA for shared traits btwn males and females
manovALLtrans <- manova(cbind(combinedtransDtG, combinedtransDt10G, combinedtransPS, combinedtransLL)~
                          combinedSex*combinedPopulation + combinedPopulation/combinedIndividual, data=combinedTransALL)
summary(manovALLtrans)

