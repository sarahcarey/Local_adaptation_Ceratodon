
#### pca plots

## analyses for McDaniel et al done with R version 3.5.3
## ggbiplot version 0.55

install_github("vqv/ggbiplot")
library(ggbiplot)

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

### averaging the traits for pca 

maleTransTraitsAvg <- aggregate(maleTransALL, by=list(maleTransALL$Population, maleTransALL$Individual), FUN=mean)
femaleTransTraitsAvg <- aggregate(femaleTransALL, by=list(femaleTransALL$Population, femaleTransALL$Individual), FUN=mean)

maleTransTraitsAvgSort <- maleTransTraitsAvg[order(maleTransTraitsAvg$Group.1),,]
femaleTransTraitsAvgSort <- femaleTransTraitsAvg[order(femaleTransTraitsAvg$Group.1),,] 

##final data for use in pca

ordinationMaleTransTraitswodevo <- cbind(maleTransTraitsAvgSort[3:6], maleTransTraitsAvgSort$Group.1)
ordinationFemaleTransTraits <- cbind(femaleTransTraitsAvgSort[3:6], femaleTransTraitsAvgSort$Group.1)
colnames(ordinationMaleTransTraitswodevo) <- c("DtG", "Dt10G", "PS", "LL", "Population")
colnames(ordinationFemaleTransTraits) <- c("DtG", "Dt10G", "PS", "LL", "Population")

ordinationMaleTransTraitswdevo <- cbind(maleTransTraitsAvgSort[3:8], maleTransTraitsAvgSort$Group.1)
colnames(ordinationMaleTransTraitswdevo) <- c("DtG", "Dt10G", "PS", "LL","DtP","NP","Population")


##### code for Gmax 

## males

## global
pcaM <- prcomp(ordinationMaleTransTraitswdevo[1:6],center=T,scale=TRUE)
summary(pcaM)
pcaM
traitloadM <- as.data.frame(pcaM$rotation)
evM <- (pcaM$sdev)^2

## by pop
maleA <- subset(ordinationMaleTransTraitswdevo, Population %in% "A")
pcaMA <- prcomp(maleA[1:6],center=T,scale=TRUE)
summary(pcaMA)
pcaMA
(pcaMA$sdev)^2
traitloadMA <- as.data.frame(pcaMA$rotation)
evMA <- (pcaMA$sdev)^2

maleD <- subset(ordinationMaleTransTraitswdevo , Population %in% "D")
pcaMD <- prcomp(maleD[1:6],center=T,scale=TRUE)
summary(pcaMD)
pcaMD
(pcaMD$sdev)^2
traitloadMD <- as.data.frame(pcaMD$rotation)
evMD <- (pcaMD$sdev)^2

maleE <- subset(ordinationMaleTransTraitswdevo , Population %in% "E")
pcaME <- prcomp(maleE[1:6],center=T,scale=TRUE)
summary(pcaME)
pcaME
(pcaME$sdev)^2
traitloadME <- as.data.frame(pcaME$rotation)
evME <- (pcaME$sdev)^2

maleN <- subset(ordinationMaleTransTraitswdevo , Population %in% "N")
pcaMN <- prcomp(maleN[1:6],center=T,scale=TRUE)
summary(pcaMN)
pcaMN
(pcaMN$sdev)^2
traitloadMN <- as.data.frame(pcaMN$rotation)
evMN <- (pcaMN$sdev)^2

maleR <- subset(ordinationMaleTransTraitswdevo , Population %in% "R")
pcaMR <- prcomp(maleR[1:6],center=T,scale=TRUE)
summary(pcaMR)
pcaMR
(pcaMR$sdev)^2
traitloadMR <- as.data.frame(pcaMR$rotation)
evMR <- (pcaMR$sdev)^2

maleS <- subset(ordinationMaleTransTraitswdevo , Population %in% "S")
pcaMS <- prcomp(maleS[1:6],center=T,scale=TRUE)
summary(pcaMS)
pcaMS
(pcaMS$sdev)^2
traitloadMS <- as.data.frame(pcaMS$rotation)
evMS <- (pcaMS$sdev)^2

maleW <- subset(ordinationMaleTransTraitswdevo , Population %in% "W")
pcaMW <- prcomp(maleW[1:6],center=T,scale=TRUE)
summary(pcaMW)
pcaMW
(pcaMW$sdev)^2
traitloadMW <- as.data.frame(pcaMW$rotation)
evMW <- (pcaMW$sdev)^2

### females

## global
pcaF <- prcomp(ordinationFemaleTransTraits[1:4],center=T,scale=T)
summary(pcaF)
pcaF
traitloadF <- as.data.frame(pcaF$rotation)
evF <- (pcaF$sdev)^2

## by pop
femaleA <- subset(ordinationFemaleTransTraits , Population %in% "A")
pcaFA <- prcomp(femaleA[1:4],center=T,scale=TRUE)
summary(pcaFA)
pcaFA
(pcaFA$sdev)^2
traitloadFA <- as.data.frame(pcaFA$rotation)
evFA <- (pcaFA$sdev)^2

femaleD <- subset(ordinationFemaleTransTraits , Population %in% "D")
pcaFD <- prcomp(femaleD[1:4],center=T,scale=TRUE)
summary(pcaFD)
pcaFD
(pcaFD$sdev)^2
traitloadFD <- as.data.frame(pcaFD$rotation)
evFD <- (pcaFD$sdev)^2

femaleE <- subset(ordinationFemaleTransTraits , Population %in% "E")
pcaFE <- prcomp(femaleE[1:4],center=T,scale=TRUE)
summary(pcaFE)
pcaFE
(pcaFE$sdev)^2
traitloadFE <- as.data.frame(pcaFE$rotation)
evFE <- (pcaFE$sdev)^2

femaleN <- subset(ordinationFemaleTransTraits , Population %in% "N")
pcaFN <- prcomp(femaleN[1:4],center=T,scale=TRUE)
summary(pcaFN)
pcaFN
(pcaFN$sdev)^2
traitloadFN <- as.data.frame(pcaFN$rotation)
evFN <- (pcaFN$sdev)^2

femaleR <- subset(ordinationFemaleTransTraits , Population %in% "R")
pcaFR <- prcomp(femaleR[1:4],center=T,scale=TRUE)
summary(pcaFR)
pcaFR
(pcaFR$sdev)^2
traitloadFR <- as.data.frame(pcaFR$rotation)
evFR <- (pcaFR$sdev)^2

femaleS <- subset(ordinationFemaleTransTraits , Population %in% "S")
pcaFS <- prcomp(femaleS[1:4],center=T,scale=TRUE)
summary(pcaFS)
pcaFS
(pcaFS$sdev)^2
traitloadFS <- as.data.frame(pcaFS$rotation)
evFS <- (pcaFS$sdev)^2

femaleW <- subset(ordinationFemaleTransTraits , Population %in% "W")
pcaFW <- prcomp(femaleW[1:4],center=T,scale=TRUE)
summary(pcaFW)
pcaFW
(pcaFW$sdev)^2
traitloadFW <- as.data.frame(pcaFW$rotation)
evFW <- (pcaFW$sdev)^2

populationVectorLoadingF <- c("A","","","",
                             "D","","","",
                             "E","","","",
                             "N","","","",
                             "R","","","",
                             "S","","","",
                             "W","","","")

populationVectorLoadingM <- c("A","","","","","",
                             "D","","","","","",
                             "E","","","","","",
                             "N","","","","","",
                             "R","","","","","",
                             "S","","","","","",
                             "W","","","","","")

populationVector <- c("A", "D", "E", "N", "R", "S", "W")

maleTraitLoadings <- as.data.frame(rbind(traitloadMA,traitloadMD,traitloadME,traitloadMN,traitloadMR,
                                         traitloadMS,traitloadMW))
maleTraitLoadingswPop <- cbind(populationVectorLoadingM,maleTraitLoadings)

femaleTraitLoadings <- as.data.frame(rbind(traitloadFA,traitloadFD,traitloadFE,traitloadFN,traitloadFR,
                                           traitloadFS,traitloadFW))
femaleTraitLoadingswPop <- cbind(populationVectorLoadingF,femaleTraitLoadings)

maleEigan <- as.data.frame(rbind(evMA,evMD,evME,evMN,evMR,
                                 evMS,evMW))
maleEiganwPop <- cbind(populationVector,maleEigan)

femaleEigan <- as.data.frame(rbind(evFA,evFD,evFE,evFN,evFR,
                                   evFS,evFW))
femaleEiganwPop <- cbind(populationVector,femaleEigan)


write.csv(traitloadM, "tables/maleTraitLoadingGlobal.csv")
write.csv(traitloadF, "tables/femaleTraitLoadingGlobal.csv")
write.csv(evM, "tables/maleEiganGlobal.csv")
write.csv(evF, "tables/femaleEiganGlobal.csv")

write.csv(maleTraitLoadingswPop, "tables/maleTraitLoadingswPop.csv")
write.csv(femaleTraitLoadingswPop, "tables/femaleTraitLoadingswPop.csv")
write.csv(maleEiganwPop, "tables/maleEiganwPop.csv")
write.csv(femaleEiganwPop, "tables/femaleEiganwPop.csv")


#########################

## code to make plots/figures

png("figures/malewDevo_PCA_2019.png", width = 8, height = 6, units = 'in', res = 300)

autoplot(pcaM, data=ordinationMaleTransTraitswdevo,loadings=TRUE,colour="gray65", size=3,
         frame=TRUE,frame.type="norm",frame.colour="Population",
         frame.level=0.68, frame.alpha=0.25,
         loadings.label=TRUE,loadings.label.size=8, loadings.colour="black", 
         loadings.label.colour="black",loadings.label.vjust=c(-0.5,1.6,-1,-1,0,1.2),
         loadings.label.hjust=c(0.5,0.5,0.5,0,-0.25,0.5),
         cex=3, main="Males") +
        theme(panel.background=element_rect(fill="white"), plot.title=element_text(size=30, hjust=0.5),
              axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
              legend.text=element_text(size=15),legend.title=element_text(size=20),
              axis.text=element_text(size=15)) +
  scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkslategrey","darkblue","purple",
                                "deeppink"))

dev.off()

png("figures/female_PCA_2019.png", width = 8, height = 6, units = 'in', res = 300)

autoplot(pcaF, data=ordinationFemaleTransTraits,loadings=TRUE,colour="gray65", size=3,
         frame=TRUE,frame.type="norm",frame.colour="Population",
         frame.level=0.68, frame.alpha=0.25,
         loadings.label=TRUE,loadings.label.size=8, loadings.colour="black", 
         loadings.label.colour="black",loadings.label.vjust=c(-0.5,1.6,-1,-1),
         loadings.label.hjust=c(0.5,0.5,0.5,0.5),
         cex=3, main="Females") +
  theme(panel.background=element_rect(fill="white"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15)) +
  scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkslategrey","darkblue","purple",
                              "deeppink"))

dev.off()


pcaMwodevo <- prcomp(ordinationMaleTransTraitswodevo[1:4],center=T,scale=TRUE)

png("figures/male_PCA_2019.png", width = 8, height = 6, units = 'in', res = 300)

autoplot(pcaMwodevo, data=ordinationMaleTransTraitswodevo,loadings=TRUE,colour="gray65", size=3,
         frame=TRUE,frame.type="norm",frame.colour="Population",
         frame.level=0.68, frame.alpha=0.25,
         loadings.label=TRUE,loadings.label.size=8, loadings.colour="black", 
         loadings.label.colour="black",loadings.label.vjust=c(-0.5,1.6,-1,-1),
         loadings.label.hjust=c(0.5,0.5,0.5,0),
         cex=3, main="Males") +
  theme(panel.background=element_rect(fill="white"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15)) +
  scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkslategrey","darkblue","purple",
                              "deeppink"))

dev.off()
