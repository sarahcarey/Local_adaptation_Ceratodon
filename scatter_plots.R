
#### scatter plots

## analyses for McDaniel et al done with R version 3.5.3
## ggplot2 version 3.1.0

install.packages("ggplot2")
library("ggplot2")

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

colnames(femaleTransALL) <- c("DtG","Dt10G","PS","LL","Population","Individual")
colnames(maleTransALL) <- c("DtG","Dt10G","PS","LL","DtP","NP","Population","Individual")

### averaging the traits for scatterplots 

maleTransTraitsAvg <- aggregate(maleTransALL, by=list(maleTransALL$Population, maleTransALL$Individual), FUN=mean)
maleTransTraitsAvgSort <- maleTransTraitsAvg[order(maleTransTraitsAvg$Group.1),,]
MaleTransTraitswdevo <- cbind(maleTransTraitsAvgSort[3:8], maleTransTraitsAvgSort$Group.1)
colnames(MaleTransTraitswdevo) <- c("DtG", "Dt10G", "PS", "LL","DtP","NP","Population")

femaleTransTraitsAvg <- aggregate(femaleTransALL, by=list(femaleTransALL$Population, femaleTransALL$Individual), FUN=mean)
femaleTransTraitsAvgSort <- femaleTransTraitsAvg[order(femaleTransTraitsAvg$Group.1),,] 
FemaleTransTraits <- cbind(femaleTransTraitsAvgSort[3:6], femaleTransTraitsAvgSort$Group.1)
colnames(FemaleTransTraits) <- c("DtG", "Dt10G", "PS", "LL", "Population")



## plot code


## Male DtG vs PS

png("figures/male_DtGvsPS_2019.png", width = 8, height = 6, units = 'in', res = 300)

p <- ggplot(MaleTransTraitswdevo, aes(x=DtG, y=log2(PS), color=Population)) + 
  geom_point(cex=5, pch=16, alpha=0.5) + 
  geom_smooth(method = "lm", aes(color=Population), se=FALSE, lwd=2) +
  theme(panel.background=element_rect(fill="gray90"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15)) + xlab("DtG") + ylab("PS") +
scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
                            "deeppink"))+
  labs(title="Males")
p

dev.off()

## DtG vs DtP

png("figures/male_DtGvsDtP_2019.png", width = 8, height = 6, units = 'in', res = 300)

p <- ggplot(MaleTransTraitswdevo, aes(x=DtG, y=DtP, color=Population)) + 
  geom_point(cex=5, pch=16, alpha=0.5) + 
  geom_smooth(method = "lm", aes(color=Population), se=FALSE, lwd=2) +
  theme(panel.background=element_rect(fill="gray90"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15)) + xlab("DtG") + ylab("DtP") +
  scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
                              "deeppink"))+
  labs(title="Males")
p

dev.off()

## NP vs DtP

png("figures/male_DtPvsNP_2019.png", width = 8, height = 6, units = 'in', res = 300)

p <- ggplot(MaleTransTraitswdevo, aes(x=DtP, y=NP, color=Population)) + 
  geom_point(cex=5, pch=16, alpha=0.5) + 
  geom_smooth(method = "lm", aes(color=Population), se=FALSE, lwd=2) +
  theme(panel.background=element_rect(fill="gray90"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15)) + xlab("DtP") + ylab("NP") +
  scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
                              "deeppink"))+
  labs(title="Males")
p

dev.off()

png("figures/male_DtGvsNP_2019.png", width = 8, height = 6, units = 'in', res = 300)

p <- ggplot(MaleTransTraitswdevo, aes(x=DtG, y=NP, color=Population)) + 
  geom_point(cex=5, pch=16, alpha=0.5) + 
  geom_smooth(method = "lm", aes(color=Population), se=FALSE, lwd=2) +
  theme(panel.background=element_rect(fill="gray90"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15)) + xlab("DtG") + ylab("NP") +
  scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
                              "deeppink"))+
  labs(title="Males")
p

dev.off()


## Female DtG vs PS

png("figures/female_DtGvsPS_2019.png", width = 8, height = 6, units = 'in', res = 300)

p <- ggplot(FemaleTransTraits, aes(x=DtG, y=log2(PS), color=Population)) + 
  geom_point(cex=5, pch=16, alpha=0.5) + 
  geom_smooth(method = "lm", aes(color=Population), se=FALSE, lwd=2) +
  theme(panel.background=element_rect(fill="gray90"), plot.title=element_text(size=30, hjust=0.5),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20),
        axis.text=element_text(size=15)) + xlab("DtG") + ylab("PS") +
  scale_color_manual(values=c("red", "darkgoldenrod", "forestgreen","darkcyan","darkblue","purple",
                              "deeppink")) +
labs(title="Females")

p

dev.off()

