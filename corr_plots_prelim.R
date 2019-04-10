
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

colnames(maleTransALL) <- c("DtG","Dt10G","PS","LL","DtP","NP","Population","Individual")
colnames(femaleTransALL) <- c("DtG","Dt10G","PS","LL","Population","Individual")

### averaging the traits for scatterplots 

maleTransTraitsAvg <- aggregate(maleTransALL, by=list(maleTransALL$Population, maleTransALL$Individual), FUN=mean)
femaleTransTraitsAvg <- aggregate(femaleTransALL, by=list(femaleTransALL$Population, femaleTransALL$Individual), FUN=mean)

maleTransTraitsAvgSort <- maleTransTraitsAvg[order(maleTransTraitsAvg$Group.1),,]
femaleTransTraitsAvgSort <- femaleTransTraitsAvg[order(femaleTransTraitsAvg$Group.1),,]


## plot code


## Male DtG vs PS for males

## PS vs DtG
p <- ggplot(maleTransTraitsAvg, aes(x =maledatatransDtG, y =log2(maledatatransPS), color=Group.1)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Group.1), se=FALSE)


## DtG vs DtP


p <- ggplot(maleTransTraitsAvg, aes(x =maledatatransDtG, y =maledatatransDtP, color=Group.1)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Group.1), se=FALSE)

## NP vs DtP


p <- ggplot(maleTransTraitsAvg, aes(x =maledatatransNP, y =maledatatransDtP, color=Group.1)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Group.1), se=FALSE)


plot(femaledata$Dt10G,femaledata$DtG)
abline(lm(femaledata$DtG~femaledata$Dt10G))


plot(maleTransALLscaled$maledatatransLL,maleTransALLscaled$maledatatransPS)
abline(lm(maleTransALLscaled$maledatatransPS~maleTransALLscaled$maledatatransLL))


plot(maleTransALLscaled$maledatatransDtG,maleTransALLscaled$maledatatransDtP)
abline(lm(maleTransALLscaled$maledatatransDtP~maleTransALLscaled$maledatatransDtG))



plot(maledata$DtG,maledata$DtP)
abline(lm(maleTransALLscaled$maledatatransDtP~maleTransALLscaled$maledatatransDtG))




p <- ggplot(maleTransALLscaled, aes(x =maledatatransDtP, y =maledatatransDtG)) + geom_point(cex=5, col="black", pch=21, fill="black", alpha=0.25)
p


#DtG vs DtP
p <- ggplot(maledata, aes(x =DtP, y =DtG)) + geom_point(cex=5, col="black", pch=21, fill="black", alpha=0.25)
p + stat_smooth(method = "lm", size = 1, col="red", aes(x =DtP, y= DtG), data=maledata)

p <- ggplot(maledata, aes(x =DtP, y =DtG, color=Population)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Population), se=FALSE)




#DtP vs NP
p <- ggplot(maledata, aes(x =DtP, y =NP)) + geom_point(cex=5, col="black", pch=21, fill="black", alpha=0.25)
p + stat_smooth(method = "lm", size = 1, col="red", aes(x =DtP, y= NP), data=maledata)

p <- ggplot(maledata, aes(x =DtP, y =NP, color=Population)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Population), se=FALSE)



#PS vs NP
p <- ggplot(maledata, aes(x =log2(PS), y =NP)) + geom_point(cex=5, col="black", pch=21, fill="black", alpha=0.25)
p + stat_smooth(method = "lm", size = 1, col="red", aes(x =log2(maledata$PS), y= maledata$NP))


p <- ggplot(maledata, aes(x =log2(PS), y =NP, color=Population)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Population), se=FALSE)


## PS vs DtP
p <- ggplot(maledata, aes(x =log2(PS), y =DtP)) + geom_point(cex=5, col="black", pch=21, fill="black", alpha=0.25)
p + stat_smooth(method = "lm", size = 1, col="red", aes(x =log2(maledata$PS), y= maledata$DtP))


p <- ggplot(maledata, aes(x =log2(PS), y =DtP, color=Population)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Population), se=FALSE)


p <- ggplot(maleTransTraitsAvg, aes(x =log2(maledatatransPS), y =maledatatransDtP, color=Group.1)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Group.1), se=FALSE)


## PS vs LL
p <- ggplot(maleTransTraitsAvg, aes(x =log2(maledatatransPS), y =log2(maledatatransLL), color=Group.1)) + geom_point(cex=5, pch=16, alpha=0.5)
p + geom_smooth(method = "lm", aes(color=Group.1), se=FALSE)




