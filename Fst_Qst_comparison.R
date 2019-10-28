
## Qst/Fst comparison

## analyses for McDaniel et al done with R version 3.5.3

#### load in data
maledata <- na.omit(read.csv("data/MaleENA.csv"))
maledata$Individual <- as.character(maledata$Individual)
maledata$Individual_unique <- as.character(maledata$Individual_unique)
maledata$Clone <- as.character(maledata$Clone)

femaledata <- na.omit(read.csv("data/FemaleENA.csv"))
femaledata$Individual <- as.character(femaledata$Individual)
femaledata$Individual_unique <- as.character(femaledata$Individual_unique)
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

### this analysis uses a unique individual number (Individual_unique)
malePopulation <- as.character(maledata$Population)
maleIndividual <- as.character(maledata$Individual)
maleIndividual_unique <- as.character(maledata$Individual_unique)
maleTransALL <- as.data.frame(cbind(maledatatransDtG,maledatatransDt10G,maledatatransPS,maledatatransLL,maledatatransDtP,maledatatransNP,malePopulation,maleIndividual, maleIndividual_unique))
maleTransALL$maledatatransDtG <- as.numeric(as.character(maleTransALL$maledatatransDtG))
maleTransALL$maledatatransDt10G <- as.numeric(as.character(maleTransALL$maledatatransDt10G))
maleTransALL$maledatatransPS <- as.numeric(as.character(maleTransALL$maledatatransPS))
maleTransALL$maledatatransLL <- as.numeric(as.character(maleTransALL$maledatatransLL))
maleTransALL$maledatatransDtP <- as.numeric(as.character(maleTransALL$maledatatransDtP))
maleTransALL$maledatatransNP <- as.numeric(as.character(maleTransALL$maledatatransNP))

femalePopulation <- as.character(femaledata$Population)
femaleIndividual <- as.character(femaledata$Individual)
femaleIndividual_unique <- as.character(femaledata$Individual_unique)
femaleTransALL <- as.data.frame(cbind(femaledatatransDtG,femaledatatransDt10G,femaledatatransPS,femaledatatransLL,femalePopulation,femaleIndividual,femaleIndividual_unique))
femaleTransALL$femaledatatransDtG <- as.numeric(as.character(femaleTransALL$femaledatatransDtG))
femaleTransALL$femaledatatransDt10G <- as.numeric(as.character(femaleTransALL$femaledatatransDt10G))
femaleTransALL$femaledatatransPS <- as.numeric(as.character(femaleTransALL$femaledatatransPS))
femaleTransALL$femaledatatransLL <- as.numeric(as.character(femaleTransALL$femaledatatransLL))


colnames(maleTransALL) <- c("DtG","Dt10G","PS","LL","DtP","NP","Population","Individual","Individual_unique")
colnames(femaleTransALL) <- c("DtG","Dt10G","PS","LL","Population","Individual","Individual_unique")

########################################


only.responses <- cbind(maleTransALL$DtG, maleTransALL$Dt10G, 
                        maleTransALL$PS,maleTransALL$LL,maleTransALL$DtP,maleTransALL$NP)
colnames(only.responses) <- c("DtG", "Dt10G", "LL", "PS", "DtP", "NP")

mu0.vec <- apply(only.responses,2,mean)

Sigma0.mat <- var(only.responses)

#### model for sampling under a normal distribution
my.rmvn <- function(n,mu.vec, cov.mat){
  
  p <- length(mu.vec);
  Tau <- chol(cov.mat);
  Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);
  out <- matrix(0,nrow=p,ncol=n);
  for(i in 1:n){
    
    Z <- Zmat[,i];
    out[,i] <- t(Tau)%*%Z + mu.vec
    
  }
  return(out)
}

nsamps <- nrow(maleTransALL)
one.rand.dataHo <- my.rmvn(n=nsamps, mu.vec=mu0.vec, cov.mat=Sigma0.mat)   

##### Parametric bootstrap under the 'Population' model:

# doing the parametric bootstrap

###### male bootstrap
B <- 2000
Qst.bootmatmale <- matrix(0,nrow=B, ncol=6)
nsamps <- nrow(maleTransALL)
ns.vec <- table(maleTransALL$Population)
ns.cum <- cumsum(ns.vec)
npops <- length(ns.vec)
nstart <- c(0,ns.cum[1:(npops-1)])+1 
nend   <- ns.cum


for(i in 1:B){
  
  # Create the (random) data under the null
  
  response.mat <- matrix(0,nrow=nsamps,ncol=6)
  colnames(response.mat) <- colnames(only.responses)
  for(j in 1:npops){
    
    nj <- ns.vec[j]
    st <- nstart[j]
    en <- nend[j]
    jth.dat <- only.responses[st:en,]
    jt.inds <- maleTransALL$Individual_unique[st:en]
    jt.fac <- as.factor(as.numeric(jt.inds))
    levs.jt <- levels(jt.fac)
    nlevs <- length(levs.jt)
    ns.perind <- table(jt.inds)
    pop.responses <- matrix(0,nrow=nj,ncol=4)

    jth.means <- apply(jth.dat,2,mean)
    jth.var <- var(jth.dat)
    jthpopdat <- t(my.rmvn(n=nj,mu.vec=jth.means, cov.mat=jth.var))
    response.mat[st:en,] <- jthpopdat
    
  }
  
  Population <- maleTransALL$Population
  Individual_unique <- maleTransALL$Individual_unique
  boot.datmale <- data.frame(Population,Individual_unique,response.mat)
  
  #manova
  boot.manovamale <- manova(cbind(DtG, Dt10G,
                              LL, PS, DtP, NP) ~
                          Population+Individual_unique, data = boot.datmale)
  
  #pull out sum of squares from manova output
  boot.sum_aovmale <- summary.aov(boot.manovamale)
  
  #extracting each trait's mean sum of squares 
  boot.SSDtGmale <- as.data.frame(boot.sum_aovmale$` Response DtG`$`Sum Sq`)
  boot.SSDt10Gmale <- as.data.frame(boot.sum_aovmale$` Response Dt10G`$`Sum Sq`)
  boot.SSPSmale <- as.data.frame(boot.sum_aovmale$` Response PS`$`Sum Sq`)
  boot.SSLLmale <- as.data.frame(boot.sum_aovmale$` Response LL`$`Sum Sq`)
  boot.SSDtPmale <- as.data.frame(boot.sum_aovmale$` Response DtP`$`Sum Sq`)
  boot.SSNPmale <- as.data.frame(boot.sum_aovmale$` Response NP`$`Sum Sq`)
  
  #haploid Qst calculation
  boot.DtGmale <- boot.SSDtGmale[1,1]/(boot.SSDtGmale[2,1]+boot.SSDtGmale[1,1])
  boot.Dt10Gmale <- boot.SSDt10Gmale[1,1]/(boot.SSDt10Gmale[2,1]+boot.SSDt10Gmale[1,1])
  boot.PSmale <- boot.SSPSmale[1,1]/(boot.SSPSmale[2,1]+boot.SSPSmale[1,1])
  boot.LLmale <-  boot.SSLLmale[1,1]/(boot.SSLLmale[2,1]+boot.SSLLmale[1,1])
  boot.DtPmale <-  boot.SSDtPmale[1,1]/(boot.SSDtPmale[2,1]+boot.SSDtPmale[1,1])
  boot.NPmale <-  boot.SSNPmale[1,1]/(boot.SSNPmale[2,1]+boot.SSNPmale[1,1])
  
  Qst.bootmatmale[i,] <- c(boot.DtGmale, boot.Dt10Gmale, boot.PSmale, boot.LLmale,
                           boot.DtPmale, boot.NPmale)
  
}


DtGmeanmale <- mean(Qst.bootmatmale[,1])
Dt10Gmeanmale <- mean(Qst.bootmatmale[,2])
PSmeanmale <- mean(Qst.bootmatmale[,3])
LLmeanmale <- mean(Qst.bootmatmale[,4])
DtPmeanmale <- mean(Qst.bootmatmale[,5])
NPmeanmale <- mean(Qst.bootmatmale[,6])

DtGsdmale <- sd(Qst.bootmatmale[,1])
Dt10Gsdmale <- sd(Qst.bootmatmale[,2])
PSsdmale <- sd(Qst.bootmatmale[,3])
LLsdmale <- sd(Qst.bootmatmale[,4])
DtPsdmale <- sd(Qst.bootmatmale[,5])
NPsdmale <- sd(Qst.bootmatmale[,6])

### converting Fst's confidence intervals in to standard deviation
Fstsd <- (0.17027-0.12379)*sqrt(7)/1.96
Fstmean <- 0.14711 

xFst <- seq(-4,4,length=2000)*Fstsd+Fstmean
hxFst <- dnorm(xFst,Fstmean,Fstsd)

xDtGmale <- seq(-4,4,length=2000)*DtGsdmale+DtGmeanmale
hxDtGmale <- dnorm(xDtGmale,DtGmeanmale,DtGsdmale)

xDt10Gmale <- seq(-4,4,length=2000)*Dt10Gsdmale+Dt10Gmeanmale
hxDt10Gmale <- dnorm(xDt10Gmale,Dt10Gmeanmale,Dt10Gsdmale)

xPSmale <- seq(-4,4,length=2000)*PSsdmale+PSmeanmale
hxPSmale <- dnorm(xPSmale,PSmeanmale,PSsdmale)

xLLmale <- seq(-4,4,length=2000)*LLsdmale+LLmeanmale
hxLLmale <- dnorm(xLLmale,LLmeanmale,LLsdmale)

xDtPmale <- seq(-4,4,length=2000)*DtPsdmale+DtPmeanmale
hxDtPmale <- dnorm(xDtPmale,DtPmeanmale,DtPsdmale)

xNPmale <- seq(-4,4,length=2000)*NPsdmale+NPmeanmale
hxNPmale <- dnorm(xNPmale,NPmeanmale,NPsdmale)

#################### females

only.responses <- cbind(femaleTransALL$DtG, femaleTransALL$Dt10G, 
                        femaleTransALL$PS,femaleTransALL$LL)
colnames(only.responses) <- c("DtG", "Dt10G", "LL", "PS")

mu0.vec <- apply(only.responses,2,mean)

Sigma0.mat <- var(only.responses)


#### model for sampling under a normal distribution

my.rmvn <- function(n,mu.vec, cov.mat){
  
  p <- length(mu.vec);
  Tau <- chol(cov.mat);
  Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);
  out <- matrix(0,nrow=p,ncol=n);
  for(i in 1:n){
    
    Z <- Zmat[,i];
    out[,i] <- t(Tau)%*%Z + mu.vec
    
  }
  return(out)
}


nsamps <- nrow(femaleTransALL)
one.rand.dataHo <- my.rmvn(n=nsamps, mu.vec=mu0.vec, cov.mat=Sigma0.mat)   

######female bootstrap over Population 

B <- 2000
Qst.bootmatfemale <- matrix(0,nrow=B, ncol=4)
nsamps <- nrow(femaleTransALL)
ns.vec <- table(femaleTransALL$Population)
ns.cum <- cumsum(ns.vec)
npops <- length(ns.vec)
nstart <- c(0,ns.cum[1:(npops-1)])+1 
nend   <- ns.cum


for(i in 1:B){
  
  # Create the (random) data under the null
  
  response.mat <- matrix(0,nrow=nsamps,ncol=4)
  colnames(response.mat) <- colnames(only.responses)
  for(j in 1:npops){
    
    nj <- ns.vec[j]
    st <- nstart[j]
    en <- nend[j]
    jth.dat <- only.responses[st:en,]
    jt.inds <- femaleTransALL$Individual_unique[st:en]
    jt.fac <- as.factor(as.numeric(jt.inds))
    levs.jt <- levels(jt.fac)
    nlevs <- length(levs.jt)
    ns.perind <- table(jt.inds)
    pop.responses <- matrix(0,nrow=nj,ncol=4)
    
    jth.means <- apply(jth.dat,2,mean)
    jth.var <- var(jth.dat)
    jthpopdat <- t(my.rmvn(n=nj,mu.vec=jth.means, cov.mat=jth.var))
    response.mat[st:en,] <- jthpopdat
    
  }
  
  Population <- femaleTransALL$Population
  Individual_unique <- femaleTransALL$Individual_unique
  boot.datfemale <- data.frame(Population,Individual_unique,response.mat)
  
  #manova
  boot.manovafemale <- manova(cbind(DtG, Dt10G,
                              LL, PS) ~
                          Population+Individual_unique, data = boot.datfemale)
  
  #pull out sum of squares from manova output
  boot.sum_aovfemale <- summary.aov(boot.manovafemale)
  
  #extracting each trait's mean sum of squares 
  boot.SSDtGfemale <- as.data.frame(boot.sum_aovfemale$` Response DtG`$`Sum Sq`)
  boot.SSDt10Gfemale <- as.data.frame(boot.sum_aovfemale$` Response Dt10G`$`Sum Sq`)
  boot.SSPSfemale <- as.data.frame(boot.sum_aovfemale$` Response PS`$`Sum Sq`)
  boot.SSLLfemale <- as.data.frame(boot.sum_aovfemale$` Response LL`$`Sum Sq`)
  
  #haploid Qst calculation
  boot.DtGfemale <- boot.SSDtGfemale[1,1]/(boot.SSDtGfemale[2,1]+boot.SSDtGfemale[1,1])
  boot.Dt10Gfemale <- boot.SSDt10Gfemale[1,1]/(boot.SSDt10Gfemale[2,1]+boot.SSDt10Gfemale[1,1])
  boot.PSfemale <- boot.SSPSfemale[1,1]/(boot.SSPSfemale[2,1]+boot.SSPSfemale[1,1])
  boot.LLfemale <-  boot.SSLLfemale[1,1]/(boot.SSLLfemale[2,1]+boot.SSLLfemale[1,1])  
  
  Qst.bootmatfemale[i,] <- c(boot.DtGfemale, boot.Dt10Gfemale, boot.PSfemale, boot.LLfemale)
  
}


DtGmeanfemale <- mean(Qst.bootmatfemale[,1])
Dt10Gmeanfemale <- mean(Qst.bootmatfemale[,2])
PSmeanfemale <- mean(Qst.bootmatfemale[,3])
LLmeanfemale <- mean(Qst.bootmatfemale[,4])

DtGsdfemale <- sd(Qst.bootmatfemale[,1])
Dt10Gsdfemale <- sd(Qst.bootmatfemale[,2])
PSsdfemale <- sd(Qst.bootmatfemale[,3])
LLsdfemale <- sd(Qst.bootmatfemale[,4])

### converting Fst's confidence intervals in to standard deviation
Fstsd <- (0.17027-0.12379)*sqrt(7)/1.96
Fstmean <- 0.14711 


xFst <- seq(-4,4,length=2000)*Fstsd+Fstmean
hxFst <- dnorm(xFst,Fstmean,Fstsd)


xDtGfemale <- seq(-4,4,length=2000)*DtGsdfemale+DtGmeanfemale
hxDtGfemale <- dnorm(xDtGfemale,DtGmeanfemale,DtGsdfemale)

xDt10Gfemale <- seq(-4,4,length=2000)*Dt10Gsdfemale+Dt10Gmeanfemale
hxDt10Gfemale <- dnorm(xDt10Gfemale,Dt10Gmeanfemale,Dt10Gsdfemale)

xPSfemale <- seq(-4,4,length=2000)*PSsdfemale+PSmeanfemale
hxPSfemale <- dnorm(xPSfemale,PSmeanfemale,PSsdfemale)

xLLfemale <- seq(-4,4,length=2000)*LLsdfemale+LLmeanfemale
hxLLfemale <- dnorm(xLLfemale,LLmeanfemale,LLsdfemale)


### plot code

## code to fill area under Fst curve
x <- seq(-3.5,3.5,length=100)*Fstsd + Fstmean
y <- dnorm(x,Fstmean,Fstsd)

## male Qst plot

png("figures/male_Qst_Fst_dist_2019.png", width = 8, height = 6, units = 'in', res = 300)

plot(xFst,hxFst,xlim=c(0,1),ylim=c(0,14),col="black",type="l",lwd=4, ylab=""
     , xlab="Among-population differentiation", main="Males", cex.axis=1.5, cex.lab=1.5,
     cex.main=2)

polygon(c( x[x>=0], -0.05 ),  c(y[x>=0],0 ), col="gray90")

lines(xFst,hxFst,col="black",lwd=5, lty=1)

lines(xPSmale,hxPSmale,col="deepskyblue",lwd=5, lty=2)
lines(xLLmale,hxLLmale,col="midnightblue",lwd=5, lty=6)

lines(xDtPmale,hxLLmale,col="turquoise",lwd=5, lty=1)
lines(xNPmale,hxLLmale,col="cadetblue",lwd=5, lty=5)
lines(xDtGmale,hxDtGmale,col="dodgerblue",lwd=5, lty=3)


legend(x=0.6,y=14,legend=c("Fst (0.147)","PS (0.137)","DtG (0.341)","LL (0.246)","DtP (0.347)","NP (0.076)"),
       col=c("black","deepskyblue","dodgerblue","midnightblue","turquoise","cadetblue"),
       lty=c(1,2,3,6,1,5), lwd=2, x.intersp=0.25, cex=1.5, title="Trait (global value)")

dev.off()

## female Qst plot

png("figures/female_Qst_Fst_dist_2019.png", width = 8, height = 6, units = 'in', res = 300)

plot(xFst,hxFst,xlim=c(0,1),ylim=c(0,14),col="black",type="l",lwd=4, ylab=""
     , xlab="Among-population differentiation", main="Females", cex.axis=1.5, cex.lab=1.5,
     cex.main=2)
polygon(c( x[x>=0], -0.05 ),  c(y[x>=0],0 ), col="gray90")

lines(xFst,hxFst,col="black",lwd=5, lty=1)

lines(xPSfemale,hxPSfemale,col="firebrick1",lwd=5, lty=2)
lines(xDtGfemale,hxDtGfemale,col="salmon",lwd=5, lty=3)
lines(xLLfemale,hxLLfemale,col="darkred",lwd=5, lty=6)

legend(x=0.6,y=14,legend=c("Fst (0.147)","PS (0.04)","DtG (0.08)","LL (0.04)"),
       col=c("black","firebrick1","salmon","darkred"),
       lty=c(1,2,3,6), lwd=2, x.intersp=0.25, cex=1.5, title="Trait (global value)")

dev.off()

