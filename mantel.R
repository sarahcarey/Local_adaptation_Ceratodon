
## Mantel test

## analyses for McDaniel et al done with R version 3.5.3
## ade4 version 1.7-13

install.packages("ade4")
library(ade4)

fst <- read.csv("data/fst.csv",header=F)
fst_female <- read.csv("data/fst_female.csv",header=F)
fst_male <- read.csv("data/fst_male.csv",header=F)

dist <- read.csv("data/dist.csv",header=F)

## all populations 
genALL <- as.dist(fst)
geoALL <- as.dist(dist)
r1 <- mantel.rtest(geoALL,genALL,nrepet=999)
r1
plot(r1)


##just the bottom 5 populations
gen <- as.dist(fst[1:5,1:5])
geo <- as.dist(dist[1:5,1:5])
r1 <- mantel.rtest(geo,gen,nrepet=999)
r1
plot(r1)


## using female sex-linked loci on all populations
genALL <- as.dist(fst_female)
geoALL <- as.dist(dist)
r1 <- mantel.rtest(geoALL,genALL,nrepet=999)
r1
plot(r1)

## using male sex-linked loci on all populations
genALL <- as.dist(fst_male)
geoALL <- as.dist(dist)
r1 <- mantel.rtest(geoALL,genALL,nrepet=999)
r1
plot(r1)


