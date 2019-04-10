
### Fligner-Killeen test for homogeneity of variance

## analyses for McDaniel et al done with R version 3.5.3


### load in data 
combined <- na.omit(read.csv("data/CombinedENA.csv"))
combined$Individual <- as.character(combined$Individual)
combined$Clone <- as.character(combined$Clone)

maledata <- na.omit(read.csv("data/MaleENA.csv"))
maledata$Individual <- as.character(maledata$Individual)
maledata$Clone <- as.character(maledata$Clone)

femaledata <- na.omit(read.csv("data/FemaleENA.csv"))
femaledata$Individual <- as.character(femaledata$Individual)
femaledata$Clone <- as.character(femaledata$Clone)

############################################

# males
fkmaleDtG <- fligner.test(maledata$DtG, maledata$Population)
fkmaleDtG
#Fligner-Killeen:med chi-squared = 56.192, df = 6, p-value = 2.661e-10
#######failed. here's the code for the plot 
plotfkmaleDtG <- plot(DtG~Population, data=maledata) 
plotfkmaleDtG

fkmaleDt10G <- fligner.test(maledata$Dt10G, maledata$Population)
fkmaleDt10G
#Fligner-Killeen:med chi-squared = 20.664, df = 6, p-value = 0.002108
#######failed. here's the code for the plot 
plotfkmaleDt10G <- plot(Dt10G~Population, data=maledata) 
plotfkmaleDt10G

fkmalePS <- fligner.test(maledata$PS, maledata$Population)
fkmalePS 
#Fligner-Killeen:med chi-squared = 5.2976, df = 6, p-value = 0.5062

fkmaleLL  <- fligner.test(maledata$LL, maledata$Population)
fkmaleLL
#Fligner-Killeen:med chi-squared = 10.033, df = 6, p-value = 0.1232

fkmaleDtP <- fligner.test(maledata$DtP, maledata$Population)
fkmaleDtP
#Fligner-Killeen:med chi-squared = 57.712, df = 6, p-value = 1.311e-10
#######failed. here's the code for the plot 
plotfkmaleDtP <- plot(DtP~Population, data=maledata) 
plotfkmaleDtP

fkmaleNP <- fligner.test(maledata$NP, maledata$Population)
fkmaleNP
#Fligner-Killeen:med chi-squared = 4.0155, df = 6, p-value = 0.6746


############################################

# females
fkfemaleDtG <- fligner.test(femaledata$DtG, femaledata$Population)
fkfemaleDtG
#Fligner-Killeen:med chi-squared = 11.349, df = 6, p-value = 0.07816

fkfemaleDt10G <- fligner.test(femaledata$Dt10G, femaledata$Population)
fkfemaleDt10G
#Fligner-Killeen:med chi-squared = 4.9972, df = 6, p-value = 0.5442

fkfemalePS <- fligner.test(femaledata$PS, femaledata$Population)
fkfemalePS 
#Fligner-Killeen:med chi-squared = 18.303, df = 6, p-value = 0.005518
#######failed. here's the code for the plot 
plotfkfemalePS <- plot(PS~Population, data=femaledata) 
plotfkfemalePS

fkfemaleLL  <- fligner.test(femaledata$LL, femaledata$Population)
fkfemaleLL
#Fligner-Killeen:med chi-squared = 6.2111, df = 6, p-value = 0.4

############################################

# sexes combined
fkcombinedDtG <- fligner.test(combined$DtG, combined$Population)
fkcombinedDtG
#Fligner-Killeen:med chi-squared = 29.399, df = 6, p-value = 5.111e-05
####failed. 
plotfkcombinedDtG <- plot(DtG~Population, data=combined) 
plotfkcombinedDtG

fkcombinedDt10G <- fligner.test(combined$Dt10G, combined$Population)
fkcombinedDt10G
#Fligner-Killeen:med chi-squared = 11.137, df = 6, p-value = 0.08425

fkcombinedPS <- fligner.test(combined$PS, combined$Population)
fkcombinedPS 
#Fligner-Killeen:med chi-squared = 24.439, df = 6, p-value = 0.0004335
####failed.
plotfkcombinedPS <- plot(PS~Population, data=combined) 
plotfkcombinedPS

fkcombinedLL  <- fligner.test(combined$LL, combined$Population)
fkcombinedLL
#Fligner-Killeen:med chi-squared = 3.8557, df = 6, p-value = 0.6962


##############################################################################

##### checking for homogeneity in variance POST TRANSFORMATION using boxcox

#### transformation 

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


## combined
combinedtransDtG <- ((combined$DtG)^(-2))
combinedtransDt10G <- ((combined$Dt10G)^(-0.25))
combinedtransPS <- (combined$PS)
combinedtransLL <- ((combined$LL)^(0.25))


### test code 

# males
fkmaleDtGtrans <- fligner.test(maledatatransDtG, maledata$Population)
fkmaleDtGtrans
#Fligner-Killeen:med chi-squared = 14.715, df = 6, p-value = 0.02259
##Still technically failed but a heck of a lot closer#########
fkmaleDtGtransplot <- plot(maledatatransDtG~maledata$Population)
fkmaleDtGtransplot

fkmaleDt10Gtrans <- fligner.test(maledatatransDt10G, maledata$Population)
fkmaleDt10Gtrans
#Fligner-Killeen:med chi-squared = 10.797, df = 6, p-value = 0.09487

fkmalePStrans <- fligner.test(maledatatransPS, maledata$Population)
fkmalePStrans 
#Fligner-Killeen:med chi-squared = 5.2976, df = 6, p-value = 0.5062

fkmaleLLtrans  <- fligner.test(maledatatransLL, maledata$Population)
fkmaleLLtrans
#Fligner-Killeen:med chi-squared = 2.4924, df = 6, p-value = 0.8693

fkmaleDtPtrans <- fligner.test(maledatatransDtP, maledata$Population)
fkmaleDtPtrans
#Fligner-Killeen:med chi-squared = 57.712, df = 6, p-value = 1.311e-10 
##doesn't improve AT ALL #######

fkmaleNPtrans <- fligner.test(maledatatransNP, maledata$Population)
fkmaleNPtrans
#Fligner-Killeen:med chi-squared = 6.2886, df = 6, p-value = 0.3917

##########################################
# females
fkfemaleDtGtrans <- fligner.test(femaledatatransDtG, femaledata$Population)
fkfemaleDtGtrans
#Fligner-Killeen:med chi-squared = 13.481, df = 6, p-value = 0.036

fkfemaleDt10Gtrans <- fligner.test(femaledatatransDt10G, femaledata$Population)
fkfemaleDt10Gtrans
#Fligner-Killeen:med chi-squared = 16.147, df = 6, p-value = 0.01299
#this means it goes from pass before transformation to no longer passing #########
plotfkfemaleDt10Gtrans <- plot(femaledatatransDt10G~femaledata$Population) 
plotfkfemaleDt10Gtrans

fkfemalePStrans <- fligner.test(femaledatatransPS, femaledata$Population)
fkfemalePStrans 
#Fligner-Killeen:med chi-squared = 17.628, df = 6, p-value = 0.007233
#######failed. here's the code for the plot######## 
plotfkfemalePStrans <- plot(femaledatatransPS~femaledata$Population) 
plotfkfemalePStrans

fkfemaleLLtrans  <- fligner.test(femaledatatransLL, femaledata$Population)
fkfemaleLLtrans
#Fligner-Killeen:med chi-squared = 9.0402, df = 6, p-value = 0.1713


# sexes combined
fkcombinedDtGtrans <- fligner.test(combinedtransDtG, combined$Population)
fkcombinedDtGtrans
#Fligner-Killeen:med chi-squared = 15.082, df = 6, p-value = 0.01963
####failed.######### 
plotfkcombinedDtGtrans <- plot(combinedtransDtG~Population, data=combined) 
plotfkcombinedDtGtrans

fkcombinedDt10Gtrans <- fligner.test(combinedtransDt10G, combined$Population)
fkcombinedDt10Gtrans
#Fligner-Killeen:med chi-squared = 20.278, df = 6, p-value = 0.00247
####failed.##########
plotfkcombinedDt10Gtrans <- plot(combinedtransDt10G~Population, data=combined) 
plotfkcombinedDt10Gtrans

fkcombinedPStrans <- fligner.test(combinedtransPS, combined$Population)
fkcombinedPStrans 
#Fligner-Killeen:med chi-squared = 24.439, df = 6, p-value = 0.0004335
####failed.##########
plotfkcombinedPStrans <- plot(combinedtransPS~Population, data=combined) 
plotfkcombinedPStrans

fkcombinedLLtrans  <- fligner.test(combinedtransLL, combined$Population)
fkcombinedLLtrans
#Fligner-Killeen:med chi-squared = 3.8276, df = 6, p-value = 0.7

