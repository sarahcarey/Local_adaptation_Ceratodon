
# beta dispersion to test for heterogeneous multivariate variances 

## analyses for McDaniel et al done with R version 3.5.3
## vegan version 2.5-4

install.packages("vegan")
library("vegan") 

dataforbetatrans <- read.csv("data/dataforbetatrans.csv")
env <- read.csv("data/dispersionEnv.csv", head = TRUE)

## getting dissimilarity indices
beta <- vegdist(dataforbetatrans, method="bray")

betadisp <- betadisper(beta, env$Sex, type="median")

## ANOVA of model
anova(betadisp)

## permutation test
permutest(betadisp)

# plotting ellipses instead of hulls

png("figures/beta_dispersion_2019.png", width = 8, height = 6, units = 'in', res = 300)

par(mar=c(5.1,5.1,4.1,2.1))
plot(betadisp, main = "Beta dispersion of shared traits", ellipse=TRUE, hull=FALSE, 
     col=c("salmon","turquoise"), label = FALSE, lwd=3, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5),
     cex=2, pch=c(17,19), cex.axis=1.5, cex.lab=2, cex.main=2)
    
dev.off()


