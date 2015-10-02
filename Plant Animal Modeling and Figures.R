# Modeling terrestrial vertebrate species richness and functional diversity as a function of site level vegetation properties and carbon storage
rm(list=ls())

#Load plant covariates from csv file
plot.VGmean <- read.csv(file="PlantDiversityCalculations.csv")
plot.VGvar <- read.csv(file="PlantDiversityVariances.csv")

# Merge response and predictor variables into a single table
CTaverages <- read.csv("CTaverages_overall.csv")
AnimalFD <- read.csv("FunctionalDiversity_Overall_20May2014.csv")
model.data <- merge(CTaverages, AnimalFD, by.x="Site.Code", by.y="Site.Code", all=TRUE)
ShannonDist <- read.csv("ShannonIndex_Distribution.csv")
model.data <- merge(model.data, ShannonDist, by.x="Site.Code", by.y="Site.Code", all=FALSE)
CT.FDisDist <- read.csv("FunctionalDiversity_Overall_Distribution_19May2014.csv")
model.data <- merge(model.data, CT.FDisDist, by.x="Site.Code", by.y="Site.Code", all=FALSE)
model.data <- merge(model.data, plot.VGmean, by.x="Site.Code", by.y="Site.Code", all=TRUE)

# Add annual rainfall from Worldclim extracted at 2.5 arc-minutes resolution
rainWC <- read.csv("SiteLatitudeLongitude.csv")
AnnualRainWC <- rainWC[,-3]
AnnualRainWC <- AnnualRainWC[,-2]
model.data <- merge(model.data, AnnualRainWC, by.x="Site.Code", by.y="Site.Code", all=FALSE)

# Define CV function
CV <- function(data){
  sd(data)/mean(data)
}

elev.range <- function(data){
  max(data) - min(data)
}
# Calculate ELEVATION CV and add elevation to model.data
elevation.data <- read.csv("CT_edgedist_elevation_final.csv") 
elevation.mean <- aggregate(elevation.data$Elevation ~ elevation.data$Site.Code, FUN=mean)
elevation.CV <- aggregate(elevation.data$Elevation ~ elevation.data$Site.Code, FUN=CV)[2]
elevation.range <- aggregate(elevation.data$Elevation ~ elevation.data$Site.Code, FUN=elev.range)[2]
elevation <- cbind(elevation.mean, elevation.CV, elevation.range)
colnames(elevation) <- c("Site.Code", "Elev.Mean", "Elev.CV", "Elev.Range")
model.data <- merge(model.data, elevation, by.x="Site.Code", by.y="Site.Code", all=FALSE)

# ADD LATITUDE for each TEAM site as a covariate
Latitude_orig <- read.csv("Latitude_MeanSiteCT.csv")
Latitude_match <- Latitude_orig
Latitude_match$Latitude <- abs(Latitude_match$Latitude)
model.data <- merge(model.data, Latitude_match, by.x="Site.Code", by.y="Site.Code", all=FALSE)

# ADD FOREST LOSS data from Alex Zvoleff's calculations
ForestLoss <- read.csv("GFC_Forest_Change_Summary.csv")
names(ForestLoss) <- c("Site.Code", "ForestLossSite", "ForestLossZOI", "PA_area", "SA_area", "ZOI_area")
model.data <- merge(model.data, ForestLoss, by.x="Site.Code", by.y="Site.Code", all=FALSE)

# PERFORM TRANSFORMATIONS OF PREDICTOR VARIABLES 
model.data$ForestLossZOI <- log(model.data$ForestLossZOI)
model.data$PA_area <- log(model.data$PA_area)
model.data$CT.FDiv <- rep(0, dim(model.data)[1])

# Format merged continuous data as data frame with numeric values and site codes as row names
ModelData <- model.data
rownames(ModelData) <- ModelData$Site.Code
ModelData <- ModelData[,-1]
MData <- matrix(NA, dim(ModelData)[1], dim(ModelData)[2])
for(i in 1:dim(ModelData)[2])  {
  MData[,i] <- as.numeric(as.character(ModelData[,i]))
  MData
}
MData <- as.data.frame(MData)
colnames(MData) <- colnames(ModelData)

# Scale predictor variables
MData <- cbind(MData[,1:17], scale(MData[,18:dim(MData)[2]]))
rownames(MData) <- model.data$Site.Code
MData <- as.data.frame(MData)

# Add categorical variables for random effects
Year <- c(2011, 2011, 2012, 2012, 2011, 2011, 2011, 2011, 2012, 2011, 2011, 2011, 2011, 2012)
Continent <- c("Asia", "America", "Africa", "America", "America", "Africa", "America", "Africa", "Asia", "Madagascar", "Africa", "America", "America", "America")
Mdata <- cbind(MData, Year, Continent)
Asia <- c(1,0,0,0,0,0,0,0,1,0,0,0,0,0)
Africa <- c(0,0,1,0,0,1,0,1,0,0,1,0,0,0)
Madagascar <- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
Mdata <- cbind(Mdata, Asia, Africa, Madagascar)

extras <- as.data.frame(cbind(ForestLoss$ForestLossZOI, ForestLoss$PA_area, Latitude_orig$Latitude))
extras <- cbind(Latitude_orig$Site.Code, extras)
extras <- extras[-9,]
extras <- extras[-6,]
colnames(extras) <- c("Site.Code2", "ForestLossZOI2", "PA_area2", "Latitude2")
 
################################## END DATA FORMATTING ###########################

# Function to visualize pairwise comparisons
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

###### MODEL Terrestrial Vertebrate SPECIES RICHNESS

library(lme4)
library(MASS) 
library(MuMIn)

fitRich <- lm(CT.median ~  V.Cstorage2 + V.TShan + V.NStemsT + WC_Bio12 + Elev.CV + ForestLossZOI + Latitude + PA_area + Africa + Asia + Madagascar, data=Mdata)
allRich.dredge <- dredge(fitRich, beta=TRUE, evaluate=TRUE, rank="AICc", trace=TRUE, extra=list("confint", "adjR^2"))
Richconfset.95p <- get.models(allRich.dredge, cumsum(weight) <= .95)
allRich <- model.avg(Richconfset.95p, beta=TRUE, fit=TRUE)
Rich95output <- model.sel(Richconfset.95p)
RichAlloutput <- model.sel(allRich.dredge)
summary(allRich)

plot(Mdata$CT.median, resid(fitRich), xlab="Global Model Residuals", ylab="Predicted Species Richness")
hist(resid(fitRich), main="", xlab="Global Model Residuals")

# VISUALIZE  FUNCTIONAL DIVERSITY
set.panel(2,1)
par(mar=c(3,2,2,1))
hist(Mdata$CT.FDis, main="Functional Dispersion")
boxplot(Mdata$CT.FDisMedian~Mdata$Continent, ylim=c(0.25,0.35))
set.panel()

###### MODEL MAMMAL Vertebrate FUNCTIONAL DIVERSITY

fitFD <- lm(CT.FDisMedian ~ V.Cstorage2 + V.TShan + V.NStemsT + WC_Bio12 + Elev.CV + ForestLossZOI + Latitude + PA_area + Africa + Asia + Madagascar, data=Mdata)
allFD.dredge <- dredge(fitFD, beta=TRUE, evaluate=TRUE, rank="AICc", trace=TRUE, extra=list("confint", "adjR^2"))
FDconfset.95p <- get.models(allFD.dredge, cumsum(weight) <= .95)
allFD <- model.avg(FDconfset.95p, beta=TRUE, fit=TRUE)
FD95output <- model.sel(FDconfset.95p)
FDAlloutput <- model.sel(allFD.dredge)
summary(allFD)

plot(resid(fitFD), Mdata$CT.median, xlab="Global Model Residuals", ylab="Predicted Functional Diversity")
hist(resid(fitFD), main="", xlab="Global Model Residuals")

###### VISUALIZE Terrestrial Vertebrate TAXONOMIC DIVERSITY
set.panel(2,1)
par(mar=c(3,2,2,1))
hist(Mdata$CT.Shannon, main="Taxonomic Diversity")
boxplot(Mdata$CT.Shannon~Mdata$Continent, ylim=c(2,4))
set.panel()

###### MODEL Terrestrial Vertebrate TAXONOMIC DIVERSITY

fitShan <- lm(Shannon.Index ~  V.Cstorage2 + V.TShan + V.NStemsT + WC_Bio12 + Elev.CV + ForestLossZOI + Latitude + PA_area + Africa + Asia + Madagascar, data=Mdata)
allShan.dredge <- dredge(fitShan, beta=TRUE, evaluate=TRUE, rank="AICc", trace=TRUE, extra=list("confint", "adjR^2"))
Shanconfset.95p <- get.models(allShan.dredge, cumsum(weight) <= .95)
allShan <- model.avg(Shanconfset.95p, beta=TRUE, fit=TRUE)
Shan95output <- model.sel(Shanconfset.95p)
ShanAlloutput <- model.sel(allShan.dredge)
summary(allShan)



############ Create Plots #########
pdf(file="AnimalDiversity_Carbon2_Plots_ByRegionColor.pdf", height=2.7)
set.panel(1,3)
par(mar=c(4, 5, 0.5, 0.5))
plot(model.data$V.Cstorage2/1000, Mdata$CT.median, las=1, ylab="Species Richness", xlab="", bty="n", xlim=c(100, 250), 
     ylim=c(10,50), cex.lab=1.4, cex.axis=1.2, pch=c(15,16,17,18)[unclass(as.factor(Mdata$Continent))], 
     col=c("blue","green3","red","black")[unclass(as.factor(Mdata$Continent))], cex=1.5)
legend("bottomleft", expression(R^2*"=-0.03, df=12, p=0.47"), cex=1.15, bty="n")
plot(model.data$V.Cstorage2/1000, Mdata$Shannon.Index, las=1, ylab="Taxonomic Diversity", xlab="", bty="n", xlim=c(100, 250), 
     ylim=c(2.2, 3.4), cex.lab=1.4, cex.axis=1.2, pch=c(15,16,17,18)[unclass(as.factor(Mdata$Continent))], 
     col=c("blue","green3","red","black")[unclass(as.factor(Mdata$Continent))], cex=1.5)
legend("bottomleft", expression(R^2*"=-0.02, df=12, p=0.40"), cex=1.15, bty="n")
plot(model.data$V.Cstorage2/1000, Mdata$CT.FDisMedian, las=1, ylab="", xlab="", bty="n", xlim=c(100, 250), 
     ylim=c(0.24,0.32), cex.lab=1.4, cex.axis=1.2, pch=c(15,16,17,18)[unclass(as.factor(Mdata$Continent))], 
     col=c("blue","green3","red","black")[unclass(as.factor(Mdata$Continent))], cex=1.5)
mtext("Trait Diversity", side=2, line=4)
legend("bottomleft", expression(R^2*"=0.05, df=12, p=0.22"), cex=1.15, bty="n")
legend("topright", legend=c("Asia", "Americas", "Africa", "Madagascar"), col=c("red", "green3", "blue","black"), 
       pch=c(17,16,15,18),  box.col="transparent", cex=1.15)
mtext("               Aboveground Carbon Storage (Mg C per ha)", side=1, outer=TRUE, line=-1.5)
dev.off()



# Map of TEAM site locations
library(maps)
LatLon <- read.csv(file="SiteLatitudeLongitude.csv")
LatLon <- LatLon[-9,]
LatLon <- LatLon[-6,]
set.panel()
pdf(file="TEAM_Map_14Sites.pdf")
par(mar=c(0,0,0,0))
map('world', interior=FALSE, xlim=c(-132, 155), ylim=c(-60, 37), col="gray60")
points(LatLon$Longitude, LatLon$Latitude, col="black", pch=c(16), cex=1)
#legend(x=-132, y=37, legend=LatLon$Site.Code, pch=c(1:14), border="transparent", col="green4", bg="white", box.col="transparent", title="TEAM Sites", title.adj=0.12, cex=0.66)
dev.off()



# Relative Variable Importance Plot
Variable <- c("Elevation CV", "Madagascar", "Africa", "Tree Diversity", "Stem Density", "Forest Loss", "PA Size", "Asia", "Latitude", "Rainfall", 
              "Carbon")
RelVar.Rich <- c(0.88, 0.71, 0.07, 0.06, 0.35, 0.05, 0.12, 0.06, 0.08, 0.08, 0.06)
RelVar.Shan <- c(0.6, 0.26, 0.15, 0.58, 0.4, 0.23, 0.17, 0.15, 0.14, 0.07, 0.07)
RelVar.FD <- c(0.1, 0.11, 0.65, 0.36, 0.27, 0.29, 0.08, 0.15, 0.12, 0.1, 0.07)
RelVar <- cbind(RelVar.Rich, RelVar.Shan, RelVar.FD)
rownames(RelVar) <- Variable

barplotdata <- t(RelVar)
colnames(barplotdata) <- Variable
pdf(file="RelativeVariableImportance_BarPlot_WCBio12_DummyVariables_4August2014.pdf", height=5)
par(mar=c(8,5,2,2))
barplot(barplotdata, beside=TRUE, horiz=FALSE, las=2, ylab="Relative Variable Importance", cex.lab=1, cex.axis=1, ylim=c(0,1))
legend("topright", pch=22, pt.cex=1.5, legend=c("Species Richness", "Taxonomic Diversity", "Trait Diversity"), col=c("black"), pt.bg=c("black", "gray", "gray92"), bty="n")
dev.off()

RelRich <- summary(allRich)[[6]]
RelShan <- summary(allShan)[[6]]
RelFD <- summary(allFD)[[6]]

RelRich[charmatch(names(RelRich), names(RelShan))]
RelShan[charmatch(names(RelFD), names(RelShan))]

## Coefficient Plot
# Code modified from https://gist.github.com/dsparks/818976
library(ggplot2)

Rich.coef <- summary(allRich)[[3]]
rownames(summary(allRich)[[3]])[1:12]
rownames(Rich.coef) <- c("(Intercept)", "Elevation CV", "Madagascar", "Stem Density",  "Africa", "PA Size", "Asia", "Tree Diversity", "Carbon", "Forest Loss", "Rainfall", "Latitude")
  
Shan.coef <- summary(allShan)[[3]]
rownames(summary(allShan)[[3]])[1:12]
rownames(Shan.coef) <- c("(Intercept)",  "Elevation CV", "Stem Density", "Tree Diversity", "Asia", "Madagascar", "Africa", "Forest Loss", 
                         "PA Size", "Latitude", "Rainfall", "Carbon")
  
FD.coef <- summary(allFD)[[3]]
rownames(summary(allFD)[[3]])[1:12]
rownames(FD.coef) <- c("(Intercept)", "Africa",  "Stem Density", "Forest Loss", "Tree Diversity", "Latitude",
                       "Asia", "Madagascar", "Elevation CV", "Carbon", "PA Size", "Rainfall")
  
#graphmodels <- list(summary(allRich)[[3]], summary(allShan)[[3]], summary(allFD)[[3]])
graphmodels <- list(Rich.coef[-1,], Shan.coef[-1,], FD.coef[-1,])#remove intercept from figure by removing its row here
names(graphmodels) <- c("Species Richness", "Taxonomic Diversity", "Trait Diversity")


CoefficientPlot <- function(models, modelnames = ""){
  # models must be a list()
  
  CoefficientTables <- graphmodels
  TableRows <- unlist(lapply(CoefficientTables, nrow))
  
  if(modelnames[1] == ""){
    ModelNameLabels <- rep(paste("Model", 1:length(TableRows)), TableRows)
  } else {
    ModelNameLabels <- rep(modelnames, TableRows)
  }
  
  MatrixofModels <- cbind(do.call(rbind, CoefficientTables), ModelNameLabels)
  MatrixofModels <- data.frame(cbind(rownames(MatrixofModels), MatrixofModels))
  colnames(MatrixofModels) <- c("IV", "Estimate", "StandardError", "AdjSE", "LowerCI", "UpperCI", "ModelName")
  #MatrixofModels$IV <- factor(MatrixofModels$IV, levels = MatrixofModels$IV)
  MatrixofModels$IV <- factor(MatrixofModels$IV, levels = c("(Intercept)", "Madagascar", "Africa", "Asia", "Carbon", "Rainfall", "Latitude", "PA Size", "Forest Loss", "Tree Diversity", "Stem Density", "Elevation CV"))
  
  MatrixofModels[, -c(1, 7)] <- apply(MatrixofModels[, -c(1, 7)], 2, function(x){as.numeric(as.character(x))})
  
  MatrixofModels$ModelName <- factor(MatrixofModels$ModelName, levels=c("Species Richness", "Taxonomic Diversity", "Trait Diversity"))
  
  OutputPlot <- qplot(IV, Estimate, ymin = LowerCI,
                      ymax = UpperCI, data = MatrixofModels, geom = "pointrange",
                      ylab = NULL, xlab = NULL)
  OutputPlot <- OutputPlot + geom_hline(yintercept = 0, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
  OutputPlot <- OutputPlot + facet_grid(~ ModelName) + coord_flip() + theme_bw()
  return(OutputPlot)
}

pdf(file="CoefficientPlot_2October2015.pdf")
CoefficientPlot(graphmodels, modelnames=c("Species Richness", "Taxonomic Diversity", "Trait Diversity"))
dev.off()


# Examine relationship between carbon storage and diversity within continents
summary(lm(Mdata$CT.median[Mdata$Continent=="Africa"] ~ Mdata$V.Cstorage2[Mdata$Continent=="Africa"]))
summary(lm(Mdata$CT.median[Mdata$Continent=="America"] ~ Mdata$V.Cstorage2[Mdata$Continent=="America"]))

summary(lm(Mdata$CT.FDisMedian[Mdata$Continent=="Africa"] ~ Mdata$V.Cstorage2[Mdata$Continent=="Africa"]))
summary(lm(Mdata$CT.FDisMedian[Mdata$Continent=="America"] ~ Mdata$V.Cstorage2[Mdata$Continent=="America"]))

summary(lm(Mdata$Shannon.Index[Mdata$Continent=="Africa"] ~ Mdata$V.Cstorage2[Mdata$Continent=="Africa"]))
summary(lm(Mdata$Shannon.Index[Mdata$Continent=="America"] ~ Mdata$V.Cstorage2[Mdata$Continent=="America"]))
