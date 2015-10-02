# Process Trait Data and Calculate Functional Diversity Metrics
# Combine both mammals and birds into single site level calculation
library(FD)
source("FunctionalTrait_HelperFunctions.R")

# Connect TEAM master species list to the PanTHERIA mammalian trait data

masterlist<-read.csv("master_species_list_updated_7April2014.csv",h=T) #master list
pantheria <- read.csv("Pantheria_Data_WR05_Aug2008.csv") # PanTHERIA data

# Manually update names in PanTheria data to match IUCN taxonomy used in TEAM master list
pantheria$Binomial <- as.character(pantheria$Binomial)
Binomial2 <- ifelse(pantheria$Binomial=="Profelis aurata", "Caracal aurata",
                    ifelse(pantheria$Binomial=="Cephalophus natalensis", "Cephalophus harveyi",
                           ifelse(pantheria$Binomial=="Eulemur rufus", "Eulemur rufifrons",
                                  ifelse(pantheria$Binomial=="Mazama gouazoubira", "Mazama nemorivaga",
                                         ifelse(pantheria$Binomial=="Catopuma temminckii", "Pardofelis temminckii",
                                                ifelse(pantheria$Binomial=="Manis tricuspis", "Phataginus tricuspis",
                                                       ifelse(pantheria$Binomial=="Manis gigantea", "Smutsia gigantea",
                                                              ifelse(pantheria$Binomial=="Muntiacus muntjak", "Muntiacus montanus", pantheria$Binomial))))))))
pantheria$Binomial <- factor(Binomial2)

#subset TEAM species list from overall PanTHERIA database
matchedlist <- pantheria[match(masterlist$Unique_Name, pantheria$Binomial),]
splist <- cbind(masterlist, matchedlist)

mammals <- splist[splist$Class=="MAMMALIA",]
birds <- splist[splist$Class=="AVES",]
birds <- birds[,1:6]

# Determine the number of species for which data are available for each trait
Nspecies <- vector(mode = "numeric", length=dim(mammals)[2])
for(i in 1:length(Nspecies)){
  Nspecies[i] <- sum(table(mammals[,i]))
}

Ntraits <- cbind(colnames(mammals), Nspecies)

# Subset trait data for which >150 mammal species on TEAM list have data

mammalianTraits <- cbind(mammals[1:4], mammals$AdultHeadBodyLen_mm, mammals$LitterSize, mammals$GR_Area_km2, as.factor(mammals$ActivityCycle), as.factor(mammals$HabitatBreadth), as.factor(mammals$DietBreadth), mammals$Guild)
names(mammalianTraits) <- c("Bin", "Mass", "Class", "Family", "BodyLength", "LitterSize", "GR_Area", "ActivityCycle",  "HabitatBreadth", "DietBreadth", "Guild")

# Use "mammalianTraits" as input data for the f.family.avg function, which requires loading the helper functions from FunctionalTrait_HelperFunctions.R

fam_avg <- f.family.avg(mammalianTraits)
fam_avg <- as.data.frame(fam_avg)

# Apply to bird data
fam_avg_bird <- f.family.avg.bird(birds)
fam_avg_bird <- as.data.frame(fam_avg_bird)

# Fill in missing values using family level averages (mode for factors; median for continuous variables)
# Replace missing values in mammalianTraits with Family-level averages (object "fam_avg") so that they are "corrected" (objects end in "_c")


Mass_c <- ifelse(is.na(mammalianTraits$Mass)==TRUE, fam_avg$Mass_m[match(mammalianTraits$Family, rownames(fam_avg))],mammalianTraits$Mass)
LitterSize_c <- ifelse(is.na(mammalianTraits$LitterSize)==TRUE, fam_avg$LitterSize_m[match(mammalianTraits$Family, rownames(fam_avg))],mammalianTraits$LitterSize)
GR_Area_c <- ifelse(is.na(mammalianTraits$GR_Area)==TRUE, fam_avg$GR_Area_m[match(mammalianTraits$Family, rownames(fam_avg))],mammalianTraits$GR_Area)
ActivityCycle_c <- ifelse(is.na(mammalianTraits$ActivityCycle)==TRUE, fam_avg$ActivityCycle_m[match(mammalianTraits$Family, rownames(fam_avg))],mammalianTraits$ActivityCycle)
Guild_c <- ifelse(is.na(mammalianTraits$Guild)==TRUE, fam_avg$Guild_m[match(mammalianTraits$Family, rownames(fam_avg))],mammalianTraits$Guild)

# Apply to bird data
Mass_c_bird <- ifelse(is.na(birds$Mass)==TRUE, fam_avg_bird$Mass_m[match(birds$Family, rownames(fam_avg_bird))],birds$Mass)
Guild_c_bird <- ifelse(is.na(birds$Guild)==TRUE, fam_avg_bird$Guild_m[match(birds$Family, rownames(fam_avg_bird))],birds$Guild)


# Create data frame out of corrected Trait values to use in functional diversity metrics; coerce factors into factors

Mtraits <- cbind(Mass_c, LitterSize_c, GR_Area_c, ActivityCycle_c, Guild_c)
rownames(Mtraits) <- mammalianTraits$Bin
Mtraits <- as.data.frame(Mtraits)
#Mtraits$HabitatBreadth_c <- as.factor(Mtraits$HabitatBreadth_c)
#Mtraits$DietBreadth_c <- as.factor(Mtraits$DietBreadth_c)
Mtraits$ActivityCycle_c <- as.factor(Mtraits$ActivityCycle_c)
Mtraits$Guild_c <- as.factor(Mtraits$Guild_c)
str(Mtraits)

# Apply to bird data
Btraits <- cbind(Mass_c_bird, Guild_c_bird, birds$Include)
Btraits <- as.data.frame(Btraits)
Btraits$Guild_c_bird <- as.factor(Btraits$Guild_c_bird)
rownames(Btraits) <- birds$Unique_Name
colnames(Btraits) <- c("Mass_c_bird", "Guild_c_bird", "Include")
Btraits <- Btraits[Btraits$Include==1,]
Btraits <- Btraits[,1:2]
BTraits <- cbind(Btraits[,1], rep(NA, dim(Btraits)[1]), rep(NA, dim(Btraits)[1]), rep(NA, dim(Btraits)[1]), Btraits[,2])
colnames(BTraits) <- colnames(Mtraits)
rownames(BTraits) <- rownames(Btraits)

# Combine mammals and birds into single object for functional diversity analyses
Mtraits <- rbind(Mtraits, BTraits)

# Limit list to species in the WPI include list
Includelist <- masterlist[masterlist$Include==1,]
IncludeTraits <- Mtraits[match(Includelist$Unique_Name, rownames(Mtraits)),]


########################## CALCULATE FD FOR ALL SITES #####################################
# Extract relevant data and create a list to loop over
# Ensure that all relevant data objects are included in the Data.Use object to loop over (and are in the desired output order - i.e. alphabetical)
CAX.Use <- allevents[allevents$Sampling.Period=="2012.01" & allevents$Site.Code=="CAX",]
PSH.Use <- allevents[allevents$Sampling.Period=="2012.01" & allevents$Site.Code=="PSH",]
YAS.Use <- allevents[allevents$Sampling.Period=="2012.01" & allevents$Site.Code=="YAS",]
BBS.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="BBS",]
BCI.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="BCI",]
BIF.Use <- allevents[allevents$Sampling.Period=="2012.01" & allevents$Site.Code=="BIF",]
COU.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="COU",]
KRP.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="KRP",]
MAS.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="MAS",]
NNN.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="NNN",]
RNF.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="RNF",]
UDZ.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="UDZ",]
VB.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="VB-",]
YAN.Use <- allevents[allevents$Sampling.Period=="2011.01" & allevents$Site.Code=="YAN",]

Data.Use <- list(BBS.Use, BCI.Use, BIF.Use, CAX.Use, COU.Use, KRP.Use, MAS.Use, NNN.Use, PSH.Use, RNF.Use, UDZ.Use, VB.Use, YAN.Use, YAS.Use)
library(FD)
CTSite_FD <- list()
B.CTSite_FD <- list()
CT.Shannon <- vector()
B.CT.Shannon <- vector()
splist <- read.csv("master_species_list_updated_7April2014.csv",h=T) #master list
sitelist <- list()
Msplist <- list()
Bsplist <- list()
BirdTraits <- list()
SiteTraits <- list()

wpi_weights <- read.csv(file="wpi_weights.csv")
wpi_use <- list()
psi_weights <- list()
psi_use <- list()
CTSite_FD <- list()


# Loop over sites
for(i in 1:length(Data.Use)){

# Extract trait data
  sitelist[[i]] <- unique(Data.Use[[i]]$bin)
  sitelist
  sitelist[[i]] <- factor(sitelist[[i]])
  
  SiteTraits[[i]] <- subset(IncludeTraits, rownames(IncludeTraits) %in%  sitelist[[i]])
  
# Extract weights and format mammal data
  #SiteTraits[[i]] <- Mtraits[match(Msplist[[i]]$Unique_Name, rownames(Mtraits)),]
  SiteTraits[[i]] <- cbind(SiteTraits[[i]][,1:3], droplevels(SiteTraits[[i]][,4:5]))
  wpi_use[[i]] <- wpi_weights[wpi_weights$Site.Code==Data.Use[[i]]$Site.Code[1], ]
  psi_weights[[i]] <- wpi_use[[i]]$psi.median[match(rownames(SiteTraits[[i]]), wpi_use[[i]]$bin)]
  names(psi_weights[[i]]) <- rownames(SiteTraits[[i]])
  psi_use[[i]] <- na.omit(psi_weights[[i]])
  SiteTraits[[i]] <- SiteTraits[[i]][match(names(psi_use[[i]]), rownames(SiteTraits[[i]])),]
  SiteTraits[[i]] <- cbind(SiteTraits[[i]][,1:3], droplevels(SiteTraits[[i]][,4:5]), psi_use[[i]])

# Calculate mammal FD and Shannon diversity indices  
  CTSite_FD[[i]] <- dbFD(SiteTraits[[i]][,1:dim(SiteTraits[[i]])[2]-1], a=psi_use[[i]][1:dim(SiteTraits[[i]])[1]], corr="cailliez", calc.FRic=FALSE)
  CT.Shannon[i] <- diversity(psi_use[[i]], index="shannon")
}  


#Extract FD calculations and create output table 
CT.nbsp <- vector()
CT.sing.sp <- vector()
CT.FEve <- vector()
CT.FDiv <- vector()
CT.FDis <- vector()
CT.RaoQ <- vector()
CWM.Mass <- vector()
CWM.Litter <- vector()
CWM.GR <- vector()
CWM.Activity <- vector()
CWM.Guild <- vector()

for(i in 1:length(CTSite_FD)){
  CT.nbsp[i] <- CTSite_FD[[i]]$nbsp
  CT.sing.sp[i] <- CTSite_FD[[i]]$sing.sp
  CT.FEve[i] <- CTSite_FD[[i]]$FEve
  CT.FDiv[i] <- CTSite_FD[[i]]$FDiv
  CT.FDis[i] <- CTSite_FD[[i]]$FDis
  CT.RaoQ[i] <- CTSite_FD[[i]]$RaoQ
  CWM.Mass[i] <- CTSite_FD[[i]]$CWM$Mass_c
  CWM.Litter[i] <- CTSite_FD[[i]]$CWM$LitterSize_c
  CWM.GR[i] <- CTSite_FD[[i]]$CWM$GR_Area_c
  CWM.Activity[i] <- CTSite_FD[[i]]$CWM$ActivityCycle_c
  CWM.Guild[i] <- CTSite_FD[[i]]$CWM$Guild_c
}

CTweighted <- cbind(CT.nbsp, CT.sing.sp, CT.FEve, CT.FDiv, CT.FDis, CT.RaoQ, CWM.Mass, CWM.Litter, CWM.GR, CT.Shannon)
CTweighted <- as.data.frame(CTweighted)
rownames(CTweighted) <- c("BBS", "BCI", "BIF", "CAX", "COU", "KRP", "MAS", "NNN", "PSH", "RNF", "UDZ", "VB", "YAN", "YAS") 


CTweighted
write.csv(CTweighted, file="FunctionalDiversity_Overall_20May2014.csv", row.names=TRUE)
# Manually add "Site.Code" as first column name for  output file
# Manuall add "-" to "VB-" Site.Code so that later merges will operate correctly