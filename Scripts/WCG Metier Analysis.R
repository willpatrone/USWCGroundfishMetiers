## Will Patrone
## University of Washington School of Aquatic and Fishery Sciences
## Metier Analysis for West Coast Groundfish

#Set path for data and random seed
setwd("C:/Users/Will Patrone/Documents/Master's Thesis/Metier Analysis")
set.seed(10403)

#Install Packages
vmstoolsPackages=c("data.table","doBy","lubridate","sf","mixtools","segmented", 
                       "cluster", "maps", "tidyverse", "devtools", "reshape2", "SOAR",
                       "FactoMineR", "writexl", "clustree", "mapdata", "PBSmapping", 
                       "writexl", "readxl", "patchwork",  "tidyverse")
for(i in vmstoolsPackages)       try(install.packages(pkgs=i,repos=getOption("repos")))
install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
install.packages("C:/Users/Will Patrone/Documents/Master's Thesis/Metier Analysis/vmstools_0.76.zip")


#Load Packages
for (pkg in vmstoolsPackages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
library(vmstools)
library(maptools)

#Load data
load("PacFIN.GRND.CompFT.10.Jul.2024.Rdata")


#Refine species list down

data = catch.pacfin[c("LANDING_YEAR", "LANDING_MONTH", "GEAR_CODE", "GEAR_NAME",
                      "PACFIN_GEAR_DESCRIPTION", "PACFIN_CATCH_AREA_DESCRIPTION",
                      "PACFIN_GROUP_CATCH_AREA_CODE", "COUNCIL_CODE", "PORT_CODE",
                      "PACFIN_PORT_NAME", "COUNTY_NAME", "COUNTY_STATE", 
                      "SUBREGION_NAME", "REGION_NAME", 
                      "PACFIN_SPECIES_CODE", "PACFIN_SPECIES_SCIENTIFIC_NAME",
                      "PACFIN_SPECIES_COMMON_NAME", "LANDED_WEIGHT_LBS", 
                      "LANDED_WEIGHT_MTONS", "PRICE_PER_POUND", "EXVESSEL_REVENUE")]

dummy = catch.pacfin[c("PACFIN_SPECIES_SCIENTIFIC_NAME",
          "PACFIN_SPECIES_COMMON_NAME",
          "LANDED_WEIGHT_MTONS")]

dummy = filter(dummy, !str_detect(PACFIN_SPECIES_SCIENTIFIC_NAME, regex("\\bSPP", ignore_case = TRUE)))
dummy = filter(dummy, !str_detect(PACFIN_SPECIES_SCIENTIFIC_NAME, "N/A"))

dummy = dummy[c("PACFIN_SPECIES_COMMON_NAME",
                       "LANDED_WEIGHT_MTONS")]

#Condense dummy to relevant species
allcatchmt= aggregate(.~PACFIN_SPECIES_COMMON_NAME, data = dummy, FUN = sum)
allcatchmt = allcatchmt[order(allcatchmt$LANDED_WEIGHT_MTONS, decreasing = TRUE),]
choosespp = filter(allcatchmt, LANDED_WEIGHT_MTONS > 5000)
sppname = unique(choosespp$PACFIN_SPECIES_COMMON_NAME)

#Refine gear types list down
data$PACFIN_GEAR_DESCRIPTION=ifelse(data$PACFIN_GEAR_DESCRIPTION %in% c("SHRIMP TRAWL, SINGLE RIGGED",
                                                                           "SHRIMP TRAWL, SINGLE OR DOUBLE RIG", 
                                                                           "DANISH/SCOTTISH SEINE (TRAWL)", 
                                                                           "BEAM TRAWL",
                                                                           "PAIR TRAWL"), "OTHER TRAWL GEAR", 
                                       data$PACFIN_GEAR_DESCRIPTION)

dummy = catch.pacfin[c("PACFIN_GEAR_DESCRIPTION", "LANDED_WEIGHT_MTONS")]

test = aggregate(.~PACFIN_GEAR_DESCRIPTION, data = dummy, FUN = sum)
test = test[order(test$LANDED_WEIGHT_MTONS, decreasing = TRUE),]
choosegear = filter(test, LANDED_WEIGHT_MTONS > 25000) 

#Extract species and gears from full data
data = filter(data, PACFIN_GEAR_DESCRIPTION %in% choosegear$PACFIN_GEAR_DESCRIPTION)
data = filter(data, LANDING_YEAR %in% 2011:2023)
data = filter(data, PACFIN_SPECIES_COMMON_NAME %in% choosespp$PACFIN_SPECIES_COMMON_NAME)
data$shotID = paste(data$LANDING_MONTH, data$PACFIN_PORT_NAME, data$PACFIN_GEAR_DESCRIPTION, sep = "_")
data = relocate(data, shotID)

#Calculate EVR and create clustering matrix
EVRbyShot = data %>% 
  dplyr::select(shotID, PACFIN_SPECIES_COMMON_NAME, PRICE_PER_POUND, LANDED_WEIGHT_LBS)

EVRbyShot$product = EVRbyShot$PRICE_PER_POUND * EVRbyShot$LANDED_WEIGHT_LBS

clusteringDF = matrix(0, nrow = length(unique(EVRbyShot$shotID)), ncol = length(unique(EVRbyShot$PACFIN_SPECIES_COMMON_NAME)),
                  dimnames = list(unique(EVRbyShot$shotID), unique(EVRbyShot$PACFIN_SPECIES_COMMON_NAME)))

for (i in 1:nrow(clusteringDF)) {
  row = EVRbyShot$shotID[i]
  col = EVRbyShot$PACFIN_SPECIES_COMMON_NAME[i]
  value = EVRbyShot$product[i]
  clusteringDF[row, col] = value
}

casting = dcast(formula = shotID ~ PACFIN_SPECIES_COMMON_NAME, data = clusteringDF, 
                fun.aggregate = sum, na.rm=TRUE)

#Run CLARA
Step1=extractTableMainSpecies(dat = casting, sppname,
                                 paramTotal=100, paramLogevent = 100)

head(Step1)
rowNamesSave=row.names(Step1)
row.names(Step1)=1:nrow(Step1)
analysisName=paste(data$LANDING_MONTH, data$PACFIN_PORT_NAME, data$COUNTY_STATE, data$PACFIN_GEAR_DESCRIPTION, sep="_")

critPca="NO_PCA"

Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)

Step3=WgetMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="euclidean",param2=NULL)

#Assing cluster number to shots
casting$clusternumber = Step3$clustering
