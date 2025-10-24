## Will Patrone
## University of Washington School of Aquatic and Fishery Sciences
## Reviewer Comments for Metier Analysis for West Coast Groundfish

##Run this to clear your environment but keep the functions and original data
rm(list=setdiff(ls(), c("catch.pacfin","getTableAfterPCA", "WgetMetierClusters")))

##Manhattan Distance - replace line 118 in WCG Metier Analysis.R
Step3=WgetMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="manhattan",param2=NULL)

##Run landed weight instead of ex-vessel revenue - replace lines 87-88 in WCG Metier Analysis.R
EVRbyShot = data %>% 
  dplyr::select(shotID, PACFIN_SPECIES_COMMON_NAME, LANDED_WEIGHT_LBS)

##Add border-lie spp - replace line 57 in WCG Metier Analysis.R
choosespp = filter(allcatchmt, LANDED_WEIGHT_MTONS > 4000)

##Add borderline gear - replace line 71 in WCG Metier Analysis.R
choosegear = filter(test, LANDED_WEIGHT_MTONS > 12000) 

## Price randomization 
#Species-specific price - place in line 84 in WCG Metier Analysis.R
randomizedPwithin = EVRbyShot %>%
  group_by(PACFIN_SPECIES_COMMON_NAME) %>%
  mutate(randomizedprice = sample(PRICE_PER_POUND, n(), replace = FALSE),
         randomizedEVR = randomizedprice * LANDED_WEIGHT_LBS
  ) %>%
  ungroup()

randomizedPwithin = randomizedPwithin[c(1,2,7)]

#Price assigned to any species - place in line 84 in WCG Metier Analysis.R
randomizedPglobal = EVRbyShot %>%
  mutate(randomizedprice = sample(PRICE_PER_POUND, n(), replace = FALSE),
         randomizedEVR = randomizedprice * LANDED_WEIGHT_LBS
  )
randomizedPglobal = randomizedPglobal[c(1,2,7)]

## Run HAC algorithm - replace line 118 in WCG Metier Analysis.R
Step3=HgetMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="hac",param1="euclidean",param2='ward')

## run Aichinson distance w/ robCompositions package - paste in line 98 in WCG Metier Analysis.R 
library(robCompositions)

#Remove zeros & keep numeric columns
castingnumeric = casting[,2:32]
detlimit = 7.64 * 10^(-6) #minimum value in dataset

zerorows = rowSums(castingnumeric == 0) == ncol(castingnumeric)
casting = casting[!zerorows, ]  

BDLresult = imputeBDLs(x = castingnumeric, maxit = 50, eps = 0.1, method = "pls", 
                       dl = rep(detlimit, ncol(castingnumeric)), R = 50, variation = FALSE)

casting2 = as.data.frame(BDLresult$x)

output = clustCoDa(x = casting2, k = 17, method = "clara", bic = TRUE, verbose = TRUE)


### Bootstrapped Jaccard with clusterboot - place in line 98 in WCG Metier Analysis.R
library(fpc)
bootJacc = clusterboot(data = casting[,2:32], clustermethod = claraCBI, k = 17, B = 100)

bootJacc$bootbrd
bootJacc$bootrecover
rowMeans(bootJacc$bootresult[,1:17]) > 0.5
rowMeans(bootJacc$bootresult[,1:17]) > 0.5 & rowMeans(bootJacc$bootresult[,1:17]) < 0.6
rowMeans(bootJacc$bootresult[,1:17]) > 0.6 & rowMeans(bootJacc$bootresult[,1:17]) < 0.75
rowMeans(bootJacc$bootresult[,1:17]) > 0.75 & rowMeans(bootJacc$bootresult[,1:17]) < 0.85
rowMeans(bootJacc$bootresult[,1:17]) < 0.85

##Time-split assignment and cluster stability across time - replace in line 76 in WCG Metier Analysis.R

# 2011-2022 in, 2023 clustered after
yearsindata = 2011:2022  
postassignyears = 2023

# 2011-2021 in, 2022-2023 clustered after
yearsindata = 2011:2021  
postassignyears = 2022:2023

# 2011-2022 in, 2023 clustered after
yearsindata = 2011:2022  
postassignyears = 2023

# 2012-2023 in, 2011 clustered after
yearsindata = 2012:2023 
postassignyears = 2011

# 2013-2023 in, 2011-2012 clustered after
yearsindata = 2013:2023 
postassignyears = 2011:2012

#Clustering
data = filter(data, LANDING_YEAR %in% yearsindata)

EVRbyShot = data %>%
  filter(data, LANDING_YEAR %in% postassignyears) %>% 
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

medoids = Step3$medoids
clusterassignments = numeric(nrow(clusteringDF))

for(i in 1:nrow(clusteringDF)) {
  distances = apply(medoids, 1, function(medoid) {
    sqrt(sum((clusteringDF[i,] - medoid)^2)) 
  })
  clusterassignments[i] = which.min(distances)
}


## Between and within variance
#Within and between variance of clusters and classification - add to end of WgetMetierClusters.R
centerOfGravityClassif=numeric()
withinVarClusters=numeric()
sizeClusti=numeric()
centerOfGravityDatLog=colMeans(datLog)
centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
for(k in 1:nbClust){  #Within variance by cluster
  clusti=datLog[which(clusters$clustering==k),]
  if(length(which(clusters$clustering==k))==1)  centerOfGravityClusti=clusti
  else centerOfGravityClusti=colMeans(clusti)
  centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
  sizeClusti[k]=length(which(clusters$clustering==k))
  if(length(which(clusters$clustering==k))==1)  withinVarClusters[k]=0
  else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))
}
#Between variance
classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
#Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100

return(list(withinVarClusterOnTot, betweenVarClassifOnTot))