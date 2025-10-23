#Testing to just keep functions and og data
rm(list=setdiff(ls(), c("catch.pacfin","getTableAfterPCA", "getMetierClusters")))

## Manhattan Distance (R3) change step 3
Step3=WgetMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="manhattan",param2=NULL)

##Run landed weight instead of ex-vessel revenue
plzwork = data %>% dplyr::select(shotID, PACFIN_SPECIES_COMMON_NAME, LANDED_WEIGHT_LBS)

## Add borderlie spp
choosespp = filter(allcatchmt, LANDED_WEIGHT_MTONS > 4000)

## add borderline gear
choosegear = filter(test, LANDED_WEIGHT_MTONS > 12000) 

## Price randomization 
#Species-specific price
df_randomizedP = plzwork %>%
  group_by(PACFIN_SPECIES_COMMON_NAME) %>%
  mutate(randomized_price = sample(PRICE_PER_POUND, n(), replace = FALSE),
    randomized_price_weight = randomized_price * EXVESSEL_REVENUE
  ) %>%
  ungroup()

df_randomizedP = df_randomizedP[c(1,2,7)]

#Price assigned to any species
df_randomizedP = plzwork %>%
  mutate(randomized_price = sample(PRICE_PER_POUND, n(), replace = FALSE),
         randomized_price_weight = randomized_price * EXVESSEL_REVENUE
  )
df_randomizedP = df_randomizedP[c(1,2,7)]

### HAC change Step3 & run HAC function
set.seed(10403)
Step3=HgetMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="hac",param1="euclidean",param2='ward')

### run Aichinson distance w/ robCompositions package 
library(robCompositions)

  # Remove zeros & keep numeric columns
castingnumeric = casting[,2:32]
detlimit = 7.64 * 10^(-6)

zero_rows = rowSums(castingnumeric == 0) == ncol(castingnumeric)
casting = casting[!zero_rows, ]  

BDLresult = imputeBDLs(x = castingnumeric, maxit = 50, eps = 0.1, method = "pls", 
                       dl = rep(detlimit, ncol(castingnumeric)), R = 50, variation = FALSE)

casting2 = as.data.frame(BDLresult$x)


output = clustCoDa(x = casting3, k = 17, method = "clara", bic = TRUE, verbose = TRUE)



### Bootstrapped Jaccard with clusterboot
library(fpc)
bootJacc = clusterboot(data = casting[,2:32], clustermethod = claraCBI, k = 17, B = 100)

bootJacc$bootbrd
bootJacc$bootrecover
rowMeans(bootJacc$bootresult[,1:17]) > 0.5
rowMeans(bootJacc$bootresult[,1:17]) > 0.5 & rowMeans(bootJacc$bootresult[,1:17]) < 0.6
rowMeans(bootJacc$bootresult[,1:17]) > 0.6 & rowMeans(bootJacc$bootresult[,1:17]) < 0.75
rowMeans(bootJacc$bootresult[,1:17]) > 0.75 & rowMeans(bootJacc$bootresult[,1:17]) < 0.85
rowMeans(bootJacc$bootresult[,1:17]) < 0.85


### Time-split assignment and cluster stability across time

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

#clustering
data = filter(data, LANDING_YEAR %in% yearsindata)

plzwork = data %>%
  filter(data, LANDING_YEAR %in% postassignyears) %>% 
  dplyr::select(shotID, PACFIN_SPECIES_COMMON_NAME, PRICE_PER_POUND, EXVESSEL_REVENUE)

plzwork$product = plzwork$PRICE_PER_POUND * plzwork$EXVESSEL_REVENUE

useme_df = matrix(0, nrow = length(unique(plzwork$shotID)), ncol = length(unique(plzwork$PACFIN_SPECIES_COMMON_NAME)),
                  dimnames = list(unique(plzwork$shotID), unique(plzwork$PACFIN_SPECIES_COMMON_NAME)))

for (i in 1:nrow(useme_df)) {
  row = plzwork$shotID[i]
  col = plzwork$PACFIN_SPECIES_COMMON_NAME[i]
  value = plzwork$product[i]
  useme_df[row, col] = value
}

casting = dcast(formula = shotID ~ PACFIN_SPECIES_COMMON_NAME, data = useme_df, 
                fun.aggregate = sum, na.rm=TRUE)

medoids = Step3$medoids
clusterassignments = numeric(nrow(useme_df))

for(i in 1:nrow(useme_df)) {
  distances = apply(medoids, 1, function(medoid) {
    sqrt(sum((useme_df[i,] - medoid)^2)) 
  })
  clusterassignments[i] = which.min(distances)
}


#### Between and within variance
# Within and between variance of clusters and classification
centerOfGravityClassif=numeric()
withinVarClusters=numeric()
sizeClusti=numeric()
centerOfGravityDatLog=colMeans(datLog)
centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
for(k in 1:nbClust){  # Within variance by cluster
  
  clusti=datLog[which(clusters$clustering==k),]
  if(length(which(clusters$clustering==k))==1)  centerOfGravityClusti=clusti
  else centerOfGravityClusti=colMeans(clusti)
  centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
  sizeClusti[k]=length(which(clusters$clustering==k))
  if(length(which(clusters$clustering==k))==1)  withinVarClusters[k]=0
  else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))
  
}
# Between variance
classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
# Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100

return(list(withinVarClusterOnTot, betweenVarClassifOnTot))