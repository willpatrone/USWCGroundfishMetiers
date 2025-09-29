WgetMetierClusters = function(datSpecies, datLog, analysisName="", methMetier="clara", param1="euclidean", param2=NULL){
  
  # Load the table linking 3A-CODE (FAO CODE of species) to the species assemblage (level 5).
  data(correspLevel7to5)
  require(lattice)
  require(clustree)  # Ensure clustree is loaded
  
  LE_ID = rownames(datSpecies)
  nbSpec = dim(datSpecies)[2]
  datSpecies = as.matrix(datSpecies, ncol=nbSpec, nrow=length(LE_ID))
  
  print("######## STEP 3 CLUSTERING ########")
  
  t1 = Sys.time()
  print(paste(" --- selected method:", methMetier, " ---"))
  
  # Data frame to store cluster assignments over iterations
  cluster_assignments = data.frame(LE_ID = LE_ID)  # To store IDs of samples
  
  ########################################################################
  # CLARA Clustering Method
  ########################################################################
  
  if (methMetier == "clara") {
    nbLog = nrow(datLog)
    propSample = 0.1
    
    clustersClara.silcoeff = numeric()
    clustersClara.silcoeff[1] = 0
    clustersClara.silcoeff[2] = 0
    clustersClara.silcoeff[3] = 0
    k = 2
    compMax = 1
    
    repeat {
      k = k + 2
      print(k)
      clustersClara = clara(datLog, k, metric = param1, stand = FALSE, samples = 5, sampsize = min(nbLog, round(propSample * nbLog + 10 * k)))
      clustersClara.silcoeff[k] = clustersClara$silinfo$avg.width
      clustersClara = clara(datLog, k + 1, metric = param1, stand = FALSE, samples = 5, sampsize = min(nbLog, round(propSample * nbLog + 10 * (k + 1))))
      clustersClara.silcoeff[k + 1] = clustersClara$silinfo$avg.width
      
      if ((clustersClara.silcoeff[k - 2] < clustersClara.silcoeff[k - 1] & clustersClara.silcoeff[k - 1] > clustersClara.silcoeff[k]) & compMax <= 2) {
        if (compMax == 2) {
          nbClust = k - 1
          print(paste("2nd max =", k - 1))
          print(paste("nbClust =", nbClust))
          break
        } else {
          compMax = compMax + 1
          print(paste("compMax1 =", compMax))
          print(paste("1st max =", k - 1))
        }
      }
      if ((clustersClara.silcoeff[k - 1] < clustersClara.silcoeff[k] & clustersClara.silcoeff[k] > clustersClara.silcoeff[k + 1]) & compMax <= 2) {
        if (compMax == 2) {
          nbClust = k
          print(paste("2nd max =", k))
          print(paste("nbClust =", nbClust))
          break
        } else {
          compMax = compMax + 1
          print(paste("compMax2 =", compMax))
          print(paste("1st max =", k))
        }
      }
      gc(reset = TRUE)
    }
    
    png(paste(analysisName, "Silhouette of the classification.png", sep = "_"), width = 1200, height = 800)
    plot(clustersClara.silcoeff, main = "Silhouette of the classification", xlab = "Number of clusters", ylab = "Silhouette")
    dev.off()
    
    cat("ClaraSilCoeff", clustersClara.silcoeff, "\n")
    
  } else {
    stop("methMetier must be clara")
  }
    
    # CLARA with optimal k
    clusters = clara(datLog, nbClust, metric = param1, stand = FALSE, samples = 5, sampsize = min(nbLog, round(propSample * nbLog + 10 * nbClust)))  # CLARA with optimal k
    summary(clusters)
    
    # Store the cluster assignments in the dataframe for the current k
    cluster_assignments[[paste("cluster", nbClust, sep = "_")]] = clusters$clustering
    
    # Prepare data for clustree (all clustering results)
    clustering_results = data.frame(LE_ID = LE_ID)
    
    # Loop through different values of k (e.g., 2 to nbClust)
    for (k_val in 2:nbClust) {
      clusters_k = clara(datLog, k_val, metric = param1, stand = FALSE, samples = 5, sampsize = min(nbLog, round(propSample * nbLog + 10 * k_val)))
      clustering_results[[as.character(k_val)]] = clusters_k$clustering
    }
    
    # Ensure column names are valid for clustree (they should just be numeric k values)
    colnames(clustering_results) <- c("LE_ID", "K1", "K2","K3", "K4","K5", "K6",
                                                   "K7", "K8","K9", "K10","K11", "K12",
                                                   "K13", "K14","K15", "K16")
    
 
    
    write.csv(clustering_results, "clustresults.csv")
    
    # Generate clustering tree using clustree package
   treee =  clustree(clustering_results, prefix = "K")  # The columns should now be numeric, without the 'cluster_k' prefix.
    
   ggsave("Cluster Tree.jpg", 
          plot = treee, 
          width = 15,  # Width in inches
          height = 10, # Height in inches
          units = "in", # Units can be 'in', 'cm', or 'mm'
          dpi = 300)   # Resolution in DPI
 
  outs = NULL
  outs$withKat = clustering_results
  outs$real = cluster_assignments
  return(outs)
  #return(cluster_assignments)
}
