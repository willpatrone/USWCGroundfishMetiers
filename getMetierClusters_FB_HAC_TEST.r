################################################################################
#  STEP 3 OF THE MULTIVARIATE CLASSIFICATION :                                 #
#         RUN THE CLUSTERING OF THE LOGEVENTS                                  #
#         4 METHODS ARE AVALAIBLE : HAC / KMEANS / PAM / CLARA                 #
################################################################################


HgetMetierClusters = function(datSpecies,datLog,analysisName="",methMetier="hac",param1="euclidean",param2='ward'){
  
  # Load the table linking 3A-CODE (FAO CODE of species) to the species assemblage (level 5).
  data(correspLevel7to5)
  require(lattice)
  require(MASS)
  
  LE_ID=rownames(datSpecies)
  nbSpec=ncol(datSpecies)
  datSpecies=as.matrix(datSpecies,ncol=nbSpec,nrow=length(LE_ID))
  
  print("######## STEP 3 CLUSTERING ########")
  
  t1=Sys.time()
  print(paste(" --- selected method :",methMetier, " ---"))
  
  
  ########################################################################################################################################   HAC
  
  if(methMetier=="hac"){
    datSpecies =as.data.frame(datSpecies)
    datLog = as.data.frame(datLog)
    
    classifWithinVar=numeric()
    classifBetweenVar=numeric()
    classifQuality=numeric()
    sampleList=numeric()
    mProfilSample=numeric()
    classifVarExplain=numeric()
    
    nbLog=nrow(datLog)
    nbDim=ncol(datLog)
    
    # Center of gravity of datLog
    centerOfGravityDatLog=colMeans(datLog)
    
    # HAC like CLARA (HAC on sample, affectation of each logevent to a cluster, quality of classification, do it 5 times, choose the sample which gives the best quality of classification)
    print("hac on subsets...")
    
    for(i in 1:5){
      
      numSample=i
      print(paste("sample",i))
      # Sample of size 10000 logevents or 30% of all logevents
      #########minsam=min(nbLog,max(10000,round(nbLog*30/100)))
      minsam=min(round(nbLog*0.7), max(1000,round(nbLog*30/100)))
      sam=sample(1:nbLog,size=minsam,replace=FALSE)
      # Record the 5 samples
      sampleList=rbind(sampleList,sam)
      outofsam=setdiff(1:nbLog,sam)
      sampleDatLog=datLog[sam,]
      sampleDatSpecies=datSpecies[sam,]
      
      # HAC on the sample
      log.hac=hcluster(sampleDatLog, method=param1, link=param2)
      inerties.vector=log.hac$height[order(log.hac$height,decreasing=TRUE)]
      nbClust=which(scree(inerties.vector)[,"epsilon"]<0)[3]
      
      # Cut the dendogram at the selected level
      sampleClusters=cutree(log.hac,k=nbClust)
      
      Store(objects())
      gc(reset=TRUE)
      
      # Add the cluster to each logevent of the sample
      sampleDatLogWithClusters=cbind(sampleDatLog,sampleClusters)
      
      sampleClusters=sampleDatLogWithClusters[,ncol(sampleDatLogWithClusters)]
      
      # Within and between variance of clusters and classification
      centerOfGravityClassif=numeric()
      withinVarClusters=numeric()
      sizeClusti=numeric()
      centerOfGravitySampleDatLog=colMeans(sampleDatLog)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravitySampleDatLog)
      for(k in 1:nbClust){  # Within variance by cluster
        
        clusti=sampleDatLogWithClusters[which(sampleClusters==k),1:nbDim]
        if(length(which(sampleClusters==k))==1)  centerOfGravityClusti=clusti
        else centerOfGravityClusti=colMeans(clusti)
        centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
        sizeClusti[k]=length(which(sampleClusters==k))
        if(length(which(sampleClusters==k))==1)  withinVarClusters[k]=0
        else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))
        
      }
      # Between variance
      classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
      # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
      withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
      betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
      classifVarExplain=c(classifVarExplain,betweenVarClassifOnTot)
      
      
      # Catch profile by cluster for each sample
      nbSpec=ncol(datSpecies)
      mprofil=numeric()
      blank=rep(00,nbSpec)
      for(k in 1:nbClust){
        mprofilclusti=colMeans(sampleDatSpecies[which(sampleClusters==k),])
        mprofil=rbind(mprofil,mprofilclusti)
      }
      mprofil=rbind(mprofil,blank)
      
      mProfilSample=rbind(mProfilSample,mprofil)
      
      
      # Graphics
      
      # Calculation of each cluster size
      sizeClusters=numeric()
      for(k in 1:nbClust){
        sizeClusters[k]=length(which(sampleClusters==k))
      }
      
      # Compute the test-values for species
      resval=test.values(sampleClusters,sampleDatSpecies)
      # Determine the target species
      target=targetspecies(resval)
      
      
      # Projections on the first factorial plans
      png(paste(analysisName,numSample,"Sample_Projections.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(2,3))
      plot(sampleDatLog[,1], sampleDatLog[,2], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
      if(dim(datLog)[2]>2) {
        plot(sampleDatLog[,2], sampleDatLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
        plot(sampleDatLog[,1], sampleDatLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
        if(dim(datLog)[2]>3) {
          plot(sampleDatLog[,1], sampleDatLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
          plot(sampleDatLog[,2], sampleDatLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
          plot(sampleDatLog[,3], sampleDatLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
        }}
      par(op)
      dev.off()
      
      
      # Catch profile of the dataset
      meanprofile=colMeans(sampleDatSpecies)
      png(paste(analysisName,numSample,"Sample_Catch profile of the sample.png",sep="_"), width = 1200, height = 800)
      op <- par(las=2)
      barplot(meanprofile, main="Catch profile of the sample", xlab="Species", ylab="Percentage of catch")
      par(op)
      mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      dev.off()
      
      
      # Catch profile by cluster
      nbSpec=ncol(sampleDatSpecies)
      summarySampleClusters=array(0,dim=c(6,nbSpec,nbClust))
      dimnames(summarySampleClusters)[[1]]=c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
      dimnames(summarySampleClusters)[[2]]=names(meanprofile)
      dimnames(summarySampleClusters)[[3]]=paste("Cluster",1:nbClust)
      for(k in 1:nbClust){
        if(sizeClusters[k]==1){
          summarySampleClusters[,,k]=apply(t(as.matrix(sampleDatSpecies[which(sampleClusters==k),])),2,
                                           function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
        }else{
          summarySampleClusters[,,k]=apply(sampleDatSpecies[which(sampleClusters==k),],2,
                                           function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
        }
      }
      # Species names for catch profile plots
      nameSpPlot=character()
      catchMeanThreshold=2
      for(k in 1:nbClust){
        namSpi=names(meanprofile[which(t(summarySampleClusters["Mean",,k])>catchMeanThreshold)])
        numSpi=which(t(summarySampleClusters["Mean",,k])>catchMeanThreshold)
        nameSpPloti=rep("",nbSpec)
        nameSpPloti[numSpi]=namSpi
        nameSpPlot=rbind(nameSpPlot,nameSpPloti)
      }
      png(paste(analysisName,numSample,"Sample_Mean profile by cluster of the sample.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
      for(k in 1:nbClust){
        op2 <- par(las=2)
        barplot(t(summarySampleClusters["Mean",,k]), names.arg=nameSpPlot[k,], xlab="Species", ylab="Percentage of catch", col="gray")
        par(op2)
        mtext(paste("Cluster",k), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      }
      par(op)
      title(main=paste("Catch profile by cluster of the sample","\n","\n",sep=""))
      dev.off()
      
      
      # Standard deviation profile by cluster
      sdprofil=matrix(0,nrow=nbClust,ncol=nbSpec)
      namSdPlot=character()
      SdThreshold=2
      for(k in 1:nbClust){
        if(length(which(sampleClusters==k))==1){ sdprofilclusti=rep(0,nbSpec)
        }else{sdprofilclusti=apply(sampleDatSpecies[which(sampleClusters==k),],2,sd)}
        namSDi=names(which(sdprofilclusti>SdThreshold))
        numSDi=which(sdprofilclusti>SdThreshold)
        namSdPloti=rep("",nbSpec)
        namSdPloti[numSDi]=namSDi
        sdprofil[k,]=sdprofilclusti
        namSdPlot=rbind(namSdPlot,namSdPloti)
      }
      rownames(sdprofil) <- 1:nrow(sdprofil)
      png(paste(analysisName,numSample,"Sample_Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
      for(k in 1:nbClust){
        op2 <- par(las=2)
        barplot(sdprofil[k,], names.arg=namSdPlot[k,], xlab="Species", ylab="Percentage of catch")
        par(op2)
        mtext(paste("Cluster",k), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      }
      par(op)
      title(main=paste("Standard deviation profile by cluster","\n","\n",sep=""))
      dev.off()
      
      
      # Number of Logevents by cluster
      x=c(1:nbClust)
      png(paste(analysisName,numSample,"Sample_Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
      coord=barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
      barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
      text(coord,sizeClusters+100,sizeClusters,font=2,xpd=NA)
      dev.off()
      
      
      # Target Species profiles (test-value)
      targetresval=numeric()
      nameTargetPlot=character()
      for(k in 1:nbClust){
        nomtargeti=as.character(target$tabnomespcib[k,which(!is.na(target$tabnumespcib[k,]))])
        numtargeti=as.numeric(target$tabnumespcib[k,which(!is.na(target$tabnumespcib[k,]))])
        nameTargetPloti=rep("",nbSpec)
        nameTargetPloti[numtargeti]=nomtargeti
        nameTargetPlot=rbind(nameTargetPlot,nameTargetPloti)
        targetresvalclusti=rep(0,nbSpec)
        targetresvalclusti[numtargeti]=resval[nomtargeti,k]
        targetresval=rbind(targetresval,targetresvalclusti)
      }
      
      png(paste(analysisName,numSample,"Sample_Profile of target species by cluster.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
      for(k in 1:nbClust){
        op2 <- par(las=2)
        barplot(targetresval[k,],names.arg=nameTargetPlot[k,], cex.names=1, xlab="Species", ylab="Test-value")
        par(op2)
        mtext(paste("Cluster",k), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      }
      par(op)
      title(main=paste("Profile of target species by cluster","\n","\n",sep=""))
      dev.off()
      
      Store(objects())
      gc(reset=TRUE)
      
    } # end of for(i in 1:5)
    
    
    
    
    # Select the sample which gives the smaller classification's quality (the best sample)
    sam=sampleList[which.max(classifVarExplain),]
    outofsam=setdiff(1:nbLog,sam)
    sampleDatLog=datLog[sam,]
    
    nbLogSample=nrow(sampleDatLog)
    nbDim=ncol(sampleDatLog)
    
    
    # HAC with the best sample
    print("Final HAC")
    log.hac=hcluster(sampleDatLog, method=param1, link=param2)
    

    # Determine the number of cluster thanks to the scree-test
    inerties.vector=log.hac$height[order(log.hac$height,decreasing=TRUE)]
    nbClust=which(scree(inerties.vector)[,"epsilon"]<0)[3]

    # Create clusters
    sampleClusters <- cutree(log.hac, k=nbClust)
    sampleClusters <- as.factor(sampleClusters)

    
    min_samples_needed <- ncol(sampleDatLog) + 1
    cluster_sizes <- table(sampleClusters)
    valid_clusters <- names(cluster_sizes)[cluster_sizes >= min_samples_needed]

    cat("Original clusters:", length(unique(sampleClusters)), "\n")
    cat("Valid clusters after filtering:", length(valid_clusters), "\n")

    if(length(valid_clusters) < 2) {
      stop("Not enough valid clusters for LDA. Try reducing the number of clusters or increasing sample size.")
    }

    # Filter data to keep only valid clusters
    valid_indices <- sampleClusters %in% valid_clusters
    sampleDatLog_filtered <- sampleDatLog[valid_indices, ]
    sampleClusters_filtered <- droplevels(sampleClusters[valid_indices])

    # Remove constant variables
    var_check <- apply(sampleDatLog_filtered, 2, var, na.rm=TRUE)
    non_constant_vars <- var_check > 1e-10
    sampleDatLog_filtered <- sampleDatLog_filtered[, non_constant_vars]

    cat("Variables after removing constants:", ncol(sampleDatLog_filtered), "\n")

    # Create final dataset for LDA
    sampleDatLogWithClusters <- cbind(sampleDatLog_filtered, sampleClusters_filtered)
    colnames(sampleDatLogWithClusters)[ncol(sampleDatLogWithClusters)] <- "sampleClusters"
  
     # Try LDA with error handling
     tryCatch({
       learning <- lda(sampleClusters ~ ., data=sampleDatLogWithClusters)
       cat("LDA successful!\n")
     }, error = function(e) {
       cat("LDA failed with error:", e$message, "\n")
       stop("LDA could not be computed. Check data structure.")
     })
  
     # # Update the other data to match filtered variables
     # otherLog <- datLog[outofsam, non_constant_vars]
     # otherLog <- as.data.frame(otherLog)
     # 
     # 
     # # Predict with error handling
     # tryCatch({
     #   pred <- predict(learning, otherLog)
     #   otherDatLogWithClusters <- cbind(otherLog, pred$class)
     #   colnames(otherDatLogWithClusters) <- colnames(sampleDatLogWithClusters)
     # }, error = function(e) {
     #   cat("Prediction failed with error:", e$message, "\n")
     #   stop("Prediction could not be computed.")
     # })
  
    if(length(outofsam) > 0) {
      otherLog <- datLog[outofsam, non_constant_vars]
      otherLog <- as.data.frame(otherLog)
      
      cat("otherLog dimensions:", dim(otherLog), "\n")
      tryCatch({
        pred <- predict(learning, otherLog)
      
        cat("pred$class levels:", levels(pred$class), "\n")
        cat("pred$class unique values:", unique(as.numeric(pred$class)), "\n")
        cat("Expected range: 1 to", nbClust, "\n")
        
        pred$class <- factor(pred$class, labels = 1:nlevels(pred$class))
        
        pred_numeric <- as.numeric(pred$class)
        otherLog_df <- as.data.frame(otherLog)
        
        otherDatLogWithClusters <- data.frame(otherLog_df, sampleClusters = pred_numeric)
        
      }, error = function(e) {
        cat("Prediction failed with error:", e$message, "\n")
        stop("Prediction could not be computed.")
      })
    
     # Rebuild complete clusters - accounting for filtered clusters
     clusters <- numeric(length=nbLog)
     clusters[sam[valid_indices]] <- as.numeric(sampleClusters_filtered)
     clusters[outofsam] <- as.numeric(pred$class)
     nbClust <- length(valid_clusters)
  
     # Within and between variance of clusters and classification
     centerOfGravityClassif=numeric()
     withinVarClusters=numeric()
     sizeClusti=numeric()
     centerOfGravityDatLog=colMeans(datLog)
     centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
    
    for(k in 1:nbClust){  # Within variance by cluster
      clusti=datLog[which(clusters==k),1:nbDim]
    
       
       if(length(which(clusters==k))==1)  centerOfGravityClusti=clusti
       else centerOfGravityClusti=colMeans(clusti)
       centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
       sizeClusti[k]=length(which(clusters==k))
       if(length(which(clusters==k))==1)  withinVarClusters[k]=0
       else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))
       
     }
     # Between variance
    classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
     
     # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
     withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
     betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
  
  
     # Calculation of each cluster size
    n=nrow(datLog)
    sizeClusters=numeric()
    for(k in 1:nbClust){
       sizeClusters[k]=length(which(clusters==k))
     }
  
     # Compute the test-values for species
     resval=test.values(clusters,datSpecies)
     # Determine the target species
     target=targetspecies(resval)
     rownames(target$tabnomespcib)=paste("Cluster",1:nbClust)
  
  
     # Compute the percentage of logevents catching each species by cluster
     mainSpecies=colnames(datSpecies)
     percLogevents=matrix(0,ncol=length(mainSpecies),nrow=nbClust,dimnames=list(paste("Cluster ",1:nbClust,sep=""),mainSpecies))
     for(i in 1:nbClust){
       percLogevents[i,]=round(sapply(mainSpecies,function(x) (sizeClusters[i]-length(which(datSpecies[clusters==i,x]==0)))/sizeClusters[i]*100),digits=1)
     }
   }
   # end of the methods
  }
  return(list(
    clusters = clusters,
    nbClust = nbClust,
    sizeClusters = sizeClusters,
    resval = resval,
    target = target,
    percLogevents = percLogevents,
    learning = learning,
  ))
} # end of the function "getMetierClusters"
