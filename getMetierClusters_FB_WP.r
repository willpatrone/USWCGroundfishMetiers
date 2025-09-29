################################################################################
#  STEP 3 OF THE MULTIVARIATE CLASSIFICATION :                                 #
#         RUN THE CLUSTERING OF THE LOGEVENTS                                  #
#         4 METHODS ARE AVALAIBLE : HAC / KMEANS / PAM / CLARA                 #
################################################################################


getMetierClusters = function(datSpecies,datLog,analysisName="",methMetier="clara",param1="euclidean",param2=NULL){
  
  # Load the table linking 3A-CODE (FAO CODE of species) to the species assemblage (level 5).
  data(correspLevel7to5)
  require(lattice)
  
  LE_ID=rownames(datSpecies)
  nbSpec=dim(datSpecies)[2]
  datSpecies=as.matrix(datSpecies,ncol=nbSpec,nrow=length(LE_ID))
  
  print("######## STEP 3 CLUSTERING ########")
  
  t1=Sys.time()
  print(paste(" --- selected method :",methMetier, " ---"))
  
        ########################################################################################################################################   CLARA
        
        if(methMetier=="clara"){
          nbLog=nrow(datLog)
          propSample=0.1
          
          # Calculation of optimal k thanks to the silhouette (second maximum)
          clustersClara.silcoeff=numeric()
          clustersClara.silcoeff[1]=0
          clustersClara.silcoeff[2]=0
          clustersClara.silcoeff[3]=0
          k=2
          compMax=1
          repeat{
            k=k+2
            print(k)
            clustersClara=clara(datLog, k, metric=param1, stand=FALSE, samples=5, sampsize=min(nbLog,round(propSample*nbLog+10*k)))
            clustersClara.silcoeff[k]=clustersClara$silinfo$avg.width
            clustersClara=clara(datLog, k+1, metric=param1, stand=FALSE, samples=5, sampsize=min(nbLog,round(propSample*nbLog+10*(k+1))))
            clustersClara.silcoeff[k+1]=clustersClara$silinfo$avg.width
            if((clustersClara.silcoeff[k-2]<clustersClara.silcoeff[k-1] & clustersClara.silcoeff[k-1]>clustersClara.silcoeff[k]) & compMax<=2){
              if(compMax==2){
                nbClust=k-1
                print(paste("2e max =",k-1))
                print(paste("nbClust =",nbClust))
                break
              } else {
                compMax=compMax+1
                print(paste("compMax1 =",compMax))
                print(paste("1er max =",k-1))
              }
            }
            if((clustersClara.silcoeff[k-1]<clustersClara.silcoeff[k] & clustersClara.silcoeff[k]>clustersClara.silcoeff[k+1]) & compMax<=2){
              if(compMax==2){
                nbClust=k
                print(paste("2e max =",k))
                print(paste("nbClust =",nbClust))
                break
              } else {
                compMax=compMax+1
                print(paste("compMax2 =",compMax))
                print(paste("1er max =",k))
              }
            }
            Store(objects())
            gc(reset=TRUE)
          }
          
          
          png(paste(analysisName,"Silhouette of the classification.png",sep="_"), width = 1200, height = 800)
          plot(clustersClara.silcoeff, main="Silhouette of the classification", xlab="Number of clusters", ylab="Silhouette")               # k optimal corresponds to maximum of silhouette's coefficients
          dev.off()
          
          Store(objects())
          gc(reset=TRUE)
          
          cat("ClaraSilCoeff",clustersClara.silcoeff,"\n")
          
          
          # CLARA with optimal k
          clusters=clara(datLog, nbClust, metric=param1, stand=FALSE, samples=5, sampsize=min(nbLog,round(propSample*nbLog+10*nbClust)))  # CLARA with optimal k
          summary(clusters)
          
          
          # # Within and between variance of clusters and classification
          # centerOfGravityClassif=numeric()
          # withinVarClusters=numeric()
          # sizeClusti=numeric()
          # centerOfGravityDatLog=colMeans(datLog)
          # centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
          # for(k in 1:nbClust){  # Within variance by cluster
          #   
          #   clusti=datLog[which(clusters$clustering==k),]
          #   if(length(which(clusters$clustering==k))==1)  centerOfGravityClusti=clusti
          #   else centerOfGravityClusti=colMeans(clusti)
          #   centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
          #   sizeClusti[k]=length(which(clusters$clustering==k))
          #   if(length(which(clusters$clustering==k))==1)  withinVarClusters[k]=0
          #   else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))
          #   
          # }
          # # Between variance
          # classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
          # # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
          # withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
          # betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
          # 
   
           
         #  # Compute the test-values for species
         #  resval=test.values(clusters$cluster,datSpecies)
         #  # Determine the target species
         #  target=targetspecies(resval)
         #  rownames(target$tabnomespcib)=paste("Cluster",1:nbClust)
         #  
         #  
         #  # Compute the percentage of logevents catching each species by cluster
         #  mainSpecies=colnames(datSpecies)
         #  percLogevents=matrix(0,ncol=length(mainSpecies),nrow=nbClust,dimnames=list(paste("Cluster ",1:nbClust,sep=""),mainSpecies))
         #  for(i in 1:nbClust){
         #    percLogevents[i,]=round(sapply(mainSpecies,function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)
         #  }
         #  
         #  
         #  # Projections on the first factorial plans
         #  png(paste(analysisName,"Projections.png",sep="_"), width = 1200, height = 800)
         #  op <- par(mfrow=c(2,3))
         #  plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
         #  if(dim(datLog)[2]>2) {
         #    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
         #    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
         #    if(dim(datLog)[2]>3) {
         #      plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
         #      plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
         #      plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
         #    }}
         #  par(op)
         #  dev.off()
         #  
         #  # Catch profile of the dataset
         #  meanprofile=colMeans(datSpecies)
         #  png(paste(analysisName,"Catch profile of the dataset.png",sep="_"), width = 1200, height = 800)
         #  op <- par(las=2)
         #  barplot(meanprofile, main="Catch profile of the dataset", xlab="Species", ylab="Percentage of catch")
         #  par(op)
         #  mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
         #  dev.off()
         #  
         #  
         #  # Catch profile by cluster
         #  nbSpec=ncol(datSpecies)
         #  summaryClusters=array(0,dim=c(6,nbSpec,nbClust))
         #  dimnames(summaryClusters)[[1]]=c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
         #  dimnames(summaryClusters)[[2]]=names(meanprofile)
         #  dimnames(summaryClusters)[[3]]=paste("Cluster",1:nbClust)
         #  for(i in 1:nbClust){
         #    if(clusters$clusinfo[i,1]==1){
         #      summaryClusters[,,i]=apply(t(as.matrix(datSpecies[which(clusters$clustering==i),])),2,
         #                                 function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
         #    }else{
         #      summaryClusters[,,i]=apply(datSpecies[which(clusters$clustering==i),],2,
         #                                 function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
         #    }
         #  }
         #  # Species names for catch profile plots
         #  nameSpPlot=character()
         #  catchMeanThreshold=2
         #  for(i in 1:nbClust){
         #    namSpi=names(meanprofile[which(t(summaryClusters["Mean",,i])>catchMeanThreshold)])
         #    numSpi=which(t(summaryClusters["Mean",,i])>catchMeanThreshold)
         #    nameSpPloti=rep("",nbSpec)
         #    nameSpPloti[numSpi]=namSpi
         #    nameSpPlot=rbind(nameSpPlot,nameSpPloti)
         #  }
         #  # Plot
         #  png(paste(analysisName,"Catch profile by cluster.png",sep="_"), width = 1200, height = 800)
         #  op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
         #  for(i in 1:nbClust){
         #    op2 <- par(las=2)
         #    barplot(t(summaryClusters["Mean",,i]), names.arg=nameSpPlot[i,], xlab="Species", ylab="Percentage of catch", col="gray")
         #    par(op2)
         #    mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
         #  }
         #  par(op)
         #  title(main=paste("Catch profile by cluster","\n","\n",sep=""))
         #  dev.off()
         #  
         #  # Standard deviation profile by cluster
         #  sdprofil=matrix(0,nrow=nbClust,ncol=nbSpec)
         #  namSdPlot=character()
         #  SdThreshold=5
         #  for(i in 1:nbClust){
         #    if(length(which(clusters$clustering==i))==1){ sdprofilclusti=rep(0,nbSpec)
         #    }else{sdprofilclusti=apply(datSpecies[which(clusters$clustering==i),],2,sd)}
         #    namSDi=names(which(sdprofilclusti>SdThreshold))
         #    numSDi=which(sdprofilclusti>SdThreshold)
         #    namSdPloti=rep("",nbSpec)
         #    namSdPloti[numSDi]=namSDi
         #    sdprofil[i,]=sdprofilclusti
         #    namSdPlot=rbind(namSdPlot,namSdPloti)
         #  }
         #  rownames(sdprofil) <- 1:nrow(sdprofil)
         #  png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
         #  op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
         #  for(i in 1:nbClust){
         #    op2 <- par(las=2)
         #    barplot(sdprofil[i,], names.arg=namSdPlot[i,], xlab="Species", ylab="Percentage of catch")
         #    par(op2)
         #    mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
         #  }
         #  par(op)
         #  title(main=paste("Standard deviation profile by cluster","\n","\n",sep=""))
         #  dev.off()
         #  
         #  
         #  # Number of Logevents by cluster
         #  x=c(1:nbClust)
         #  png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
         #  coord=barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
         #  barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
         #  text(coord,clusters$clusinfo[,1]+5,clusters$clusinfo[,1],font=2,xpd=NA)
         #  dev.off()
         #  
         #  
         #  # Profile of test-values by cluster
         #  targetresval=matrix(0,nrow=nbClust,ncol=nbSpec)
         #  colnames(targetresval)=colnames(datSpecies)
         #  rownames(targetresval)=1:nbClust
         #  nameTargetPlot=matrix(NA,nrow=nbClust,ncol=nbSpec)
         #  for(i in 1:nbClust){
         #    nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
         #    numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
         #    nameTargetPlot[i,numtargeti]=nomtargeti
         #    targetresval[i,numtargeti]=resval[nomtargeti,i]
         #  }
         #  
         #  png(paste(analysisName,"Profile of test-values by cluster.png",sep="_"), width = 1200, height = 800)
         #  op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
         #  for(i in 1:nbClust){
         #    op2 <- par(las=2)
         #    barplot(targetresval[i,],names.arg=nameTargetPlot[i,], xlab="Species", ylab="Test-value")
         #    par(op2)
         #    mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
         #  }
         #  par(op)
         #  title(main=paste("Profile of test-values by cluster","\n","\n",sep=""))
         #  dev.off()
         #  
         #  #
         #  # # Descriptive tables of the clusters
         #  clusterDesc=matrix(0,nrow=9,ncol=nbClust)
         #  for(i in 1:nbClust){
         #    clusterDesc[,i]=c(i,
         #                      length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<50))+1,
         #                      length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<90))+1,
         #                      length(which(t(summaryClusters["Median",,i])>50)),
         #                      length(which(resval[,i]>1.96)),
         #                      length(which(resval[,i]>3.29)),
         #                      length(which(apply(datSpecies,2,function (x) (clusters$clusinfo[i,1]-length(which(x[clusters$clustering==i]==0)))/clusters$clusinfo[i,1]*100)>50)),
         #                      length(which(apply(datSpecies,2,function (x) (clusters$clusinfo[i,1]-length(which(x[clusters$clustering==i]==0)))/clusters$clusinfo[i,1]*100)>90)),
         #                      clusters$clusinfo[i,1])
         #  }
         #  rownames(clusterDesc)=c("Number of species",
         #                          "to have 50% of catch", "to have 90% of catch",
         #                          "with a median higher than 50",
         #                          "with a test-value > 1.96", "with a test-value > 3.29",
         #                          "catch in 50% of the logevents", "catch in 90% of the logevents",
         #                          "Clusters size")
         #  colnames(clusterDesc)=1:nbClust
         #  clusterDesc2=as.data.frame(clusterDesc)
         # 
         #  #
         # 
         #  #
         #  # # Summary tables of the clusters
         #  namesSpecies=matrix(NA,nrow=nbClust,ncol=10)
         #  namesCapt=matrix(NA,nrow=nbClust,ncol=5)
         #  nbSpeciesCatch = min(5,dim(t(summaryClusters["Mean",,]))[[2]])
         #  namesTarget=matrix(NA,nrow=nbClust,ncol=5)
         #  nbSpeciesVT = min(5,dim(target$tabnomespcib)[[2]])
         #  tabLibname=matrix(NA,nrow=nbClust,ncol=10)
         #  listLibname=list()
         # 
         #  for(i in 1:nbClust){
         #    namesCapt[i,]=colnames(t(summaryClusters["Mean",,i]))[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)][1:nbSpeciesCatch]
         #    a=as.data.frame(t(summaryClusters["Mean",target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])],i]))
         #    colnames(a)= target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])]
         #    if(length(a)!=0){
         #      namesTarget[i,1:length(target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])])]=colnames(a[order(a,decreasing=TRUE)])
         #    }
         #    namesSpecies[i,1:length(union(namesCapt[i,],namesTarget[i,]))]=union(namesCapt[i,],namesTarget[i,])
         #  }
         # 
         #  for(i in 1:nbClust){
         #    listLibname[[i]]=lapply(as.list(namesSpecies[i,]), function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0) "NA"
         #                            else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"French_name"])
         #    tabLibname[i,]=unlist(lapply(listLibname[[i]], function(x) as.character(unlist(x))))
         #  }
         # 
         #  tabPropCatch=matrix(NA,nrow=nbClust,ncol=10)
         #  tabTestVal=matrix(NA,nrow=nbClust,ncol=10)
         #  tabPropLog=matrix(NA,nrow=nbClust,ncol=10)
         # 
         #  for(i in 1:nbClust){
         #    print("-----------------------------------------------------------------")
         #    print(paste("Cluster",i))
         #    propCatch=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
         #    tabPropCatch[i,1:length(propCatch)]=propCatch
         #    print(propCatch)
         #    testVal=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) resval[x,i]),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
         #    tabTestVal[i,1:length(testVal)]=testVal
         #    print(testVal)
         #    propLog=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)>=0.1)]
         #    tabPropLog[i,1:length(propLog)]=propLog
         #    print(propLog)
         #  }
         # 
         #  tabClusters=array(0,dim=c(10,5,nbClust))
         #  dimnames(tabClusters)[[2]]=c("Libname","FAO","Test-value","% Catch","% Logevents")
         #  dimnames(tabClusters)[[3]]=paste("Cluster",1:nbClust)
         #  for(i in 1:nbClust){
         #    tabClusters[,,i]=cbind(tabLibname[i,],namesSpecies[i,],tabTestVal[i,],tabPropCatch[i,],tabPropLog[i,])
         #  }
         # 
         #  sizeTabClusters=numeric()
         #  for(i in 1:nbClust){
         #    sizeTabClusters[i]=min(length(namesSpecies[i,!is.na(namesSpecies[i,])]),length(tabPropCatch[i,!is.na(tabPropCatch[i,])]),length(tabTestVal[i,!is.na(tabTestVal[i,])]),length(tabPropLog[i,!is.na(tabPropLog[i,])]))
         #  }
         # 
         # 
         # ## Target Species
         # ## Intersection of species from tabClusters having : - % Cumulated Catch > thresholdCatch
         #                                               ##     - Test-value > thresholdTestValue
         #  #                                                   - % Logevents > thresholdLogevents
         #  thresholdCatch=75
         #  thresholdTestValue=3
         #  thresholdLogevents=30
         # 
         #  sppCumCatch=list()
         #  sppTestValue=list()
         #  sppLogevents=list()
         #  targetSpeciesByCluster=list()
         # 
         #  for (i in 1:nbClust){
         #    percCatchCum=cumsum(as.numeric(tabClusters[,"% Catch",i]))
         #    nbSpSel=length(which(percCatchCum<thresholdCatch))+1
         #    sppCumCatch[[i]]=tabClusters[1:nbSpSel,"FAO",i]
         # 
         #    sppTestValue[[i]]=tabClusters[which(as.numeric(tabClusters[,"Test-value",i])>thresholdTestValue),"FAO",i]
         # 
         #    sppLogevents[[i]]=tabClusters[which(as.numeric(tabClusters[,"% Logevents",i])>thresholdLogevents),"FAO",i]
         # 
         #    targetSpeciesByCluster[[i]]=intersect(sppCumCatch[[i]],sppTestValue[[i]])
         #    targetSpeciesByCluster[[i]]=intersect(targetSpeciesByCluster[[i]],sppLogevents[[i]])
         #  }
         # 
         #  # List of metiers (level 7)
         #  listMetiersL7=list()
         #  for (i in 1:nbClust){
         #    metiersClusteri=targetSpeciesByCluster[[i]]
         #    metiersClusteri=as.character(unique(unlist(metiersClusteri)))
         #    metiersClusteri=paste(unlist(strsplit(metiersClusteri," ")),collapse=" ")
         #    listMetiersL7[[i]]=metiersClusteri
         #  }
         # 
         #  # Metier (level 7) of each logevent
         #  metierByLogeventL7=unlist(sapply(clusters$clustering,function(x) listMetiersL7[[x]]))
         # 
         # 
         # 
         #  # Create csv tables
         #  write.table(clusterDesc2,file="descClusters.csv",col.names=FALSE,sep=";")
         # 
         #  dfClust=data.frame()
         #  dfClust=paste("Clust ",1:nbClust,sep="")
         #  for(i in 1:nbClust){
         #    write.table(dfClust[i],file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")
         #    tabClusti=as.data.frame(tabClusters[1:sizeTabClusters[i],,i])
         #    write.table(tabClusti,file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")
         #  }
         # 
         # 
         #  LE_ID_clust=data.frame(LE_ID=LE_ID,clust=metierByLogeventL7)
         #  print(" --- end of step 3 ---")
         #  print(Sys.time()-t1)
         # 
         #  return(list(LE_ID_clust=LE_ID_clust, clusters=clusters,
         #              betweenVarClassifOnTot=betweenVarClassifOnTot, nbClust=nbClust,
         #              summaryClusters=summaryClusters, testValues=resval,
         #              testValuesSpecies=target$tabnomespcib, percLogevents=percLogevents,
         #              descClusters=clusterDesc2, tabClusters=tabClusters,
         #              targetSpecies=targetSpeciesByCluster))

        }  else stop("methMetier must be hac, kmeans, pam or clara")
  # end of the methods
  
  
} # end of the function "getMetierClusters"
