BFN3<-function(DATA,nrreplics,updown,nrruns,cutvalue)
{
  cutval <- cutvalue
	nrprobesets<-nrow(DATA)
	nrarrays<-ncol(DATA)
	nrconditions<-length(nrreplics)
	target <- matrix(rep(0,nrprobesets*3),nrow=nrprobesets,ncol=3)

	#Preliminary: Calculation of condition-specific averages for all probesets
	lower<-1
	AV1<-matrix(c(rep(0,nrprobesets*nrarrays)),nrow=nrprobesets,ncol=nrarrays)
	AV2<-matrix(c(rep(0,nrprobesets*nrconditions)),nrow=nrprobesets,ncol=nrconditions)
	AV3<-matrix(c(rep(0,nrprobesets*(nrconditions-3))),nrow=nrprobesets,ncol=(nrconditions-3))
	for (i in (1:nrconditions)){
    		upper<-lower+nrreplics[i]-1
    		v<-rep(1,nrreplics[i])
    		average<-(DATA[,lower:upper])%*%(v*(t(v)%*%v)^(-1))
    		AV1[,lower:upper]<-average%*%t(v)
    		AV2[,i]<-average
    		lower<-upper+1
	}

	#1. Initial values
	sigma2_0<-rowSums((AV1-DATA)^2)/nrarrays
	SIGMA2_SAMPLED<-sigma2_0%*%t(rep(1,nrconditions))
  s_d <- matrix(rep(0,nrprobesets),nrow=nrprobesets)
  for(i in 1:nrprobesets){
        target[i,1] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[1],arr.ind=TRUE)[1]
        target[i,2] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[2],arr.ind=TRUE)[1]
        target[i,3] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[3],arr.ind=TRUE)[1]
        if((target[i,1] > target[i,2]) && (target[i,2] > target[i,3])){
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- temp2[1,-(target[i,3])]
             AV3[i,] <- t(temp3)
        }
        else if((target[i,1] < target[i,2]) && (target[i,2] < target[i,3])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- temp2[1,-(target[i,1])]
             AV3[i,] <- t(temp3)
        }
        else if((target[i,1] > target[i,2]) && (target[i,2] < target[i,3]) && (target[i,1] > target[i,3])){
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- temp2[1,-(target[i,2])]
             AV3[i,] <- t(temp3)
        }
        else if((target[i,1] > target[i,2]) && (target[i,2] < target[i,3]) && (target[i,1] < target[i,3])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- temp2[1,-(target[i,2])]
             AV3[i,] <- t(temp3)
        }
        else if((target[i,1] < target[i,2]) && (target[i,2] > target[i,3]) && (target[i,1] > target[i,3])){
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- temp2[1,-(target[i,3])]
             AV3[i,] <- t(temp3)
        }
        else{
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- temp2[1,-(target[i,1])]
             AV3[i,] <- t(temp3)
        }
  }
  s_d <- sd(t(AV3),na.rm=FALSE)

	#2. Gibbs Sampler
	indicsum<-rep(0,nrprobesets)
	NRREPLICS<-rep(1,nrprobesets)%*%t(nrreplics)
	MU_SAMPLED2<-matrix(rep(0,nrprobesets*nrarrays),nrow=nrprobesets)
	#SIGNDIFF <- matrix(rep(0,nrprobesets*(nrconditions-3)),nrow=nrprobesets,ncol=(nrconditions-3))
	SIGNDIFF <- matrix(rep(0,nrprobesets),nrow=nrprobesets)
	for (i in (1:nrruns)){
        #cat('Average: \n',AV2,'\n')
    		SEM<-sqrt(SIGMA2_SAMPLED/NRREPLICS)
    		#cat('Standard Deviation for run ',i,': \n',SEM,'\n')
    		MU_SAMPLED<-AV2+(SEM*matrix(rnorm(nrprobesets*nrconditions),nrow=nrprobesets))
    		if(nrruns == 1){
        }
        else{
             for(k in 1:nrprobesets){
                   if((target[k,1] > target[k,2]) && (target[k,2] > target[k,3])){
                                  temp <- AV2[k,-(target[k,1])]
                                  temp1 <- t(temp)
                                  temp2 <- t(temp1[1,-(target[k,2])])
                                  temp3 <- temp2[1,-(target[k,3])]
                                  AV3[k,] <- t(temp3)
                   }
                   else if((target[k,1] < target[k,2]) && (target[k,2] < target[k,3])){
                                  temp <- AV2[k,-(target[k,3])]
                                  temp1 <- t(temp)
                                  temp2 <- t(temp1[1,-(target[k,2])])
                                  temp3 <- temp2[1,-(target[k,1])]
                                  AV3[k,] <- t(temp3)
                  }
                  else if((target[k,1] > target[k,2]) && (target[k,2] < target[k,3]) && (target[k,1] > target[k,3])){
                                  temp <- AV2[k,-(target[k,1])]
                                  temp1 <- t(temp)
                                  temp2 <- t(temp1[1,-(target[k,3])])
                                  temp3 <- temp2[1,-(target[k,2])]
                                  AV3[k,] <- t(temp3)
                  }
                  else if((target[k,1] > target[k,2]) && (target[k,2] < target[k,3]) && (target[k,1] < target[k,3])){
                                  temp <- AV2[k,-(target[k,3])]
                                  temp1 <- t(temp)
                                  temp2 <- t(temp1[1,-(target[k,1])])
                                  temp3 <- temp2[1,-(target[k,2])]
                                  AV3[k,] <- t(temp3)
                  }
                  else if((target[k,1] < target[k,2]) && (target[k,2] > target[k,3]) && (target[k,1] > target[k,3])){
                                  temp <- AV2[k,-(target[k,2])]
                                  temp1 <- t(temp)
                                  temp2 <- t(temp1[1,-(target[k,1])])
                                  temp3 <- temp2[1,-(target[k,3])]
                                  AV3[k,] <- t(temp3)
                  }
                  else{
                                  temp <- AV2[k,-(target[k,2])]
                                  temp1 <- t(temp)
                                  temp2 <- t(temp1[1,-(target[k,3])])
                                  temp3 <- temp2[1,-(target[k,1])]
                                  AV3[k,] <- t(temp3)
                  }
              }
             s_d <- sd(t(AV3),na.rm=FALSE)
        }
    		#cat('Samples Mean per tissue for run ',i,': \n',MU_SAMPLED,'\n')
    		lower<-1
    		for (j in (1:nrconditions)){
        		upper<-lower+nrreplics[j]-1
        		MU_SAMPLED2[,lower:upper]<-MU_SAMPLED[,j]%*%t(rep(1,nrreplics[j]))
        		lower<-upper+1
    		}
    		#cat('Samples Mean tissue with replications for run ',i,': \n',MU_SAMPLED2,'\n')

    TEMP_MU_SAMPLED <- MU_SAMPLED
    for(k in (1:nrprobesets)){
        if((target[k,1] > target[k,2]) && (target[k,2] > target[k,3])){
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			      TEMP1 <- t(TEMP)
			      TEMP2 <- TEMP1[1,-(target[k,2])]
			      TEMP2 <- t(TEMP2)
			      TEMP3 <- TEMP2[1,-(target[k,3])]
			      TEMP3 <- t(TEMP3)
    	      #SIGNDIFF[k,] <-updown*sign((MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k]))%*%t(rep(1,nrconditions-3)) - TEMP3)
    	      SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k] + max(TEMP3)[1]))
    	  }
       else if((target[k,1] < target[k,2]) && (target[k,2] < target[k,3])){
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			      TEMP1 <- t(TEMP)
			      TEMP2 <- TEMP1[1,-(target[k,2])]
			      TEMP2 <- t(TEMP2)
			      TEMP3 <- TEMP2[1,-(target[k,1])]
			      TEMP3 <- t(TEMP3)
    	      #SIGNDIFF[k,] <-updown*sign((MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k]))%*%t(rep(1,nrconditions-3)) - TEMP3)
    	      SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k] + max(TEMP3)[1]))
        }
       else if((target[k,1] > target[k,2]) && (target[k,2] < target[k,3]) && (target[k,1] > target[k,3])){
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			      TEMP1 <- t(TEMP)
			      TEMP2 <- TEMP1[1,-(target[k,3])]
			      TEMP2 <- t(TEMP2)
			      TEMP3 <- TEMP2[1,-(target[k,2])]
			      TEMP3 <- t(TEMP3)
    	      #SIGNDIFF[k,] <-updown*sign((MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k]))%*%t(rep(1,nrconditions-3)) - TEMP3)
    	      SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k] + max(TEMP3)[1]))
        }
       else if((target[k,1] > target[k,2]) && (target[k,2] < target[k,3]) && (target[k,1] < target[k,3])){
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			      TEMP1 <- t(TEMP)
			      TEMP2 <- TEMP1[1,-(target[k,1])]
			      TEMP2 <- t(TEMP2)
			      TEMP3 <- TEMP2[1,-(target[k,2])]
			      TEMP3 <- t(TEMP3)
    	      #SIGNDIFF[k,] <-updown*sign((MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k]))%*%t(rep(1,nrconditions-3)) - TEMP3)
    	      SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k] + max(TEMP3)[1]))
        }
       else if((target[k,1] < target[k,2]) && (target[k,2] > target[k,3]) && (target[k,1] > target[k,3])){
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			      TEMP1 <- t(TEMP)
			      TEMP2 <- TEMP1[1,-(target[k,1])]
			      TEMP2 <- t(TEMP2)
			      TEMP3 <- TEMP2[1,-(target[k,3])]
			      TEMP3 <- t(TEMP3)
    	      #SIGNDIFF[k,] <-updown*sign((MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k]))%*%t(rep(1,nrconditions-3)) - TEMP3)
    	      SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k] + max(TEMP3)[1]))
        }
        else{
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			      TEMP1 <- t(TEMP)
			      TEMP2 <- TEMP1[1,-(target[k,3])]
			      TEMP2 <- t(TEMP2)
			      TEMP3 <- TEMP2[1,-(target[k,1])]
			      TEMP3 <- t(TEMP3)
    	      #SIGNDIFF[k,] <-updown*sign((MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k]))%*%t(rep(1,nrconditions-3)) - TEMP3)
    	      SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,3]] - (cutval * s_d[k] + max(TEMP3)[1]))
        }
  	}
		#SIGNDIFF<-updown*sign((MU_SAMPLED[,target[1]] - MU_SAMPLED[,target[2]])%*%t(rep(1,nrconditions))-MU_SAMPLED)
		#cat('Signdifference for run ',i,': \n',SIGNDIFF,'\n')
		
    #indic<-rowSums(SIGNDIFF,2)==(nrconditions-3)
    indic <- (SIGNDIFF) == 1
		indicsum<-indicsum+indic
		#cat('Indicsum for run ',i,': \n',indicsum,'\n')
		#cat('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%','\n')
		#df<-nrarrays-3
		df<-nrarrays-2
    nscale<-rowSums((MU_SAMPLED2-DATA)^2)
    chisqvals<-2*rgamma(nrprobesets, df/2, rate = 1, scale = 1)
   	sigma2_sampled<-nscale*chisqvals^(-1)
    SIGMA2_SAMPLED<-sigma2_sampled%*%t(rep(1,nrconditions))
	}
	#print(SIGNDIFF)
  #print(target)
  #print(indicsum)
	f1<-indicsum/nrruns
	f2<-1-f1
  #BFN3<-nrconditions*f1/((nrconditions/(nrconditions-3))*f2)
	BFN3<- (((nrconditions)^2+2)*(nrconditions-3) * f1) / (6*f2)
}
