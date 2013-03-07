#BAYES FACTOR CALCULATOR FOR 4-SELECTIVE CASE
BFN4<-function(DATA,nrreplics,updown,nrruns,cutvalue)
{
  cutval <- cutvalue
	nrprobesets<-nrow(DATA)
	nrarrays<-ncol(DATA)
	nrconditions<-length(nrreplics)
	#Top four high expression tissues are considered for each probeset
	target <- matrix(rep(0,nrprobesets*4),nrow=nrprobesets,ncol=4)
	#Preliminary: Calculation of condition-specific averages for all probesets
	lower<-1
	AV1<-matrix(c(rep(0,nrprobesets*nrarrays)),nrow=nrprobesets,ncol=nrarrays)
	AV2<-matrix(c(rep(0,nrprobesets*nrconditions)),nrow=nrprobesets,ncol=nrconditions)
	AV3<-matrix(c(rep(0,nrprobesets*(nrconditions-4))),nrow=nrprobesets,ncol=(nrconditions-4))
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
        target[i,4] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[4],arr.ind=TRUE)[1]
        if((target[i,1] > target[i,2]) && (target[i,2] > target[i,3]) && (target[i,3] > target[i,4])){
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- t(temp2[1,-(target[i,3])])
             temp4 <- temp3[1,-(target[i,4])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,1] < target[i,2]) && (target[i,2] < target[i,3]) && (target[i,3] < target[i,4])){
             temp <- AV2[i,-(target[i,4])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- t(temp2[1,-(target[i,2])])
             temp4 <- temp3[1,-(target[i,1])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,4] > target[i,2]) && (target[i,2] > target[i,3]) && (target[i,3] > target[i,1])){
             temp <- AV2[i,-(target[i,4])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- t(temp2[1,-(target[i,3])])
             temp4 <- temp3[1,-(target[i,1])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,2] > target[i,4]) && (target[i,4] > target[i,3]) && (target[i,3] > target[i,1])){
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,4])])
             temp3 <- t(temp2[1,-(target[i,3])])
             temp4 <- temp3[1,-(target[i,1])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,2] > target[i,3]) && (target[i,3] > target[i,4]) && (target[i,4] > target[i,1])){
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- t(temp2[1,-(target[i,4])])
             temp4 <- temp3[1,-(target[i,1])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,3] > target[i,4]) && (target[i,4] > target[i,2]) && (target[i,2] > target[i,1])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,4])])
             temp3 <- t(temp2[1,-(target[i,2])])
             temp4 <- temp3[1,-(target[i,1])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,3] > target[i,2]) && (target[i,2] > target[i,4]) && (target[i,4] > target[i,1])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- t(temp2[1,-(target[i,4])])
             temp4 <- temp3[1,-(target[i,1])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,4] > target[i,3]) && (target[i,3] > target[i,1]) && (target[i,1] > target[i,2])){
             temp <- AV2[i,-(target[i,4])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- t(temp2[1,-(target[i,1])])
             temp4 <- temp3[1,-(target[i,2])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,3] > target[i,4]) && (target[i,4] > target[i,1]) && (target[i,1] > target[i,2])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,4])])
             temp3 <- t(temp2[1,-(target[i,1])])
             temp4 <- temp3[1,-(target[i,2])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,4] > target[i,2]) && (target[i,2] > target[i,1]) && (target[i,1] > target[i,3])){
             temp <- AV2[i,-(target[i,4])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- t(temp2[1,-(target[i,1])])
             temp4 <- temp3[1,-(target[i,3])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,3] > target[i,2]) && (target[i,2] > target[i,1]) && (target[i,1] > target[i,4])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- t(temp2[1,-(target[i,1])])
             temp4 <- temp3[1,-(target[i,4])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,2] > target[i,4]) && (target[i,4] > target[i,1]) && (target[i,1] > target[i,3])){
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,4])])
             temp3 <- t(temp2[1,-(target[i,1])])
             temp4 <- temp3[1,-(target[i,3])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,2] > target[i,3]) && (target[i,3] > target[i,1]) && (target[i,1] > target[i,4])){
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- t(temp2[1,-(target[i,1])])
             temp4 <- temp3[1,-(target[i,4])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,4] > target[i,1]) && (target[i,1] > target[i,3]) && (target[i,3] > target[i,2])){
             temp <- AV2[i,-(target[i,4])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- t(temp2[1,-(target[i,3])])
             temp4 <- temp3[1,-(target[i,2])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,4] > target[i,1]) && (target[i,1] > target[i,2]) && (target[i,2] > target[i,3])){
             temp <- AV2[i,-(target[i,4])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- t(temp2[1,-(target[i,2])])
             temp4 <- temp3[1,-(target[i,3])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,3] > target[i,1]) && (target[i,1] > target[i,4]) && (target[i,4] > target[i,2])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- t(temp2[1,-(target[i,4])])
             temp4 <- temp3[1,-(target[i,2])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,2] > target[i,1]) && (target[i,1] > target[i,4]) && (target[i,4] > target[i,3])){
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- t(temp2[1,-(target[i,4])])
             temp4 <- temp3[1,-(target[i,3])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,2] > target[i,1]) && (target[i,1] > target[i,3]) && (target[i,3] > target[i,4])){
             temp <- AV2[i,-(target[i,2])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- t(temp2[1,-(target[i,3])])
             temp4 <- temp3[1,-(target[i,4])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,3] > target[i,1]) && (target[i,1] > target[i,2]) && (target[i,2] > target[i,4])){
             temp <- AV2[i,-(target[i,3])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,1])])
             temp3 <- t(temp2[1,-(target[i,2])])
             temp4 <- temp3[1,-(target[i,4])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,1] > target[i,4]) && (target[i,4] > target[i,3]) && (target[i,3] > target[i,2])){
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,4])])
             temp3 <- t(temp2[1,-(target[i,3])])
             temp4 <- temp3[1,-(target[i,2])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,1] > target[i,4]) && (target[i,4] > target[i,2]) && (target[i,2] > target[i,3])){
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,4])])
             temp3 <- t(temp2[1,-(target[i,2])])
             temp4 <- temp3[1,-(target[i,3])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,1] > target[i,2]) && (target[i,2] > target[i,4]) && (target[i,4] > target[i,3])){
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,2])])
             temp3 <- t(temp2[1,-(target[i,4])])
             temp4 <- temp3[1,-(target[i,3])]
             AV3[i,] <- t(temp4)
        }
        else if((target[i,1] > target[i,3]) && (target[i,3] > target[i,2]) && (target[i,2] > target[i,4])){
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- t(temp2[1,-(target[i,2])])
             temp4 <- temp3[1,-(target[i,4])]
             AV3[i,] <- t(temp4)
        }
        else{
             temp <- AV2[i,-(target[i,1])]
             temp1 <- t(temp)
             temp2 <- t(temp1[1,-(target[i,3])])
             temp3 <- t(temp2[1,-(target[i,4])])
             temp4 <- temp3[1,-(target[i,2])]
             AV3[i,] <- t(temp4)
        }
  }
  s_d <- sd(t(AV3),na.rm=FALSE)

	#2. Gibbs Sampler
	indicsum<-rep(0,nrprobesets)
	NRREPLICS<-rep(1,nrprobesets)%*%t(nrreplics)
	MU_SAMPLED2<-matrix(rep(0,nrprobesets*nrarrays),nrow=nrprobesets)
	SIGNDIFF <- matrix(rep(0,nrprobesets),nrow=nrprobesets)
	for (i in (1:nrruns)){
        	SEM<-sqrt(SIGMA2_SAMPLED/NRREPLICS)
    		MU_SAMPLED<-AV2+(SEM*matrix(rnorm(nrprobesets*nrconditions),nrow=nrprobesets))
    		if(nrruns == 1){
        }
        else{
             for(k in 1:nrprobesets){
                   if((target[k,1] > target[k,2]) && (target[k,2] > target[k,3]) && (target[k,3] > target[k,4])){
                                   temp <- AV2[k,-(target[k,1])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,2])])
                                   temp3 <- t(temp2[1,-(target[k,3])])
                                   temp4 <- temp3[1,-(target[k,4])]
                                   AV3[k,] <- t(temp4)
                   }
                   else if((target[k,1] < target[k,2]) && (target[k,2] < target[k,3]) && (target[k,3] < target[k,4])){
                                   temp <- AV2[k,-(target[k,4])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,3])])
                                   temp3 <- t(temp2[1,-(target[k,2])])
                                   temp4 <- temp3[1,-(target[k,1])]
                                   AV3[k,] <- t(temp4)
                   }
                   else if((target[k,4] > target[k,2]) && (target[k,2] > target[k,3]) && (target[k,3] > target[k,1])){
                                   temp <- AV2[k,-(target[k,4])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,2])])
                                   temp3 <- t(temp2[1,-(target[k,3])])
                                   temp4 <- temp3[1,-(target[k,1])]
                                   AV3[k,] <- t(temp4)
                   }
                   else if((target[k,2] > target[k,4]) && (target[k,4] > target[k,3]) && (target[k,3] > target[k,1])){
                                   temp <- AV2[k,-(target[k,2])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,4])])
                                   temp3 <- t(temp2[1,-(target[k,3])])
                                   temp4 <- temp3[1,-(target[k,1])]
                                   AV3[k,] <- t(temp4)
                   }
                   else if((target[k,2] > target[k,3]) && (target[k,3] > target[k,4]) && (target[k,4] > target[k,1])){
                                   temp <- AV2[k,-(target[k,2])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,3])])
                                   temp3 <- t(temp2[1,-(target[k,4])])
                                   temp4 <- temp3[1,-(target[k,1])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,3] > target[k,4]) && (target[k,4] > target[k,2]) && (target[k,2] > target[k,1])){
                                   temp <- AV2[k,-(target[k,3])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,4])])
                                   temp3 <- t(temp2[1,-(target[k,2])])
                                   temp4 <- temp3[1,-(target[k,1])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,3] > target[k,2]) && (target[k,2] > target[k,4]) && (target[k,4] > target[k,1])){
                                   temp <- AV2[k,-(target[k,3])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,2])])
                                   temp3 <- t(temp2[1,-(target[k,4])])
                                   temp4 <- temp3[1,-(target[k,1])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,4] > target[k,3]) && (target[k,3] > target[k,1]) && (target[k,1] > target[k,2])){
                                   temp <- AV2[k,-(target[k,4])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,3])])
                                   temp3 <- t(temp2[1,-(target[k,1])])
                                   temp4 <- temp3[1,-(target[k,2])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,3] > target[k,4]) && (target[k,4] > target[k,1]) && (target[k,1] > target[k,2])){
                                   temp <- AV2[k,-(target[k,3])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,4])])
                                   temp3 <- t(temp2[1,-(target[k,1])])
                                   temp4 <- temp3[1,-(target[k,2])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,4] > target[k,2]) && (target[k,2] > target[k,1]) && (target[k,1] > target[k,3])){
                                   temp <- AV2[k,-(target[k,4])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,2])])
                                   temp3 <- t(temp2[1,-(target[k,1])])
                                   temp4 <- temp3[1,-(target[k,3])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,3] > target[k,2]) && (target[k,2] > target[k,1]) && (target[k,1] > target[k,4])){
                                   temp <- AV2[k,-(target[k,3])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,2])])
                                   temp3 <- t(temp2[1,-(target[k,1])])
                                   temp4 <- temp3[1,-(target[k,4])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,2] > target[k,4]) && (target[k,4] > target[k,1]) && (target[k,1] > target[k,3])){
                                   temp <- AV2[k,-(target[k,2])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,4])])
                                   temp3 <- t(temp2[1,-(target[k,1])])
                                   temp4 <- temp3[1,-(target[k,3])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,2] > target[k,3]) && (target[k,3] > target[k,1]) && (target[k,1] > target[k,4])){
                                   temp <- AV2[k,-(target[k,2])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,3])])
                                   temp3 <- t(temp2[1,-(target[k,1])])
                                   temp4 <- temp3[1,-(target[k,4])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,4] > target[k,1]) && (target[k,1] > target[k,3]) && (target[k,3] > target[k,2])){
                                   temp <- AV2[k,-(target[k,4])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,1])])
                                   temp3 <- t(temp2[1,-(target[k,3])])
                                   temp4 <- temp3[1,-(target[k,2])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,4] > target[k,1]) && (target[k,1] > target[k,2]) && (target[k,2] > target[k,3])){
                                   temp <- AV2[k,-(target[k,4])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,1])])
                                   temp3 <- t(temp2[1,-(target[k,2])])
                                   temp4 <- temp3[1,-(target[k,3])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,3] > target[k,1]) && (target[k,1] > target[k,4]) && (target[k,4] > target[k,2])){
                                   temp <- AV2[k,-(target[k,3])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,1])])
                                   temp3 <- t(temp2[1,-(target[k,4])])
                                   temp4 <- temp3[1,-(target[k,2])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,2] > target[k,1]) && (target[k,1] > target[k,4]) && (target[k,4] > target[k,3])){
                                   temp <- AV2[k,-(target[k,2])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,1])])
                                   temp3 <- t(temp2[1,-(target[k,4])])
                                   temp4 <- temp3[1,-(target[k,3])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,2] > target[k,1]) && (target[k,1] > target[k,3]) && (target[k,3] > target[k,4])){
                                   temp <- AV2[k,-(target[k,2])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,1])])
                                   temp3 <- t(temp2[1,-(target[k,3])])
                                   temp4 <- temp3[1,-(target[k,4])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,3] > target[k,1]) && (target[k,1] > target[k,2]) && (target[k,2] > target[k,4])){
                                   temp <- AV2[k,-(target[k,3])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,1])])
                                   temp3 <- t(temp2[1,-(target[k,2])])
                                   temp4 <- temp3[1,-(target[k,4])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,1] > target[k,4]) && (target[k,4] > target[k,3]) && (target[k,3] > target[k,2])){
                                   temp <- AV2[k,-(target[k,1])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,4])])
                                   temp3 <- t(temp2[1,-(target[k,3])])
                                   temp4 <- temp3[1,-(target[k,2])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,1] > target[k,4]) && (target[k,4] > target[k,2]) && (target[k,2] > target[k,3])){
                                   temp <- AV2[k,-(target[k,1])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,4])])
                                   temp3 <- t(temp2[1,-(target[k,2])])
                                   temp4 <- temp3[1,-(target[k,3])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,1] > target[k,2]) && (target[k,2] > target[k,4]) && (target[k,4] > target[k,3])){
                                   temp <- AV2[k,-(target[k,1])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,2])])
                                   temp3 <- t(temp2[1,-(target[k,4])])
                                   temp4 <- temp3[1,-(target[k,3])]
                                   AV3[k,] <- t(temp4)
                  }
                  else if((target[k,1] > target[k,3]) && (target[k,3] > target[k,2]) && (target[k,2] > target[k,4])){
                                   temp <- AV2[k,-(target[k,1])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,3])])
                                   temp3 <- t(temp2[1,-(target[k,2])])
                                   temp4 <- temp3[1,-(target[k,4])]
                                   AV3[k,] <- t(temp4)
                  }
                  else{
                                   temp <- AV2[k,-(target[k,1])]
                                   temp1 <- t(temp)
                                   temp2 <- t(temp1[1,-(target[k,3])])
                                   temp3 <- t(temp2[1,-(target[k,4])])
                                   temp4 <- temp3[1,-(target[k,2])]
                                   AV3[k,] <- t(temp4)
                  }
              }
             s_d <- sd(t(AV3),na.rm=FALSE)
        }
    		lower<-1
    		for (j in (1:nrconditions)){
        		upper<-lower+nrreplics[j]-1
        		MU_SAMPLED2[,lower:upper]<-MU_SAMPLED[,j]%*%t(rep(1,nrreplics[j]))
        		lower<-upper+1
    		}
    TEMP_MU_SAMPLED <- MU_SAMPLED
    for(k in (1:nrprobesets)){
        if((target[k,1] > target[k,2]) && (target[k,2] > target[k,3]) && (target[k,3] > target[k,4])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,2])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,3])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,4])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,1] < target[k,2]) && (target[k,2] < target[k,3]) && (target[k,3] < target[k,4])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,4])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,3])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,2])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,1])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,4] > target[k,2]) && (target[k,2] > target[k,3]) && (target[k,3] > target[k,1])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,4])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,2])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,3])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,1])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,2] > target[k,4]) && (target[k,4] > target[k,3]) && (target[k,3] > target[k,1])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,4])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,3])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,1])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,2] > target[k,3]) && (target[k,3] > target[k,4]) && (target[k,4] > target[k,1])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,3])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,4])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,1])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,3] > target[k,4]) && (target[k,4] > target[k,2]) && (target[k,2] > target[k,1])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,4])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,2])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,1])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,3] > target[k,2]) && (target[k,2] > target[k,4]) && (target[k,4] > target[k,1])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,2])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,4])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,1])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,4] > target[k,3]) && (target[k,3] > target[k,1]) && (target[k,1] > target[k,2])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,4])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,3])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,1])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,2])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,3] > target[k,4]) && (target[k,4] > target[k,1]) && (target[k,1] > target[k,2])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,4])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,1])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,2])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,4] > target[k,2]) && (target[k,2] > target[k,1]) && (target[k,1] > target[k,3])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,4])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,2])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,1])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,3])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,3] > target[k,2]) && (target[k,2] > target[k,1]) && (target[k,1] > target[k,4])){
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,2])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,1])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,4])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,2] > target[k,4]) && (target[k,4] > target[k,1]) && (target[k,1] > target[k,3])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,4])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,1])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,3])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,2] > target[k,3]) && (target[k,3] > target[k,1]) && (target[k,1] > target[k,4])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,3])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,1])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,4])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,4] > target[k,1]) && (target[k,1] > target[k,3]) && (target[k,3] > target[k,2])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,4])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,1])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,3])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,2])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,4] > target[k,1]) && (target[k,1] > target[k,2]) && (target[k,2] > target[k,3])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,4])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,1])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,2])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,3])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,3] > target[k,1]) && (target[k,1] > target[k,4]) && (target[k,4] > target[k,2])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,1])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,4])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,2])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,2] > target[k,1]) && (target[k,1] > target[k,4]) && (target[k,4] > target[k,3])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,1])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,4])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,3])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,2] > target[k,1]) && (target[k,1] > target[k,3]) && (target[k,3] > target[k,4])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,2])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,1])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,3])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,4])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,3] > target[k,1]) && (target[k,1] > target[k,2]) && (target[k,2] > target[k,4])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,3])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,1])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,2])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,4])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,1] > target[k,4]) && (target[k,4] > target[k,3]) && (target[k,3] > target[k,2])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,4])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,3])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,2])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,1] > target[k,4]) && (target[k,4] > target[k,2]) && (target[k,2] > target[k,3])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,4])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,2])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,3])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,1] > target[k,2]) && (target[k,2] > target[k,4]) && (target[k,4] > target[k,3])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,2])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,4])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,3])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else if((target[k,1] > target[k,3]) && (target[k,3] > target[k,2]) && (target[k,2] > target[k,4])){
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,3])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,2])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,4])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
        else{
             TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
			       TEMP1 <- t(TEMP)
			       TEMP2 <- TEMP1[1,-(target[k,3])]
			       TEMP2 <- t(TEMP2)
			       TEMP3 <- TEMP2[1,-(target[k,4])]
			       TEMP3 <- t(TEMP3)
			       TEMP4 <- TEMP3[1,-(target[k,2])]
			       TEMP4 <- t(TEMP4)
			       SIGNDIFF[k] <- updown*sign(MU_SAMPLED[k,target[k,4]] - (cutval * s_d[k] + max(TEMP4)[1]))
        }
  	}
		
    indic <- (SIGNDIFF) == 1
    indicsum<-indicsum+indic
    df<-nrarrays-2
    nscale<-rowSums((MU_SAMPLED2-DATA)^2)
    chisqvals<-2*rgamma(nrprobesets, df/2, rate = 1, scale = 1)
    sigma2_sampled<-nscale*chisqvals^(-1)
    SIGMA2_SAMPLED<-sigma2_sampled%*%t(rep(1,nrconditions))
  }
    f1<-indicsum/nrruns
    f2<-1-f1
    BFN4 <- (( ( (nrconditions^3) - (2*(nrconditions^2)) + (3*nrconditions) + 6 )*(nrconditions-4) ) * f1 / (24*f2))
}
