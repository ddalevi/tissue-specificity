#BAYES FACTOR FOR TISSUE SPECIFIC CASE
BF<-function(DATA,nrreplics,updown,nrruns,cutvalue)
{
	cutval <- cutvalue
	nrprobesets<-nrow(DATA)
	nrarrays<-ncol(DATA)
	nrconditions<-length(nrreplics)
	#IN TISSUE SPECIFIC CASE TOP MOST HIGHLY EXPRESSED TISSUE IS CONSIDERED FOR ALL PROBESETS
	target <- matrix(rep(0,nrprobesets*1),nrow=nrprobesets,ncol=1)
	SIGNDIFF <- matrix(rep(0,nrprobesets),nrow=nrprobesets)

	#Preliminary: Calculation of condition-specific averages for all probesets
	lower<-1
	AV1<-matrix(c(rep(0,nrprobesets*nrarrays)),nrow=nrprobesets,ncol=nrarrays)
	AV2<-matrix(c(rep(0,nrprobesets*nrconditions)),nrow=nrprobesets,ncol=nrconditions)
	AV3<-matrix(c(rep(0,nrprobesets*(nrconditions-1))),nrow=nrprobesets,ncol=(nrconditions-1))
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
        temp <- AV2[i,-(target[i,1])]
        AV3[i,] <- t(temp)
  }
  s_d <- sd(t(AV3),na.rm=FALSE)
	#2. Gibbs Sampler to calculate f1
	indicsum<-rep(0,nrprobesets)
	NRREPLICS<-rep(1,nrprobesets)%*%t(nrreplics)
	MU_SAMPLED2<-matrix(rep(0,nrprobesets*nrarrays),nrow=nrprobesets)
	for (i in (1:nrruns)){
    		SEM<-sqrt(SIGMA2_SAMPLED/NRREPLICS)
    		MU_SAMPLED<-AV2+(SEM*matrix(rnorm(nrprobesets*nrconditions),nrow=nrprobesets))
    		if(nrruns == 1){
        }
        else{
             for(k in 1:nrprobesets){
                   temp <- MU_SAMPLED[k,-(target[k,1])]
                   AV3[k,] <- t(temp)
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
            TEMP <- TEMP_MU_SAMPLED[k,-(target[k,1])]
            TEMP <- t(TEMP)
            SIGNDIFF[k] <-updown*sign(MU_SAMPLED[k,target[k,1]] - (cutval * s_d[k] + max(TEMP)[1]))
        }
    		indic <- (SIGNDIFF) == 1
    		indicsum<-indicsum+indic
    		df<-nrarrays-2
    		nscale<-rowSums((MU_SAMPLED2-DATA)^2)
    		chisqvals<-2*rgamma(nrprobesets, df/2, rate = 1, scale = 1)
   		  sigma2_sampled<-nscale*chisqvals^(-1)
    		SIGMA2_SAMPLED<-sigma2_sampled%*%t(rep(1,nrconditions))
	}
	#CALCULATION OF BAYES FACTOR
	f1<-indicsum/nrruns
	f2<-1-f1
	BF<-nrconditions*f1/((nrconditions/(nrconditions-1))*f2)
}
