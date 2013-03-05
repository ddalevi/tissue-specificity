#####Caculate SPM from TiSGeD database###############################
#Please source this file before using ROKUspm.r######################

spm <- function(x){ #expression level for one gene

    le<-length(x)   
    sum<-0
    spm<-array(0,dim=c(1,le))
    for(i in 1:le){
       sum<-sum+x[i]^2
    } 
    for(j in 1:le){
        spm[j]<-x[j]^2/sum
    }
    return (spm)
}


############Return the largest SPM value for each probesets#######

SPMsummary <- function(DATA){

    spmall <- rep(0,nrow(DATA))
    for (i in 1:nrow(DATA)){
          spmall[i] <- max(spm(DATA[i,]))
         #spmall[i] <- sort(spm(DATA[i,]),decreasing=T)[2]
    }
    plot(sort(spmall))
    SPMsummary <- c(mean(spmall),quantile(spmall,prob=c(0,.2,.5,.8,1)))
    return (SPMsummary)
}
