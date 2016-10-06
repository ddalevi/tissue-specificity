#####combine ROKU and SPM method ############################
#############################################################

#Please source the "shannon.entropy" and "spm" function first. 
#Or use source(".../ROKUspm.r",chdir =TRUE)
   source("Shannon_Entropy.r")
   source("spm.r")
   
################ROKU#########################################
##log(n!)####################################################
stirling <- function(x)
{
  1/2*log(2*pi)+(1/2+x)*log(x)-x
}


##U statistic################################################

Ut<-function(x,s){  
n <- length(x)
#sigma <- sqrt(sum((x-mean(x))^2)/n)
ut <- n*log(sd(x))+sqrt(2)*stirling(n)/n*s
return(ut)
}


Us<-function(x,s)   #return the outlier Ut matrix
{
 n <-length(x)
 nmax <- min(length(x)/2,s)
       #default maximum outlier number is 5

 ordered <- sort(x)  #lose the original index
 t <- scale(ordered) #normalized 
 
 Us<-matrix(10,nrow=nmax+1,ncol=nmax+1)
 for (i in 0:nmax){
      for(j in 0:i){
          Us[i+1,j+1]<-Ut(t[(1+j):(n+j-i)],i)
      }
 }
 return (Us)
}

####outlier detection##########################################

outlier<-function(x,s){ #s outlier number 
 ordered <- sort(x,index.return=TRUE)$x
 orderseries <- sort(x,index.return=TRUE)$ix
 Um<-Us(x,s)   
 pos<-which(Um==min(Um),arr.ind=TRUE) #the position of min U in Us matrix 
 #i=pos[1,1]-1 j=pos[1,2]-1 #assume only one Ut has the min value.
 nlow<-as.integer(pos[1,2]-1)           #number of lowest outlier
 nhigh<-as.integer(pos[1,1]-pos[1,2])   #number of highest outlier
 noutlier<-as.integer(pos[1,1]-1)       #i number of outlier
 count<-c(noutlier,nlow,nhigh) 
 return(count)
}

#####ROKUspm###################################################
#Optimized value for Shannon entropy and two spms (Tspm1 for 
#tissue-specific and Tspm2 for 2-selective) are needed.

ROKUspm <- function(DATA,Tentropy,Tspm1,Tspm2){

   nout <- array(0,nrow(DATA))
   for(i in 1:nrow(DATA)){
     entropy <- shannon.entropy(DATA[i,])
     spm <- spm(DATA[i,])
     deout <- outlier(DATA[i,],15)[3]
     spm1 <- max(spm)
     spm2 <- sort(spm, decreasing = TRUE)[2]
     
     if(entropy <= Tentropy & deout>=1 ){
        if (spm2>=Tspm2&deout>=2){
            nout[i]=2
            #tissuenames[i,1] <-tname[i,1] 
        }
        else if(spm1>=Tspm1){
            nout[i]=1
           #tissuenames[i,1:2] <-tname[i,1:2]
        }
#        else if(deout>=3){nout[i] <- deout}
     }#if entropy
   }#for
  return (nout)
}#function
