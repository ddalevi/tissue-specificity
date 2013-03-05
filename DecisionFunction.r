#summary for DECISION function, input whole dataset

DECISION<-function(EXPRS,NUM,WEIGHT,GT=0.18,STL=-5,STH=-20){   #only apply to over expression
                               #EXPRS expression level, vector
                               #NUM -st largest to be concerted
  n<-length(EXPRS)+1-NUM     #exclude the larger values
  intensity<-EXPRS/max(EXPRS) #?median correction
  sortindex<-sort(intensity,index.return=TRUE)$ix[1:n]
  sorted<-sort(intensity,index.return=TRUE)$x[1:n]  #x in the paper


#Step 3, Discodancy statistical test  
  GAP <- sorted[n]-sorted[n-1]
  TAU <- GAP/(max(intensity)-sorted[1])
  LSP <- (n-2)* log10(1-TAU)

#Step 4, Adjust due to intensity baseline
  Q<-WEIGHT[sortindex[1:n]]
  baseline<-sum((Q*sorted)[1:(n-1)])/sum(Q[1:(n-1)])
  c<-mean(sorted)  
  b<-0.8
  lambda<-(1+(baseline/c)^b)^(-1)
  LSPad<-(n-2)*log(1-lambda*TAU)
  
##Step 5, Minimum intensity GAP criterion 
  #GT<-0.2
  if(GAP<GT){
    g<-0
  }
  else{
   g<-(GAP-GT)/(1-GT)
  }
  #STH <- -5 
  #STL <- -20
  if(LSPad > STH)  {s<-0}
  else if(LSPad<STL){s<-1}
  else s<-(LSPad-STH)/(STL-STH)

#Step 6, Overall Confidence in selective expression determination
  DECISION<-1-((1-s)^1.5*(1-g)^1.5*((0.3*(1-g)+0.7*(1-s))/((1-g)+(1-s)))^1.5)^(2/9)
  return(DECISION)
  #return(c(GAP,g,LSPad,s,DECISION))
  
}##

Dtest<-function(DATA,THRESHOLD,OUTNUM,GT=0.18,STL=-5,STH=-20){

  Dmatrix<-matrix(-1,nrow=nrow(DATA),ncol=OUTNUM)
  tnr<-array(0,nrow(DATA))
  
  for(i in 1:nrow(DATA)){
     for(j in 1:OUTNUM){
         Dmatrix[i,j]<-DECISION(DATA[i,],j,rep(1,ncol(DATA)),GT,STL,STH)
         if(Dmatrix[i,j]>=THRESHOLD) {tnr[i]<-j}  
     }
    #Dmatrix[i,] <- sort(Dmatrix[i,])
  }#for
  
  return(tnr)
}
