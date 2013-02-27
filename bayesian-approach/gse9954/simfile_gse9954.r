###########################################################################################################
#READING THE ANNOTATION GPL1261.annot FILE FOR MOUSE GSE9954 DATASET TO GET THE PROBESETS WITH FULL ANNOTATION.
#PROBESETS WITH NO REFERENCE TO GENE SYMBOLS ARE NOT CONSIDERED.
DATA <- strsplit(readLines("GPL1261.annot",n=-1,ok=TRUE,warn=TRUE),"\n")
ndata <- matrix(rep(0,length(DATA)*2),nrow=length(DATA),ncol=2)
for(i in 1:length(DATA)){
     if(regexpr("^([0-9A])",DATA[[i]],perl=TRUE) == TRUE){
          if(strsplit(DATA[[i]],"\t")[[1]][3] == ""){
              ndata[i,2] <- as.matrix(strsplit(DATA[[i]],"\t")[[1]][1])
              ndata[i,1] <- Inf
          }
          else{
              ndata[i,2] <- as.matrix(strsplit(DATA[[i]],"\t")[[1]][1])
              ndata[i,1] <- as.matrix(strsplit(DATA[[i]],"\t")[[1]][3])
          }
     }
}
ndata <- ndata[-which(ndata == "Inf", arr.ind=TRUE)[,1],]
ndata <- ndata[-which(ndata == "0", arr.ind=TRUE)[,1],]
rm(DATA)
#COMBINING ALL PROBESETS OF A GENE
ldata <- split(ndata,1:nrow(ndata))
for(i in 1:length(ldata)){
      count <- 3
      for(j in i:length(ldata)){
          if((ldata[[i]][1] == ldata[[j]][1]) && (ldata[[i]][1] != "NULL") && (i!=j)){
              ldata[[i]][count] <- ldata[[j]][2]
              ldata[[j]][1] <- "NULL"
              count <- count + 1
          }
      }
}
cmax <- 2
for(i in 1:length(ldata)){
      if(cmax < length(ldata[[i]])){
          cmax <- length(ldata[[i]])
      }
}
tdata <- matrix(rep(0,length(ldata)*cmax),nrow=length(ldata),ncol=cmax)
for(i in 1:length(ldata)){
      if(strsplit(ldata[[i]]," ")[[1]][1] != "NULL"){
          for(j in 1:length(ldata[[i]])){
              tdata[i,j] <- as.matrix(strsplit(ldata[[i]]," ")[[j]][1])
          }
      }
      else{
          tdata[i,1] <- Inf
      }
}
tdata <- tdata[-(which(tdata=="Inf",arr.ind=TRUE)[,1]),]
rm(ldata)
rm(cmax)
rm(ndata)
m_data <- split(tdata,1:nrow(tdata))
rm(tdata)
#PRINT THE FINAL PROBESETS TO A FILE
sink("probe_number_mouse.txt")
for(i in 1:length(m_data)){
      for(j in 2:length(m_data[[i]])){
          if(m_data[[i]][j] != "0"){
              cat(m_data[[i]][j],'\n')
          }
      }
}
sink()
###########################################################################################################

###########################################################################################################
#READ THE GENE EXPRESSION DATA
data = read.table("GSE9954_mini.txt")
source("../commonfiles/BFN4.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/BFN3.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/BFN.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/BF.r") #EDIT THE DIRECTORY PATH
#REMOVE THE CELL LINE TISSUES FROM GENE EXPRESSION DATA
data$GSM252122<-NULL
data$GSM252123<-NULL
data$GSM252124<-NULL
data$GSM252131<-NULL
data$GSM252132<-NULL
data$GSM252133<-NULL
###########################################################################################################
#PERFORM MODIFIED BAYES FACTOR TO FIND TISSUE SPECIFIC, 2-SELECTIVE, 3-SELECTIVE AND 4-SELECTIVE PROBESETS
#IN THE GENE EXPRESSION DATA. NRREPLIC IS THE NUMBER OF REPLICATES FOR TISSUES IN THE DATA.
nrreplic = c(3,3,4,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3)
ndata <- matrix(rep(0,nrow(data)*sum(nrreplic)),nrow=nrow(data),ncol=sum(nrreplic))
for(i in 1:nrow(data)){
      m1<-mean(c(data$GSM252077[i],data$GSM252078[i],data$GSM252079[i]))
      m2<-mean(c(data$GSM252096[i],data$GSM252097[i],data$GSM252098[i],data$GSM252099[i],data$GSM252100[i]))
      if(m1>m2){
                ndata[i,1:32] <- as.matrix(data[i,1:32])
                ndata[i,33:59] <- as.matrix(data[i,38:64])
      }
      else{
                ndata[i,1:13] <- as.matrix(data[i,1:13])
                temp <- sort(c(data$GSM252096[i],data$GSM252097[i],data$GSM252098[i],data$GSM252099[i],data$GSM252100[i]),decreasing=TRUE)
                ndata[i,14] <- as.matrix(temp[1])
                ndata[i,15] <- as.matrix(temp[2])
                ndata[i,16] <- as.matrix(temp[5])
                ndata[i,17:32] <- as.matrix(data[i,17:32])
                ndata[i,33:59] <- as.matrix(data[i,38:64])
      }
}
#READ THE PROBESETS THAT ARE RELEVANT
probe_number = as.matrix(read.table("probe_number_mouse.txt",comment.char=""))
###########################################################################################################
#GET THE GENE EXPRESSION DATA FOR THE RELEVANT PROBESETS
DATA = matrix(rep(0,length(probe_number)*ncol(ndata)),nrow=length(probe_number),ncol=ncol(ndata))
row.names(DATA) <- probe_number
for(i in 1:length(probe_number)){
	for(j in 1:length(row.names(data))){
		if(probe_number[i] == row.names(data)[j]){
			DATA[i,] <- ndata[j,]
		}
	}
}
###########################################################################################################
DATA = log2(DATA)
nrruns = 10000
dat_row <- nrow(DATA)
dat_col <- ncol(DATA)
#GET THE TOP 4 HIGHLY EXPRESSED TISSUES FOR ALL PROBESETS
tissue_specific <- matrix(rep(0,dat_row*length(nrreplic)), nrow = dat_row, ncol = length(nrreplic))
row.names(tissue_specific) <- probe_number
max_tissue <- matrix(rep(0,dat_row*4),nrow=dat_row,ncol=4)
lower <- 1
AV1<-matrix(c(rep(0,nrow(DATA)*ncol(DATA))),nrow=nrow(DATA),ncol=ncol(DATA))
AV2<-matrix(c(rep(0,nrow(DATA)*length(nrreplic))),nrow=nrow(DATA),ncol=length(nrreplic))
for (i in (1:length(nrreplic))){
    		upper<-lower+nrreplic[i]-1
    		v<-rep(1,nrreplic[i])
    		average<-(DATA[,lower:upper])%*%(v*(t(v)%*%v)^(-1))
    		AV1[,lower:upper]<-average%*%t(v)
    		AV2[,i]<-average
    		lower<-upper+1
}
for(i in 1:dat_row){
      max_tissue[i,1] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[1],arr.ind=TRUE)[1]
      max_tissue[i,2] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[2],arr.ind=TRUE)[1]
      max_tissue[i,3] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[3],arr.ind=TRUE)[1]
      max_tissue[i,4] <- which(AV2[i,] == sort(AV2[i,],decreasing=TRUE)[4],arr.ind=TRUE)[1]
}
#RUN THE MODIFIED BAYES FACTOR. NRRUNS IS NUMBER OF TIMES TO RUN BAYES FACTOR FOR EACH PROBESET.
#"1" IN THE FUNCTION IS FOR UP REGULATED (HIGHLY EXPRESSED) GENES. THE LAST PARAMETER IN THE
#BAYES FACTOR IS FOR VARIANCE BETWEEN TISSUES.
bfn4 <- BFN4(DATA,nrreplic,1,nrruns,3.5)
bfn3 <- BFN3(DATA,nrreplic,1,nrruns,3.5)
bfn1 <- BFN(DATA,nrreplic,1,nrruns,3.5)
bf1 <- BF(DATA,nrreplic,1,nrruns,1.79)
#ASSIGN PROBESETS TO THE APPROPRIATE TISSUE SELECTIVE CATEGORY BASED ON THE THRESHOLD
for(i in 1:dat_row){
      if(bfn4[i] > 15500){
              tissue_specific[i,max_tissue[i,1]] = 4
              tissue_specific[i,max_tissue[i,2]] = 4
              tissue_specific[i,max_tissue[i,3]] = 4
              tissue_specific[i,max_tissue[i,4]] = 4
      }
      if(bfn3[i] > 3872 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 3
              tissue_specific[i,max_tissue[i,2]] = 3
              tissue_specific[i,max_tissue[i,3]] = 3
      }
      if(bfn1[i] > 680 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 2
              tissue_specific[i,max_tissue[i,2]] = 2
      }
      if(bf1[i] > 32 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
             tissue_specific[i,max_tissue[i,1]] = 1
      }
}
#PRINT THE RESULTS IN TO A FILE FOR FURTHER PROCESSING
k <- 1
sink("mouse_genes_tissue_specificity_final.txt")
for(i in 1:length(m_data)){
      for(j in 2:length(m_data[[i]])){
          if(m_data[[i]][j] != "0"){
              cat(m_data[[i]][1],'\t',m_data[[i]][j],'\t',tissue_specific[k,],'\n')
              k <- k+1
          }
      }
}
sink()
