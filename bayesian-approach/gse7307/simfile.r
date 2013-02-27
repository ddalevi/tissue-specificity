###########################################################################################################
#READING THE ANNOTATION GPL570.annot FILE FOR HUMAN GSE7307 DATASET TO GET THE PROBESETS WITH FULL ANNOTATION.
#PROBESETS WITH NO REFERENCE TO GENE SYMBOLS ARE NOT CONSIDERED.
DATA <- strsplit(readLines("GPL570.annot",n=-1,ok=TRUE,warn=TRUE),"\n")
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
rm(ndata)
h_gse_data <- split(tdata,1:nrow(tdata))
rm(tdata)
#PRINT THE FINAL PROBESETS TO A FILE
sink("probe_number_human_gse7307.txt")
for(i in 1:length(h_gse_data)){
      for(j in 2:length(h_gse_data[[i]])){
          if(h_gse_data[[i]][j] != "0"){
              cat(h_gse_data[[i]][j],'\n')
          }
      }
}
sink()
###########################################################################################################

###########################################################################################################
#READ THE GENE EXPRESSION DATA
data<-read.table("gse7307_grouped.txt")
source("../commonfiles/BFN4.r") # EDIT THE DIRECTORY PATH
source("../commonfiles/BFN3.r") # EDIT THE DIRECTORY PATH
source("../commonfiles/BFN.r") # EDIT THE DIRECTORY PATH
source("../commonfiles/BF.r") # EDIT THE DIRECTORY PATH
#READ THE PROBESETS THAT ARE RELEVANT
probe_number = as.matrix(read.table("probe_number_human_gse7307.txt"))
###########################################################################################################
#GET THE GENE EXPRESSION DATA FOR THE RELEVANT PROBESETS
DATA = matrix(rep(0,length(probe_number)*ncol(data)),nrow=length(probe_number),ncol=ncol(data))
row.names(DATA) <- probe_number
for(i in 1:length(probe_number)){
	for(j in 1:length(row.names(data))){
		if(probe_number[i] == row.names(data)[j]){
			temp<- as.matrix(data[j,])
			DATA[i,] <- t(temp)
		}
	}
}
###########################################################################################################
#PERFORM MODIFIED BAYES FACTOR TO FIND TISSUE SPECIFIC, 2-SELECTIVE, 3-SELECTIVE AND 4-SELECTIVE PROBESETS
#IN THE GENE EXPRESSION DATA. NRREPLIC IS THE NUMBER OF REPLICATES FOR TISSUES IN THE DATA.
nrreplic = c(3,4,5,2,3,1,4,6,3,23,4,3,1,1,4,3,4,3,6,22,8,5,1,1,6,5,7,1,9,5,1,6,2,5,4,3,4,5,1,4,3,4)
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
bfn4 <- BFN4(DATA,nrreplic,1,nrruns,1.79)
bfn3 <- BFN3(DATA,nrreplic,1,nrruns,1.79)
bfn1 <- BFN(DATA,nrreplic,1,nrruns,1.64)
bf1 <- BF(DATA,nrreplic,1,nrruns,1.79)
#ASSIGN PROBESETS TO THE APPROPRIATE TISSUE SELECTIVE CATEGORY BASED ON THE THRESHOLD
for(i in 1:dat_row){
      if(bfn4[i] > 261168){
              tissue_specific[i,max_tissue[i,1]] = 4
              tissue_specific[i,max_tissue[i,2]] = 4
              tissue_specific[i,max_tissue[i,3]] = 4
              tissue_specific[i,max_tissue[i,4]] = 4
      }
      if(bfn3[i] > 26784 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 3
              tissue_specific[i,max_tissue[i,2]] = 3
              tissue_specific[i,max_tissue[i,3]] = 3
      }
      if(bfn1[i] > 2007 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 2
              tissue_specific[i,max_tissue[i,2]] = 2
      }
      if(bf1[i] > 38 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 1
      }
}
#PRINT THE RESULTS IN TO A FILE FOR FURTHER PROCESSING
k <- 1
sink("human_genes_tissue_specificity_gse7307_final.txt")
for(i in 1:length(h_gse_data)){
      for(j in 2:length(h_gse_data[[i]])){
          if(h_gse_data[[i]][j] != "0"){
              cat(h_gse_data[[i]][1],'\t',h_gse_data[[i]][j],'\t',tissue_specific[k,],'\n')
              k <- k+1
          }
      }
}
sink()
