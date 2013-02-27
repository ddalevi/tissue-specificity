###########################################################################################################
#READING THE ANNOTATION GPL85.annot FILE FOR RAT GSE952 DATASET TO GET THE PROBESETS WITH FULL ANNOTATION.
#PROBESETS WITH NO REFERENCE TO GENE SYMBOLS ARE NOT CONSIDERED.
DATA <- strsplit(readLines("GPL85.annot",n=-1,ok=TRUE,warn=TRUE),"\n")
ndata <- matrix(rep(0,length(DATA)*2),nrow=length(DATA),ncol=2)
for(i in 1:length(DATA)){
     if(regexpr("^([A-HJ-Zr])",DATA[[i]],perl=TRUE) == TRUE){
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
rm(ldata,cmax,ndata)
r_data <- split(tdata,1:nrow(tdata))
rm(tdata)
#PRINT THE PROBESETS TO A FILE
sink("probe_number_rat.txt")
for(i in 1:length(r_data)){
      for(j in 2:length(r_data[[i]])){
          if(r_data[[i]][j] != "0"){
              cat(r_data[[i]][j],'\n')
          }
      }
}
sink()
###########################################################################################################

###########################################################################################################
#READ THE GENE EXPRESSION DATA
data <- read.table("GSE952_full_series.txt")
source("../commonfiles/BFN4.r") # EDIT THE DIRECTORY PATH
source("../commonfiles/BFN3.r") # EDIT THE DIRECTORY PATH
source("../commonfiles/BFN.r") # EDIT THE DIRECTORY PATH
source("../commonfiles/BF.r") # EDIT THE DIRECTORY PATH
#REMOVE THE CELL LINE TISSUES FROM GENE EXPRESSION DATA
data$GSM15217 <- NULL
data$GSM15218 <- NULL
#READ THE PROBESETS THAT ARE RELEVANT
probe_number = as.matrix(read.table("probe_number_rat.txt",comment.char=""))
###########################################################################################################
#GET THE GENE EXPRESSION DATA FOR THE RELEVANT PROBESETS
DATA = matrix(rep(0,length(probe_number)*ncol(data)),nrow=length(probe_number))
rownames(DATA) <- probe_number
colnames(DATA) <- colnames(data)
for(i in 1:length(probe_number)){
	for(j in 1:length(row.names(data))){
		if(probe_number[i] == row.names(data)[j]){
			temp<- as.matrix(data[j,])
			DATA[i,] <- t(temp)
		}
	}
}
###########################################################################################################
rm(data)
#LOG2 TRANFORMATION OF DATA BEFORE PERFORMING BAYES FACTOR AND REMOVING PROBESETS WITH -Inf VALUES
data <- log2(DATA)
temp <- row.names(data)[which(data == -Inf,arr.ind=TRUE)[,1]]
data <- data[-which(data == -Inf, arr.ind=TRUE)[,1],]
data <- as.data.frame(data)

for(i in 1:length(r_data)){
      for(j in 2:length(r_data[[i]])){
          for(k in 1:length(temp)){
              if(r_data[[i]][j] == temp[k]){
                   r_data[[i]][j] <- "0"
              }
          }
      }
}
#PRINT THE FINAL PROBESETS TO A FILE
sink("probe_number_rat_final.txt")
for(i in 1:length(r_data)){
      for(j in 2:length(r_data[[i]])){
          if(r_data[[i]][j] != "0"){
              cat(r_data[[i]][j],'\n')
          }
      }
}
sink()
#PERFORM MODIFIED BAYES FACTOR TO FIND TISSUE SPECIFIC, 2-SELECTIVE, 3-SELECTIVE AND 4-SELECTIVE PROBESETS
#IN THE GENE EXPRESSION DATA. NRREPLIC IS THE NUMBER OF REPLICATES FOR TISSUES IN THE DATA.
nrreplic <- c(2,2,6,2,2,2,2,2,2,2,2,2)
ndata <- matrix(rep(0,nrow(data)*sum(nrreplic)),nrow = nrow(data), ncol = sum(nrreplic))
row.names(ndata) <- row.names(data)

#Amygdala and Amygdala,central nucleus are in one group and dorsal root ganglion fischer is removed
for(i in 1:nrow(data)){
      m1 <- mean(c(as.matrix(data$GSM15169[i]),as.matrix(data$GSM15170[i]),as.matrix(data$GSM15171[i]),as.matrix(data$GSM15172[i]),as.matrix(data$GSM15173[i]),as.matrix(data$GSM15174[i])))
      s1 <- sd(c(as.matrix(data$GSM15169[i]),as.matrix(data$GSM15170[i]),as.matrix(data$GSM15171[i]),as.matrix(data$GSM15172[i]),as.matrix(data$GSM15173[i]),as.matrix(data$GSM15174[i])))
      m2 <- mean(c(as.matrix(data$GSM15125[i]),as.matrix(data$GSM15221[i]),as.matrix(data$GSM15222[i])))
      s2 <- sd(c(as.matrix(data$GSM15125[i]),as.matrix(data$GSM15221[i]),as.matrix(data$GSM15222[i])))
      m3 <- mean(as.matrix(data$GSM15179[i]))
      m4 <- mean(c(as.matrix(data$GSM15229[i]),as.matrix(data$GSM15230[i])))
      s4 <- sd(c(as.matrix(data$GSM15229[i]),as.matrix(data$GSM15230[i])))
      m5 <- mean(c(as.matrix(data$GSM15183[i]),as.matrix(data$GSM15184[i]),as.matrix(data$GSM15185[i])))
      s5 <- sd(c(as.matrix(data$GSM15183[i]),as.matrix(data$GSM15184[i]),as.matrix(data$GSM15185[i])))
      m6 <- mean(c(as.matrix(data$GSM15223[i]),as.matrix(data$GSM15224[i])))
      s6 <- sd(c(as.matrix(data$GSM15223[i]),as.matrix(data$GSM15224[i])))
      m7 <- mean(c(as.matrix(data$GSM15195[i]),as.matrix(data$GSM15196[i])))
      s7 <- sd(c(as.matrix(data$GSM15195[i]),as.matrix(data$GSM15196[i])))
      m8 <- mean(c(as.matrix(data$GSM15199[i]),as.matrix(data$GSM15200[i])))
      s8 <- sd(c(as.matrix(data$GSM15199[i]),as.matrix(data$GSM15200[i])))
      m9 <- mean(c(as.matrix(data$GSM15225[i]),as.matrix(data$GSM15226[i])))
      s9 <- sd(c(as.matrix(data$GSM15225[i]),as.matrix(data$GSM15226[i])))
      m10 <- mean(c(as.matrix(data$GSM15227[i]),as.matrix(data$GSM15228[i])))
      s10 <- sd(c(as.matrix(data$GSM15227[i]),as.matrix(data$GSM15228[i])))
      m11 <- mean(c(as.matrix(data$GSM15219[i]),as.matrix(data$GSM15220[i])))
      s11 <- sd(c(as.matrix(data$GSM15219[i]),as.matrix(data$GSM15220[i])))
      m12 <- mean(c(as.matrix(data$GSM15231[i]),as.matrix(data$GSM15232[i])))
      s12 <- sd(c(as.matrix(data$GSM15231[i]),as.matrix(data$GSM15232[i])))
      m13 <- mean(c(as.matrix(data$GSM15233[i]),as.matrix(data$GSM15234[i])))
      s13 <- sd(c(as.matrix(data$GSM15233[i]),as.matrix(data$GSM15234[i])))
      m14 <- mean(c(as.matrix(data$GSM15151[i]),as.matrix(data$GSM15152[i]),as.matrix(data$GSM15153[i]),as.matrix(data$GSM15154[i]),as.matrix(data$GSM15155[i]),as.matrix(data$GSM15156[i])))
      s14 <- sd(c(as.matrix(data$GSM15151[i]),as.matrix(data$GSM15152[i]),as.matrix(data$GSM15153[i]),as.matrix(data$GSM15154[i]),as.matrix(data$GSM15155[i]),as.matrix(data$GSM15156[i])))
      m15 <- mean(c(as.matrix(data$GSM15235[i]),as.matrix(data$GSM15236[i])))
      s15 <- sd(c(as.matrix(data$GSM15235[i]),as.matrix(data$GSM15236[i])))
      m <- which(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15) == max(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)),arr.ind=TRUE)
      s3 <- mean(c(s1,s2,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15))
      switch(m[1],{ndata[i,1] <- as.matrix(sort(c(data$GSM15169[i],data$GSM15170[i],data$GSM15171[i],data$GSM15172[i],data$GSM15173[i],data$GSM15174[i]),decreasing=TRUE)[1]);ndata[i,2] <- as.matrix(sort(c(data$GSM15169[i],data$GSM15170[i],data$GSM15171[i],data$GSM15172[i],data$GSM15173[i],data$GSM15174[i]),decreasing=TRUE)[6])},
			            {ndata[i,1] <- as.matrix(sort(c(data$GSM15125[i],data$GSM15221[i],data$GSM15222[i]),decreasing=TRUE)[1]);ndata[i,2] <- as.matrix(sort(c(data$GSM15125[i],data$GSM15221[i],data$GSM15222[i]),decreasing=TRUE)[3])},
			            {ndata[i,1] <- as.matrix(data$GSM15179[i]) + s3;ndata[i,2] <- as.matrix(data$GSM15179[i]) - s3},
			            {ndata[i,1] <- as.matrix(data$GSM15229[i]);ndata[i,2] <- as.matrix(data$GSM15230[i])},
			            {ndata[i,1] <- as.matrix(sort(c(data$GSM15183[i],data$GSM15184[i],data$GSM15185[i]),decreasing=TRUE)[1]);ndata[i,2] <- as.matrix(sort(c(data$GSM15183[i],data$GSM15184[i],data$GSM15185[i]),decreasing=TRUE)[3])},
			            {ndata[i,1] <- as.matrix(data$GSM15223[i]);ndata[i,2] <- as.matrix(data$GSM15224[i])},
			            {ndata[i,1] <- as.matrix(data$GSM15195[i]);ndata[i,2] <- as.matrix(data$GSM15196[i])},
                  {ndata[i,1] <- as.matrix(data$GSM15199[i]);ndata[i,2] <- as.matrix(data$GSM15200[i])},
			            {ndata[i,1] <- as.matrix(data$GSM15225[i]);ndata[i,2] <- as.matrix(data$GSM15226[i])},
			            {ndata[i,1] <- as.matrix(data$GSM15227[i]);ndata[i,2] <- as.matrix(data$GSM15228[i])},
			            {ndata[i,1] <- as.matrix(data$GSM15219[i]);ndata[i,2] <- as.matrix(data$GSM15220[i])},
			            {ndata[i,1] <- as.matrix(data$GSM15231[i]);ndata[i,2] <- as.matrix(data$GSM15232[i])},
			            {ndata[i,1] <- as.matrix(data$GSM15233[i]);ndata[i,2] <- as.matrix(data$GSM15234[i])},
			            {ndata[i,1] <- as.matrix(sort(c(data$GSM15151[i],data$GSM15152[i],data$GSM15153[i],data$GSM15154[i],data$GSM15155[i],data$GSM15156[i]),decreasing=TRUE)[1]);ndata[i,2] <- as.matrix(sort(c(data$GSM15151[i],data$GSM15152[i],data$GSM15153[i],data$GSM15154[i],data$GSM15155[i],data$GSM15156[i]),decreasing=TRUE)[6])},
			            {ndata[i,1] <- as.matrix(data$GSM15235[i]);ndata[i,2] <- as.matrix(data$GSM15236[i])})
      ndata[i,3]	<- as.matrix(data$GSM15215[i]); ndata[i,4]	<- as.matrix(data$GSM15216[i]); ndata[i,5]	<- as.matrix(data$GSM15138[i])
      ndata[i,6]	<- as.matrix(data$GSM15139[i]); ndata[i,7]	<- as.matrix(data$GSM15140[i]); ndata[i,8]	<- as.matrix(data$GSM15141[i])
      ndata[i,9]	<- as.matrix(data$GSM15142[i]); ndata[i,10]	<- as.matrix(data$GSM15143[i]); ndata[i,11]	<- as.matrix(data$GSM15237[i])
      ndata[i,12]	<- as.matrix(data$GSM15238[i]); ndata[i,13]	<- as.matrix(data$GSM15193[i]); ndata[i,14]	<- as.matrix(data$GSM15194[i])
      ndata[i,15]	<- as.matrix(data$GSM15203[i]); ndata[i,16]	<- as.matrix(data$GSM15204[i]); ndata[i,17]	<- as.matrix(data$GSM15205[i])
      ndata[i,18]	<- as.matrix(data$GSM15206[i]); ndata[i,19]	<- as.matrix(data$GSM15207[i]); ndata[i,20]	<- as.matrix(data$GSM15208[i])
      ndata[i,21]	<- as.matrix(data$GSM15209[i]); ndata[i,22]	<- as.matrix(data$GSM15210[i]); ndata[i,23]	<- as.matrix(data$GSM15211[i])
      ndata[i,24]	<- as.matrix(data$GSM15212[i]); ndata[i,25]	<- as.matrix(data$GSM15213[i]); ndata[i,26]	<- as.matrix(data$GSM15214[i])
      ndata[i,27]	<- as.matrix(data$GSM15201[i]); ndata[i,28]	<- as.matrix(data$GSM15202[i]);
      
}

rm(DATA)
DATA <- ndata
nrruns = 10000
dat_row <- nrow(DATA)
dat_col <- ncol(DATA)
#GET THE TOP 4 HIGHLY EXPRESSED TISSUES FOR ALL PROBESETS
tissue_specific <- matrix(rep(0,dat_row*length(nrreplic)), nrow = dat_row, ncol = length(nrreplic))
row.names(tissue_specific) <- row.names(data)
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
bfn4 <- BFN4(DATA,nrreplic,1,nrruns,2)
bfn3 <- BFN3(DATA,nrreplic,1,nrruns,2)
bfn1 <- BFN(DATA,nrreplic,1,nrruns,2)
bf1 <- BF(DATA,nrreplic,1,nrruns,1.6)
#ASSIGN PROBESETS TO THE APPROPRIATE TISSUE SELECTIVE CATEGORY BASED ON THE THRESHOLD
for(i in 1:dat_row){
      if(bfn4[i] > 1153){
              tissue_specific[i,max_tissue[i,1]] = 4
              tissue_specific[i,max_tissue[i,2]] = 4
              tissue_specific[i,max_tissue[i,3]] = 4
              tissue_specific[i,max_tissue[i,4]] = 4
      }
      if(bfn3[i] > 657 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 3
              tissue_specific[i,max_tissue[i,2]] = 3
              tissue_specific[i,max_tissue[i,3]] = 3
      }
      if(bfn1[i] > 152 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 2
              tissue_specific[i,max_tissue[i,2]] = 2
      }
      if(bf1[i] > 24 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
             tissue_specific[i,max_tissue[i,1]] = 1
      }
}
#PRINT THE RESULTS IN TO A FILE FOR FURTHER PROCESSING
k <- 1
sink("rat_genes_tissue_specificity_final.txt")
for(i in 1:length(r_data)){
      for(j in 2:length(r_data[[i]])){
          if(r_data[[i]][j] != "0"){
              cat(r_data[[i]][1],'\t',r_data[[i]][j],'\t',tissue_specific[k,],'\n')
              k <- k+1
          }
      }
}
sink()




