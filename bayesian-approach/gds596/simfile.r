###########################################################################################################
#READING THE ANNOTATION GPL96.annot FILE FOR HUMAN GDS596 DATASET TO GET THE PROBESETS WITH FULL ANNOTATION.
#PROBESETS WITH NO REFERENCE TO GENE SYMBOLS ARE NOT CONSIDERED.
DATA <- strsplit(readLines("GPL96.annot",n=-1,ok=TRUE,warn=TRUE),"\n")
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
rm(ldata,cmax,ndata)
h_gds_data <- split(tdata,1:nrow(tdata))
rm(tdata)
#PRINT THE FINAL PROBESETS TO A FILE
sink("probe_number_human_gds596.txt")
for(i in 1:length(h_gds_data)){
      for(j in 2:length(h_gds_data[[i]])){
          if(h_gds_data[[i]][j] != "0"){
              cat(h_gds_data[[i]][j],'\n')
          }
      }
}
sink()
###########################################################################################################

###########################################################################################################
#READ THE GENE EXPRESSION DATA
data = read.table("gds596_original.txt")
source("../commonfiles/BFN4.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/BFN3.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/BFN.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/BF.r") #EDIT THE DIRECTORY PATH
#REMOVE THE CELL LINE TISSUES FROM GENE EXPRESSION DATA
data$GSM18883 <- NULL
data$GSM18884 <- NULL
data$GSM18871 <- NULL
data$GSM18872 <- NULL
data$GSM18881 <- NULL
data$GSM18882 <- NULL
data$GSM18869 <- NULL
data$GSM18870 <- NULL
data$GSM18885 <- NULL
data$GSM18886 <- NULL
data$GSM18877 <- NULL
data$GSM18878 <- NULL
data$GSM18875 <- NULL
data$GSM18876 <- NULL
data$GSM18907 <- NULL
data$GSM18908 <- NULL
data$GSM18879 <- NULL
data$GSM18880 <- NULL
data$GSM18973 <- NULL
data$GSM18974 <- NULL
data$GSM18873 <- NULL
data$GSM18874 <- NULL
data$GSM18963 <- NULL
data$GSM18964 <- NULL
data$GSM18945 <- NULL
data$GSM18946 <- NULL
data$GSM18905 <- NULL
data$GSM18906 <- NULL
data$GSM18965 <- NULL
data$GSM18966 <- NULL
ndata <- matrix(rep(0,nrow(data)*60),nrow=nrow(data),ncol=60)
row.names(ndata) <- row.names(data)
for(i in 1:nrow(data)){
      m1<-mean(c(data$GSM18981[i],data$GSM18982[i]))
      m2<-mean(c(data$GSM18983[i],data$GSM18984[i]))
      m3<-mean(c(data$GSM18985[i],data$GSM18986[i]))
      m4<-mean(c(data$GSM18987[i],data$GSM18988[i]))
      m5<-mean(c(data$GSM18989[i],data$GSM18990[i]))
      m<-which(c(m1,m2,m3,m4,m5) == max(c(m1,m2,m3,m4,m5)),arr.ind=TRUE)
      switch(m[1],{ndata[i,1] <- as.matrix(data$GSM18981[i]);ndata[i,2] <- as.matrix(data$GSM18982[i])},{ndata[i,1] <- as.matrix(data$GSM18983[i]);ndata[i,2] <- as.matrix(data$GSM18984[i])},{ndata[i,1] <- as.matrix(data$GSM18985[i]);ndata[i,2] <- as.matrix(data$GSM18986[i])},{ndata[i,1] <- as.matrix(data$GSM18987[i]);ndata[i,2] <- as.matrix(data$GSM18988[i])},{ndata[i,1] <- as.matrix(data$GSM18989[i]);ndata[i,2] <- as.matrix(data$GSM18990[i])})
      m1 <- mean(c(data$GSM18959[i],data$GSM18960[i]))
      m2 <- mean(c(data$GSM19015[i],data$GSM19016[i]))
      m<-which(c(m1,m2) == max(c(m1,m2)),arr.ind=TRUE)
      switch(m[1],{ndata[i,3] <- as.matrix(data$GSM18959[i]);ndata[i,4] <- as.matrix(data$GSM18960[i])},{ndata[i,3] <- as.matrix(data$GSM19015[i]);ndata[i,4] <- as.matrix(data$GSM19016[i])})
      m1 <- mean(c(data$GSM18977[i],data$GSM18978[i]))
      m2 <- mean(c(data$GSM18979[i],data$GSM18980[i]))
      m<-which(c(m1,m2) == max(c(m1,m2)),arr.ind=TRUE)
      switch(m[1],{ndata[i,5] <- as.matrix(data$GSM18977[i]);ndata[i,6] <- as.matrix(data$GSM18978[i])},{ndata[i,5] <- as.matrix(data$GSM18979[i]);ndata[i,6] <- as.matrix(data$GSM18980[i])})
      m1<-mean(c(data$GSM19013[i],data$GSM19014[i]))
      m2<-mean(c(data$GSM18971[i],data$GSM18972[i]))
      m3<-mean(c(data$GSM18969[i],data$GSM18970[i]))
      m<-which(c(m1,m2,m3) == max(c(m1,m2,m3)),arr.ind=TRUE)
      switch(m[1],{ndata[i,7] <- as.matrix(data$GSM19013[i]);ndata[i,8] <- as.matrix(data$GSM19014[i])},{ndata[i,7] <- as.matrix(data$GSM18971[i]);ndata[i,8] <- as.matrix(data$GSM18972[i])},{ndata[i,7] <- as.matrix(data$GSM18969[i]);ndata[i,8] <- as.matrix(data$GSM18970[i])})
      m1 <- mean(c(data$GSM18995[i],data$GSM18996[i]))
      m2 <- mean(c(data$GSM18947[i],data$GSM18948[i]))
      m<-which(c(m1,m2) == max(c(m1,m2)),arr.ind=TRUE)
      switch(m[1],{ndata[i,9] <- as.matrix(data$GSM18995[i]);ndata[i,10] <- as.matrix(data$GSM18996[i])},{ndata[i,9] <- as.matrix(data$GSM18947[i]);ndata[i,10] <- as.matrix(data$GSM18948[i])})
      m1<-mean(c(data$GSM18897[i],data$GSM18898[i]))
      m2<-mean(c(data$GSM18893[i],data$GSM18894[i]))
      m3<-mean(c(data$GSM18887[i],data$GSM18888[i]))
      m<-which(c(m1,m2,m3) == max(c(m1,m2,m3)),arr.ind=TRUE)
      switch(m[1],{ndata[i,11] <- as.matrix(data$GSM18897[i]);ndata[i,12] <- as.matrix(data$GSM18898[i])},{ndata[i,11] <- as.matrix(data$GSM18893[i]);ndata[i,12] <- as.matrix(data$GSM18894[i])},{ndata[i,11] <- as.matrix(data$GSM18887[i]);ndata[i,12] <- as.matrix(data$GSM18888[i])})
      m1<-mean(c(data$GSM18903[i],data$GSM18904[i]))
      m2<-mean(c(data$GSM18895[i],data$GSM18896[i]))
      m3<-mean(c(data$GSM18891[i],data$GSM18892[i]))
      m4<-mean(c(data$GSM18889[i],data$GSM18890[i]))
      m<-which(c(m1,m2,m3,m4) == max(c(m1,m2,m3,m4)),arr.ind=TRUE)
      switch(m[1],{ndata[i,13] <- as.matrix(data$GSM18903[i]);ndata[i,14] <- as.matrix(data$GSM18904[i])},{ndata[i,13] <- as.matrix(data$GSM18895[i]);ndata[i,14] <- as.matrix(data$GSM18896[i])},{ndata[i,13] <- as.matrix(data$GSM18891[i]);ndata[i,14] <- as.matrix(data$GSM18892[i])},{ndata[i,13] <- as.matrix(data$GSM18889[i]);ndata[i,14] <- as.matrix(data$GSM18890[i])})
      m1 <- mean(c(data$GSM18917[i],data$GSM18918[i]))
      m2 <- mean(c(data$GSM18915[i],data$GSM18916[i]))
      m<-which(c(m1,m2) == max(c(m1,m2)),arr.ind=TRUE)
      switch(m[1],{ndata[i,15] <- as.matrix(data$GSM18917[i]);ndata[i,16] <- as.matrix(data$GSM18918[i])},{ndata[i,15] <- as.matrix(data$GSM18915[i]);ndata[i,16] <- as.matrix(data$GSM18916[i])})
      m1 <- mean(c(data$GSM18919[i],data$GSM18920[i]))
      m2 <- mean(c(data$GSM18913[i],data$GSM18914[i]))
      m3 <- mean(c(data$GSM18925[i],data$GSM18926[i]))
      m4 <- mean(c(data$GSM19019[i],data$GSM19020[i]))
      m5 <- mean(c(data$GSM19021[i],data$GSM19022[i]))
      m6 <- mean(c(data$GSM18941[i],data$GSM18942[i]))
      m7 <- mean(c(data$GSM18937[i],data$GSM18938[i]))
      m8 <- mean(c(data$GSM18933[i],data$GSM18934[i]))
      m9 <- mean(c(data$GSM18935[i],data$GSM18936[i]))
      m10 <- mean(c(data$GSM18921[i],data$GSM18922[i]))
      m11 <- mean(c(data$GSM18927[i],data$GSM18928[i]))
      m12 <- mean(c(data$GSM18931[i],data$GSM18932[i]))
      m13 <- mean(c(data$GSM18923[i],data$GSM18924[i]))
      m14 <- mean(c(data$GSM18911[i],data$GSM18912[i]))
      m15 <- mean(c(data$GSM18939[i],data$GSM18940[i]))
      m16 <- mean(c(data$GSM18929[i],data$GSM18930[i]))
      m <- which(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16) == max(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)),arr.ind=TRUE)
      switch(m[1],{ndata[i,17] <- as.matrix(data$GSM18919[i]);ndata[i,18] <- as.matrix(data$GSM18920[i])},{ndata[i,17] <- as.matrix(data$GSM18913[i]);ndata[i,18] <- as.matrix(data$GSM18914[i])},{ndata[i,17] <- as.matrix(data$GSM18925[i]);ndata[i,18]<-as.matrix(data$GSM18926[i])},{ndata[i,17] <- as.matrix(data$GSM19019[i]);ndata[i,18] <- as.matrix(data$GSM19020[i])},{ndata[i,17] <- as.matrix(data$GSM19021[i]);ndata[i,18] <- as.matrix(data$GSM19022[i])},{ndata[i,17] <- as.matrix(data$GSM18941[i]);ndata[i,18] <- as.matrix(data$GSM18942[i])},{ndata[i,17] <- as.matrix(data$GSM18937[i]);ndata[i,18] <- as.matrix(data$GSM18938[i])},{ndata[i,17] <- as.matrix(data$GSM18933[i]);ndata[i,18] <- as.matrix(data$GSM18934[i])},{ndata[i,17] <- as.matrix(data$GSM18935[i]);ndata[i,18] <- as.matrix(data$GSM18936[i])},{ndata[i,17] <- as.matrix(data$GSM18921[i]);ndata[i,18] <- as.matrix(data$GSM18922[i])},{ndata[i,17] <- as.matrix(data$GSM18927[i]);ndata[i,18] <- as.matrix(data$GSM18928[i])},
      {ndata[i,17] <- as.matrix(data$GSM18931[i]);ndata[i,18] <- as.matrix(data$GSM18932[i])},{ndata[i,17] <- as.matrix(data$GSM18923[i]);ndata[i,18] <- as.matrix(data$GSM18924[i])},{ndata[i,17] <- as.matrix(data$GSM18911[i]);ndata[i,18] <- as.matrix(data$GSM18912[i])},{ndata[i,17] <- as.matrix(data$GSM18939[i]);ndata[i,18] <- as.matrix(data$GSM18940[i])},{ndata[i,17] <- as.matrix(data$GSM18929[i]);ndata[i,18] <- as.matrix(data$GSM18930[i])})
      m1<-mean(c(data$GSM19003[i],data$GSM19004[i]))
      m2<-mean(c(data$GSM19009[i],data$GSM19010[i]))
      m3<-mean(c(data$GSM19011[i],data$GSM19012[i]))
      m4<-mean(c(data$GSM19005[i],data$GSM19006[i]))
      m<-which(c(m1,m2,m3,m4) == max(c(m1,m2,m3,m4)),arr.ind=TRUE)
      switch(m[1],{ndata[i,19] <- as.matrix(data$GSM19003[i]);ndata[i,20] <- as.matrix(data$GSM19004[i])},{ndata[i,19] <- as.matrix(data$GSM19009[i]);ndata[i,20] <- as.matrix(data$GSM19010[i])},{ndata[i,19] <- as.matrix(data$GSM19011[i]);ndata[i,20] <- as.matrix(data$GSM19012[i])},{ndata[i,19] <- as.matrix(data$GSM19005[i]);ndata[i,20] <- as.matrix(data$GSM19006[i])})
      m1 <- mean(c(data$GSM18951[i],data$GSM18952[i]))
      m2 <- mean(c(data$GSM19007[i],data$GSM19008[i]))
      m<-which(c(m1,m2) == max(c(m1,m2)),arr.ind=TRUE)
      switch(m[1],{ndata[i,21] <- as.matrix(data$GSM18951[i]);ndata[i,22] <- as.matrix(data$GSM18952[i])},{ndata[i,21] <- as.matrix(data$GSM19007[i]);ndata[i,22] <- as.matrix(data$GSM19008[i])})
      ndata[i,23] <- as.matrix(data$GSM18975[i]); ndata[i,24] <- as.matrix(data$GSM18976[i])
      ndata[i,25] <- as.matrix(data$GSM18999[i]); ndata[i,26] <- as.matrix(data$GSM19000[i])
      ndata[i,27] <- as.matrix(data$GSM18909[i]); ndata[i,28] <- as.matrix(data$GSM18910[i])
      ndata[i,29] <- as.matrix(data$GSM18865[i]); ndata[i,30] <- as.matrix(data$GSM18866[i])
      ndata[i,31] <- as.matrix(data$GSM18955[i]); ndata[i,32] <- as.matrix(data$GSM18956[i])
      ndata[i,33] <- as.matrix(data$GSM18997[i]); ndata[i,34] <- as.matrix(data$GSM18998[i])
      ndata[i,35] <- as.matrix(data$GSM18967[i]); ndata[i,36] <- as.matrix(data$GSM18968[i])
      ndata[i,37] <- as.matrix(data$GSM18957[i]); ndata[i,38] <- as.matrix(data$GSM18958[i])
      ndata[i,39] <- as.matrix(data$GSM18991[i]); ndata[i,40] <- as.matrix(data$GSM18992[i])
      ndata[i,41] <- as.matrix(data$GSM19001[i]); ndata[i,42] <- as.matrix(data$GSM19002[i])
      ndata[i,43] <- as.matrix(data$GSM18943[i]); ndata[i,44] <- as.matrix(data$GSM18944[i])
      ndata[i,45] <- as.matrix(data$GSM18899[i]); ndata[i,46] <- as.matrix(data$GSM18900[i])
      ndata[i,47] <- as.matrix(data$GSM19017[i]); ndata[i,48] <- as.matrix(data$GSM19018[i])
      ndata[i,49] <- as.matrix(data$GSM18901[i]); ndata[i,50] <- as.matrix(data$GSM18902[i])
      ndata[i,51] <- as.matrix(data$GSM18993[i]); ndata[i,52] <- as.matrix(data$GSM18994[i])
      ndata[i,53] <- as.matrix(data$GSM18867[i]); ndata[i,54] <- as.matrix(data$GSM18868[i])
      ndata[i,55] <- as.matrix(data$GSM18961[i]); ndata[i,56] <- as.matrix(data$GSM18962[i])
      ndata[i,57] <- as.matrix(data$GSM18953[i]); ndata[i,58] <- as.matrix(data$GSM18954[i])
      ndata[i,59] <- as.matrix(data$GSM18949[i]); ndata[i,60] <- as.matrix(data$GSM18950[i])
}
#READ THE PROBESETS THAT ARE RELEVANT
probe_number = as.matrix(read.table("probe_number_human_gds596.txt",comment.char=""))
###########################################################################################################
#GET THE GENE EXPRESSION DATA FOR THE RELEVANT PROBESETS
DATA = matrix(rep(0,length(probe_number)*ncol(ndata)),nrow=length(probe_number))
row.names(DATA) <- probe_number
for(i in 1:length(probe_number)){
	for(j in 1:length(row.names(ndata))){
		if(probe_number[i] == row.names(ndata)[j]){
       DATA[i,] <- ndata[j,]
		}
	}
}
###########################################################################################################
#PERFORM MODIFIED BAYES FACTOR TO FIND TISSUE SPECIFIC, 2-SELECTIVE, 3-SELECTIVE AND 4-SELECTIVE PROBESETS
#IN THE GENE EXPRESSION DATA. NRREPLIC IS THE NUMBER OF REPLICATES FOR TISSUES IN THE DATA.
nrreplic = rep(2,30)
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
bf1 <- BF(DATA,nrreplic,1,nrruns,1.6)
#ASSIGN PROBESETS TO THE APPROPRIATE TISSUE SELECTIVE CATEGORY BASED ON THE THRESHOLD
for(i in 1:dat_row){
      # 63943 - 4 selective
      if(bfn4[i] > 27404){
              tissue_specific[i,max_tissue[i,1]] = 4
              tissue_specific[i,max_tissue[i,2]] = 4
              tissue_specific[i,max_tissue[i,3]] = 4
              tissue_specific[i,max_tissue[i,4]] = 4
      }
      # 16236 - 3 selective
      if(bfn3[i] > 4059 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 3
              tissue_specific[i,max_tissue[i,2]] = 3
              tissue_specific[i,max_tissue[i,3]] = 3
      }
      # 1736 - 2 selective
      if(bfn1[i] > 434 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 2
              tissue_specific[i,max_tissue[i,2]] = 2
      }
      if(bf1[i] > 27 && (tissue_specific[i,max_tissue[i,4]] == 0) && (tissue_specific[i,max_tissue[i,3]] == 0) && (tissue_specific[i,max_tissue[i,2]] == 0) && (tissue_specific[i,max_tissue[i,1]] == 0)){
              tissue_specific[i,max_tissue[i,1]] = 1
      }
}
#PRINT THE RESULTS IN TO A FILE FOR FURTHER PROCESSING
k <- 1
sink("human_genes_tissue_specificity_gds596_final.txt")
for(i in 1:length(h_gds_data)){
      for(j in 2:length(h_gds_data[[i]])){
          if(h_gds_data[[i]][j] != "0"){
              cat(h_gds_data[[i]][1],'\t',h_gds_data[[i]][j],'\t',tissue_specific[k,],'\n')
              k <- k+1
          }
      }
}
sink()
