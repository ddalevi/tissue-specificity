DATA <- strsplit(readLines("mouse_genes_tissue_specificity_final.txt",n=-1,ok=TRUE,warn=TRUE),"\n")
ndata <- matrix(rep(0,length(DATA)*3),nrow=length(DATA))
for(i in 1:length(DATA)){
      temp <- strsplit(DATA[[i]],"\t")
      if(length(temp[[1]]) > 1){
           if(max(strsplit(temp[[1]][3]," ")[[1]]) == "0"){
                 ndata[i,1] <- as.matrix(temp[[1]][1])
                 ndata[i,2] <- as.matrix(temp[[1]][2])
                 ndata[i,3] <- 0
           }
           if(max(strsplit(temp[[1]][3]," ")[[1]]) == "1"){
                 ndata[i,1] <- as.matrix(temp[[1]][1])
                 ndata[i,2] <- as.matrix(temp[[1]][2])
                 ndata[i,3] <- 1
          }
          if(max(strsplit(temp[[1]][3]," ")[[1]]) == "2"){
                 ndata[i,1] <- as.matrix(temp[[1]][1])
                 ndata[i,2] <- as.matrix(temp[[1]][2])
                 ndata[i,3] <- 2
          }
          if(max(strsplit(temp[[1]][3]," ")[[1]]) == "3"){
                 ndata[i,1] <- as.matrix(temp[[1]][1])
                 ndata[i,2] <- as.matrix(temp[[1]][2])
                 ndata[i,3] <- 3
          }
          if(max(strsplit(temp[[1]][3]," ")[[1]]) == "4"){
                 ndata[i,1] <- as.matrix(temp[[1]][1])
                 ndata[i,2] <- as.matrix(temp[[1]][2])
                 ndata[i,3] <- 4
          }
      }
      else{
          ndata[i,1] <- as.matrix(temp[[1]][1])
          ndata[i,2] <- "NULL"
          ndata[i,3] <- 0
      }
}
ndata <- ndata[-which(ndata == "0",arr.ind=TRUE)[,1],]
m_specific <- matrix(rep(0,length(which(ndata == "1",arr.ind=TRUE)[,1])*2),nrow=length(which(ndata == "1",arr.ind=TRUE)[,1]))
m_2sel <- matrix(rep(0,length(which(ndata == "2",arr.ind=TRUE)[,1])*2),nrow=length(which(ndata == "2",arr.ind=TRUE)[,1]))
m_3sel <- matrix(rep(0,length(which(ndata == "3",arr.ind=TRUE)[,1])*2),nrow=length(which(ndata == "3",arr.ind=TRUE)[,1]))
m_4sel <- matrix(rep(0,length(which(ndata == "4",arr.ind=TRUE)[,1])*2),nrow=length(which(ndata == "4",arr.ind=TRUE)[,1]))

k<-1
l<-1
m<-1
n<-1
for(i in 1:nrow(ndata)){
      if(ndata[i,3] == "1"){
          m_specific[k,1] <- ndata[i,1]
          m_specific[k,2] <- ndata[i,2]
          k<-k+1
      }
      else if(ndata[i,3] == "2"){
          m_2sel[l,1] <- ndata[i,1]
          m_2sel[l,2] <- ndata[i,2]
          l<-l+1
      }
      else if(ndata[i,3] == "3"){
          m_3sel[m,1] <- ndata[i,1]
          m_3sel[m,2] <- ndata[i,2]
          m<-m+1
      }
      else{
          m_4sel[n,1] <- ndata[i,1]
          m_4sel[n,2] <- ndata[i,2]
          n<-n+1
      }
}

write(m_specific[!duplicated(m_specific[,1]),1],"../tissue_selective/m9954_specific.txt") #EDIT THE DIRECTORY PATH
write(m_2sel[!duplicated(m_2sel[,1]),1],"../tissue_selective/m9954_2sel.txt") #EDIT THE DIRECTORY PATH
write(m_3sel[!duplicated(m_3sel[,1]),1],"../tissue_selective/m9954_3sel.txt") #EDIT THE DIRECTORY PATH
write(m_4sel[!duplicated(m_4sel[,1]),1],"../tissue_selective/m9954_4sel.txt") #EDIT THE DIRECTORY PATH


