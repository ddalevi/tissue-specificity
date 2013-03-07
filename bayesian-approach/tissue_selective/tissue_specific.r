#CALCULATES THE OVERLAP OF TISSUE SPECIFIC GENES ACROSS FIVE DATASETS.
source("../commonfiles/RMDUP.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/RMDUP1.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/OVP.r") #EDIT THE DIRECTORY PATH
h3113 <- as.matrix(read.table("h3113_specific.txt"))
h596 <- as.matrix(read.table("h596_specific.txt"))
h7307 <- as.matrix(read.table("h7307_specific.txt"))
m9954 <- as.matrix(read.table("m9954_specific.txt"))
r952 <- as.matrix(read.table("r952_specific.txt"))
#STORING THE OVERLAP OF TISSUE SPECIFIC GENES IN 5 DIMENSIONAL MATRIX FOR FIVE DATASETS
specific_overlap <- matrix(rep(0,5*5),nrow=5,ncol=5)

h3113_specific <- RMDUP(h3113)
h596_specific <- RMDUP(h596)
h7307_specific <- RMDUP(h7307)
m9954_specific <- RMDUP(m9954)
r952_specific <- RMDUP(r952)

temp <- list(h3113_specific,h596_specific,h7307_specific,m9954_specific,r952_specific)
#THE LOWER DIAGNOL OF MATRIX IS FILLED WITH OVERLAP BETWEEN DATASETS
for(i in 1:nrow(specific_overlap)){
      for(j in 1:ncol(specific_overlap)){
            if(i == j){
                 specific_overlap[i,j] <- 1
            }
            if(i > j){
                 specific_overlap[j,i] <- OVP(as.matrix(temp[[i]]),as.matrix(temp[[j]]))
            }
      }
}
#TJE UPPER DIAGONAL OF MATRIX IS FILLED WITH PERCENTAGE OVERLAP BETWEEN DATASETS
for(i in 1:nrow(specific_overlap)){
      for(j in 1:ncol(specific_overlap)){
            if(i < j){
                 specific_overlap[j,i] <-  ( specific_overlap[i,j] / ( length(temp[[i]]) + length(temp[[j]]) - specific_overlap[i,j] ) ) * 100
            }
      }
}

temp1 <- rbind(h3113_specific,h596_specific,h7307_specific,m9954_specific,r952_specific)
temp2 <- RMDUP1(temp1)
temp1 <- temp2
rm(temp2)
#CLUSTERING OF THE DATASETS
specific_clust <- matrix(rep(0,length(temp1)*5),nrow=length(temp1),ncol=5)
row.names(specific_clust) <- temp1
colnames(specific_clust) <- c("Human(GDS3113)","Human(GDS596)","Human(GSE7307)","Mouse(GSE9954)","Rat(GSE952)")

for(i in 1:length(temp1)){
      for(j in 1:length(temp)){
            for(k in 1:length(temp[[j]])){
                  if( toupper(temp1[i]) == toupper(temp[[j]][k]) ){
                      specific_clust[i,j] <- 1
                  }
            }
      }
}

dis <- dist(t(specific_clust),method="binary")
jpeg("specific_clust.jpg",quality = 100)
plot(hclust(dis),axes=TRUE,main="Specific Cluster",xlab=NULL,ylab="Height",sub=NULL)
dev.off()




