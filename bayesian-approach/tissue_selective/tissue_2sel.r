#CALCULATES THE OVERLAP OF 2-SELECTIVE GENES ACROSS FIVE DATASETS.
source("../commonfiles/RMDUP.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/RMDUP1.r") #EDIT THE DIRECTORY PATH
source("../commonfiles/OVP.r") #EDIT THE DIRECTORY PATH
h3113 <- as.matrix(read.table("h3113_2sel.txt"))
h596 <- as.matrix(read.table("h596_2sel.txt"))
h7307 <- as.matrix(read.table("h7307_2sel.txt"))
m9954 <- as.matrix(read.table("m9954_2sel.txt"))
r952 <- as.matrix(read.table("r952_2sel.txt"))
#STORING THE OVERLAP OF 2-SELECTIVE GENES IN 5 DIMENSIONAL MATRIX FOR FIVE DATASETS
sel2_overlap <- matrix(rep(0,5*5),nrow=5,ncol=5)

h3113_2sel <- RMDUP(h3113)
h596_2sel <- RMDUP(h596)
h7307_2sel <- RMDUP(h7307)
m9954_2sel <- RMDUP(m9954)
r952_2sel <- RMDUP(r952)

temp <- list(h3113_2sel,h596_2sel,h7307_2sel,m9954_2sel,r952_2sel)
#THE LOWER DIAGONAL OF MATRIX IS FILLED WITH OVERLAP BETWEEN DATASETS
for(i in 1:nrow(sel2_overlap)){
      for(j in 1:ncol(sel2_overlap)){
            if(i == j){
                 sel2_overlap[i,j] <- 1
            }
            if(i > j){
                 sel2_overlap[j,i] <- OVP(as.matrix(temp[[i]]),as.matrix(temp[[j]]))
            }
      }
}
#THE UPPER DIAGONAL OF MATRIX IS FILLED WITH PERCENTAGE OVERLAP BETWEEN DATASETS
for(i in 1:nrow(sel2_overlap)){
      for(j in 1:ncol(sel2_overlap)){
            if(i < j){
                 sel2_overlap[j,i] <- ( sel2_overlap[i,j] / ( length(temp[[i]]) + length(temp[[j]]) - sel2_overlap[i,j] ) ) * 100
            }
      }
}
temp1 <- rbind(h3113_2sel,h596_2sel,h7307_2sel,m9954_2sel,r952_2sel)
temp2 <- RMDUP1(temp1)
temp1 <- temp2
rm(temp2)
#CLUSTERING OF THE DATASETS
sel2_clust <- matrix(rep(0,length(temp1)*5),nrow=length(temp1),ncol=5)
row.names(sel2_clust) <- temp1
colnames(sel2_clust) <- c("Human(GDS3113)","Human(GDS596)","Human(GSE7307)","Mouse(GSE9954)","Rat(GSE952)")

for(i in 1:length(temp1)){
      for(j in 1:length(temp)){
            for(k in 1:length(temp[[j]])){
                  if( toupper(temp1[i]) == toupper(temp[[j]][k]) ){
                      sel2_clust[i,j] <- 1
                  }
            }
      }
}

dis <- dist(t(sel2_clust),method="binary")
jpeg("sel2_clust.jpg",quality = 100)
plot(hclust(dis),axes=TRUE,main="2-Selective Cluster",xlab=NULL,ylab="Height",sub=NULL)
dev.off()
