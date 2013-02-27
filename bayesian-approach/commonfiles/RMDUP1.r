RMDUP1 <- function(DATA){
 tdata <- as.matrix(DATA[!duplicated(DATA)])
 for(i in 1:nrow(tdata)){
       for(j in 1:nrow(tdata)){
             if( (i != j) && (toupper(tdata[i]) == toupper(tdata[j])) && (tdata[i] != "Inf") && (tdata[j] !="Inf") ){
                 tdata[j] = "Inf"
             }
       }
 }
 tdata <- tdata[-which(tdata == "Inf",arr.ind=TRUE)[,1],]
 RMDUP1 <- tdata
}
