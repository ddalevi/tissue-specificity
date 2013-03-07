#REMOVES THE DUPLICATE GENE SYMBOLS IN A DATASET AND RETURNS UNIQUE GENE SYMBOL LIST
RMDUP <- function(DATA){
 tdata <- DATA
 temp <- matrix()
 for(i in 1:nrow(tdata)){
      if(regexpr("^(.*)///",tdata[i],perl=TRUE) == TRUE){
           temp <- rbind(temp,tdata[i])
           tdata[i] <- "Inf"

      }
 }
 temp <- temp[-1,]
 if( length(which(tdata == "Inf",arr.ind=TRUE)) != 0){
     tdata <- tdata[-which(tdata == "Inf",arr.ind=TRUE)[,1],]
 }
 if(length(temp!=0)){
  temp <- strsplit(temp,"///")
  for(i in 1:length(temp)){
      for(j in 1:length(temp[[i]])){
            for(k in 1:length(tdata)){
                  if(temp[[i]][j] == tdata[k]){
                        temp[[i]][j] <- "NULL"
                  }
            }
      }
  }
  for(i in 1:length(temp)){
      for(j in 1:length(temp)){
            for(k in 1:length(temp[[i]])){
                  for(l in 1:length(temp[[j]])){
                        if((temp[[i]][k] != "NULL") && (temp[[j]][l] != "NULL") && (temp[[i]][k] == temp[[j]][l]) && (i!=j)){
                            temp[[i]][k] <- "NULL"
                        }
                  }
            }
      }
  }
  for(i in 1:length(temp)){
      count <- 0
      for(j in 1:length(temp[[i]])){
            if(temp[[i]][j] == "NULL"){
                  count <- count + 1
            }
      }
      if(count != 0){
            temp1 <- temp[[i]]
            temp1 <- temp1[-(which(temp1 == "NULL",arr.ind=TRUE))]
            temp[[i]] <- temp1
      }
  }

  store <- matrix()
  for(i in 1:length(temp)){
      if(length(temp[[i]]) == 0){
            store <- rbind(store,i)
      }
  }
 
  store <- store[-1]
  if(length(store) > 0){
      for(i in length(store):1){
            temp[[store[i]]] <- NULL
      }
  }
  tdata <- as.matrix(tdata)
  for(i in 1:length(temp)){
      for(j in 1:length(temp[[i]])){
            tdata <- rbind(tdata,as.matrix(temp[[i]][j]))
      }
  }
 }
 RMDUP <- tdata
}
