#OVERLAP CALCULATOR
OVP <- function(DATA1,DATA2){
 tdata1 <- DATA1
 tdata2 <- DATA2
 count <- 0
 for(i in 1:nrow(tdata1)){
      for(j in 1:nrow(tdata2)){
           if( toupper(tdata1[i]) == toupper(tdata2[j]) ){
                count <- count + 1
           }
      }
 }
 OVP <- count
}
