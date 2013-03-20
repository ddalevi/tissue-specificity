#######Shannon Entropy########################################
#Please source this file before using ROKUspm.r or ROKU.r#####

tukey.biweight <- function(x, c=5, epsilon=0.0001)
  {
    m <- median(x)
    s <- median(abs(x - m))
    u <- (x - m) / (c * s + epsilon)
    w <- rep(0, length(x))
    i <- abs(u) <= 1
    w[i] <- ((1 - u^2)^2)[i]
    t.bi <- sum(w * x) / sum(w)
    return(t.bi)
  }

probablity.tbw <- function(x)
{
    xa<- tukey.biweight(x)
    xt<-abs(x-xa)
    s<-sum(xt)
    pt<-xt/s
    return(pt)
}

shannon.entropy <- function(x)
{     
      p<-probablity.tbw(x)
  if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
}

#######output a histogram of entropy distribution########### 
ENTROPYsummary <- function(DATA){

    enall <- rep(0,nrow(DATA))
    for (i in 1:nrow(DATA)){
         enall[i] <- shannon.entropy(DATA[i,])
    }
    plot(sort(enall), main="ENTROPY BioGPS", xlab = "num of genes", ylab = "entropy")
    ENTROPYsummary <- c(mean(enall),quantile(enall,prob=c(0,.25,.5,.8,1)))
    windows()
    hist(enall,100,main = "Histogram of BioGPS", xlab="Entropy")
    return (ENTROPYsummary)
}

