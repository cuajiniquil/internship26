library("ape")


#================================
#1: initial tests
#================================

trs <- rcoal(1000)
notU <- read.tree(text = "(((Pan:5,Homo:5):2, Chimp:6, Gorilla:7):4,Worm:10);")

is.ultrametric(trs)
is.ultrametric(notU)

#================================
#2: tree traversing benchmarks
#================================

benchmarkSys <- function(fnctn, times) {
  ttaken <- rep(0,times)
  #+1 to avoid returning the "warmup" iteration
  for (t in 1:(times+1)) {
    tstart <- Sys.time()
    fnctn
    tend <- Sys.time()
    ttaken[t] <- (tend - tstart) * 10^6
  }
  #plot(ttaken)
  return(ttaken[2:times+1])
}

#two things to fix:
#1 - the first iteration always takes much longer (see warmup workarounds)
#2 - autoplot doesn't work since ggplot's 'aes()' is now depracated (see proposed function change in githun issues)
benchmarkMicro <- function(fnctn,times){
  library("microbenchmark")
  #library("ggplot2")
  total <- microbenchmark(fnctn, times=times)
  ggplot2::autoplot(total)
  return(total)
}

set.seed(273)
sstart <- Sys.time()
sres <- benchmarkSys(is.ultrametric(rcoal(100000)),1000)
send  <-Sys.time()


apedf <- data.frame(
  traverse = "tree traversal",
  time = sres
)

#================================
# 3: plotting and statistics
#================================

ggplot(apedf,aes(x = traverse, y = time)) + geom_violin() + labs(x="Density", y="Time (ns)")
mx <- which.max(sres)
sres[mx]
