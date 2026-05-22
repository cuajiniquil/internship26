library("ape")
library("ggplot2")

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

benchmarkSys <- function(trset, fn, times) {
  ttaken <- rep(0, times)
  #warmup:
  fn(trset[[1]])
  for (t in 1:times) {
    tstart <- Sys.time()
    fn(trset[[t]])
    tend <- Sys.time()
    ttaken[t] <- (tend - tstart) * 10^6
  }
  return(ttaken)
}

set.seed(273)
trset1 = lapply(1:1000, function(i) rcoal(1000))
trset2 = lapply(1:1000, function(i) rcoal(10000))
trset3 = lapply(1:1000, function(i) rcoal(100000))

start <- list()
totalt <- list()

start[[1]] <- Sys.time()
res1 <- benchmarkSys(trset1, is.ultrametric, 1000)
totalt[[1]]  <- Sys.time() - start[[1]]

start[[2]] <- Sys.time()
res2 <- benchmarkSys(trset2, is.ultrametric, 1000)
totalt[[2]]  <- Sys.time() - start[[2]]

start[[3]] <- Sys.time()
res3 <- benchmarkSys(trset3, is.ultrametric, 1000)
totalt[[3]]  <- Sys.time() - start[[3]]

apedf <- data.frame(
  traverse = rep(c("1k","100k","10M"), each = 1000),
  time = c(res1,res2,res3)
)

#================================
# 3: plotting and statistics
#================================
cat("\n=== ape Benchmark Statistics ===\n")
cat("Tree size: 1000\n")
summary(apedf$time[apedf$traverse =="1k"])
cat("Tree size: 100k\n")
summary(apedf$time[apedf$traverse =="100k"])
cat("Tree size: 10M\n")
summary(apedf$time[apedf$traverse =="10M"])
cat("\n=== Total run times ===\n")
totalt

ggplot(apedf,aes(x = traverse, y = time)) + geom_violin() + labs(x="Tree size", y="Time (micros)")

ggsave("ape_traverse.png", width = 8, height = 6)

