library("ape")
library("ggplot2")

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
setwd(args[1])

#================================
#1: initial tests
#================================

trs <- rcoal(1000)
notU <- read.tree(text = "(((Pan:5,Homo:5):2, Chimp:6, Gorilla:7):4,Worm:10);")

is.ultrametric(trs)
is.ultrametric(notU)

#================================
#2: data generation (in comments, eventually make different file)
#================================

#set.seed(273)

#trsetu1 = lapply(1:1000, function(i) rcoal(1000))
#trsetu10 = lapply(1:1000, function(i) rcoal(10000))
#trseth1 = lapply(1:1000, function(i) rtree(1000))
#trseth10 = lapply(1:1000, function(i) rtree(10000))

#write.tree(trsetu1, file = "ultra1k.nwk")
#write.tree(trsetu10, file = "ultra10k.nwk")
#write.tree(trseth1, file = "hetero1k.nwk")
#write.tree(trseth10, file = "hetero10k.nwk")

#================================
#3: tree traversal benchmarks
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

ultra1k <- read.tree("ultra1k.nwk")
ultra10k <- read.tree("ultra10k.nwk")
hetero1k <- read.tree("hetero1k.nwk")
hetero10k <- read.tree("hetero10k.nwk")

start <- list()
totalt <- list()

start[[1]] <- Sys.time()
res1 <- benchmarkSys(ultra1k, is.ultrametric, 1000)
totalt[[1]]  <- Sys.time() - start[[1]]

start[[2]] <- Sys.time()
res2 <- benchmarkSys(ultra10k, is.ultrametric, 1000)
totalt[[2]]  <- Sys.time() - start[[2]]

start[[3]] <- Sys.time()
res3 <- benchmarkSys(hetero1k, is.ultrametric, 1000)
totalt[[3]]  <- Sys.time() - start[[3]]

start[[4]] <- Sys.time()
res4 <- benchmarkSys(hetero10k, is.ultrametric, 1000)
totalt[[4]]  <- Sys.time() - start[[4]]


apedf <- data.frame(
  traverse = rep(c("Ultram 1k","Ultram 10k","Heteroc 1k","Heteroc 10k"), each = 1000),
  time = c(res1,res2,res3,res4)
)


#================================
# 4: plotting and statistics
#================================

cat("\n=== ape Benchmark Statistics ===\n")
cat("Ultrametric 1k tips\n")
summary(apedf$time[apedf$traverse =="Ultram 1k"])
cat("Ultrametric 10k tips\n")
summary(apedf$time[apedf$traverse =="Ultram 10k"])
cat("Heterochronic 1k tips\n")
summary(apedf$time[apedf$traverse =="Heteroc 1k"])
cat("Heterochronic 10k tips\n")
summary(apedf$time[apedf$traverse =="Heteroc 10k"])
cat("\n=== Total run times for each type ===\n")
totalt

ggplot(apedf,aes(x = traverse, y = time)) + geom_violin() + labs(x="Tree size", y="Time (micros)")

ggsave("ape_traverse.png", width = 8, height = 6)


apemin <- apedf[apedf$time < 5000,]

ggplot(apemin,aes(x = traverse, y = time)) + geom_violin() + labs(x="Tree size", y="Time (micros)")
ggsave("ape_traverse_manualtrim.png", width = 8, height = 6)


cat("Test if these results are significant:")

#what tests could be done? ASK AGAIN BUT OUT OF CONTEXT (give generic description)
#kruskal.test(time ~ traverse, data = apedf)
