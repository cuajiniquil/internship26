library("castor")
library("ggplot2")

#SETWD AS BASH ARG
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

#================================
#1: initial tests : ADAPT IT TO CASTOR
#================================

#trs <- generate_gene_tree_msc(1000)
Ult <- read_tree(string = "((A:5,B:5):2,C:7);")
notU <- read_tree(string = "(((Pan:5,Homo:5):2, Chimp:6, Gorilla:7):4,Worm:10);")

plot(Ult)

#================================
#2: tree traversal benchmarks 
#================================

benchmarkSys <- function(trset, fn, times) {
  ttaken <- rep(0, times)
  #warmup?
  fn(trset[[1]])
  
  for (t in 1:times) {
    tstart <- Sys.time()
    fn(trset[[t]])
    tend <- Sys.time()
    ttaken[t] <- (tend - tstart) * 10^6
  }
  return(ttaken)
}

#had to add this fn in order to preserve the benchmarksys structure; 
#might slow it down - test against manual gad2tip? => res: still fast but does have a bigger gap for warmup
postorderfn <- function(tr){
  get_all_distances_to_tip(tr, focal_tip = 1)
}

#explanation in Running BEAST GUI... convo
isultramfn <- function(tr) {
  totdist <- get_all_distances_to_root(tr)
  return(diff(range(totdist)) < 1*10^-6) #10^-6 is the biggest specificity possible  
}


isultramfn(notU)
isultramfn(Ult)


ultra1k <- readLines("ultra1k.nwk")
ultra1k <- lapply(ultra1k, function(l) read_tree(string = l))

start <- list()
totalt <- list()

start[[1]] <- Sys.time()
res1 <- benchmarkSys(ultra1k, get_all_distances_to_root, 1000)
totalt[[1]]  <- Sys.time() - start[[1]]

start[[2]] <- Sys.time()
res2 <- benchmarkSys(ultra1k, postorderfn, 1000)
totalt[[2]]  <- Sys.time() - start[[2]]

start[[3]] <- Sys.time()
res3 <- rep(0,1000)
  #benchmarkSys(ultra1k, is.ultrametric, 1000)
totalt[[3]]  <- Sys.time() - start[[3]]

start[[4]] <- Sys.time()
res4 <- rep(0,1000)
  #benchmarkSys(ultra1k, is.ultrametric, 1000)
totalt[[4]]  <- Sys.time() - start[[4]]


castordf <- data.frame(
  traverse = rep(c("Preorder","Postorder","N/A (LEVELO)","N/A (ULTRAM)"), each = 1000),
  time = c(res1,res2,res3,res4)
)


#================================
# 3: plotting and statistics
#================================

cat("\n=== ape Benchmark Statistics ===\n")
cat("Preorder traversal \n")
summary(castordf$time[castordf$traverse =="Ultram 1k"])
cat("Postorder traversal \n")
summary(castordf$time[castordf$traverse =="Ultram 10k"])
cat("Future Levelorder\n")
summary(castordf$time[castordf$traverse =="Heteroc 1k"])
cat("Future isUltrametric\n")
summary(castordf$time[castordf$traverse =="Heteroc 10k"])
cat("\n=== Total run times for each type ===\n")
totalt

ggplot(castordf,aes(x = traverse, y = time)) + geom_violin() + labs(x="Tree size", y="Time (micros)")

ggsave("castor_traverse.png", width = 8, height = 6)


apemin <- castordf[castordf$time < 500,]

ggplot(apemin,aes(x = traverse, y = time)) + geom_violin() + labs(x="Castor Traversal Methods Comparison", y="Time (micros)")
ggsave("castor_traverse_manualtrim.png", width = 8, height = 6)