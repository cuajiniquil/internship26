#A survey of phylogenetic tools' functionalities:

---

### 1. Tree


### 2. Tree traversing benchmarks

#### 2.1. Pseudocode

Measures how long it takes to traverse n trees of size m 

    benchmark(f,t):
    for i in range t
        times[t] = measure_time(f)
    end for
    return times
    end benchmark
    
    benchmark(isUltrametric(ultrametricTree(m)),n)

#### 2.2. ape

The first result is a "warmup" that skews the average. 
> sres[1]
[1] 1051343
> sres[2]
[1] 1.192093

Fix: do 1:t+1 for loop and return 2:t+1 results. Ther


A similar issue with no current fix has been found in microbenchmarking (https://github.com/joshuaulrich/microbenchmark/issues/53).

> mres$time[2:21] - sres[1:20]
 [1] 239.046326   9.046326   9.046326   8.807907   9.046326  19.046326   9.284744  10.284744  19.284744
[10]   8.569489   8.807907   9.284744   9.046326   9.284744   9.046326   9.046326   8.807907  19.046326
[19]   9.046326   9.284744

For these reasons, the function used from here on out is benchmarkSys.

for sres

> max(sres)
[1] 5.245209

### (3. phyl2vec?)
