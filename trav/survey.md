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
    bonus:
    benchmark(preorder(ultrametricTree(m)),n)    
    benchmark(postorder(ultrametricTree(m)),n)
    benchmark(levelorder(ultrametricTree(m)),n)
    
#### 2.2. ape

The first result is a "warmup" that skews the average. 
> sres[1]
[1] 1051343
> sres[2]
[1] 1.192093

Fix: do 1:t+1 for loop and return 2:t+1 results. Ther


A similar issue with no current fix has been found in microbenchmarking (https://github.com/joshuaulrich/microbenchmark/issues/53). Additionally, 

> mres$time[2:21] - sres[1:20]
 [1] 239.046326   9.046326   9.046326   8.807907   9.046326  19.046326   9.284744  10.284744  19.284744
[10]   8.569489   8.807907   9.284744   9.046326   9.284744   9.046326   9.046326   8.807907  19.046326
[19]   9.046326   9.284744

For these reasons, the function used from here on out is benchmarkSys.

for sres

> max(sres)
[1] 5.245209

Results summary:

=== ape Benchmark Statistics ===

Ultrametric 1k tips
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  151.2   203.0   213.3   259.6   223.6 45334.8 

Ultrametric 10k tips
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1299    1330    1341    1350    1357    2021 
   
Heterochronic 1k tips
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  157.6   170.5   175.2   204.2   183.2 19997.1 

Heterochronic 10k tips
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1231    1368    1386    1413    1403   23197 

=== Total run times for each type ===

[[1]]
Time difference of 0.2854817 secs

[[2]]
Time difference of 1.370536 secs

[[3]]
Time difference of 0.2203369 secs

[[4]]
Time difference of 1.43401 secs


#### 2.3. ggtree

A visualizer based on ape, therefore recycles its methods (no need to benchmark tree traversing)

#### 2.4. dendropy

Generating initial test trees...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: False

Test tree root: <Node object at 0x7d0515f52b10: 'None' (None)>
Test tree number of leaves: 1000

Extracting trees from nwk files...

Starting benchmarks...
Benchmark completed in 1.43 seconds

=== dendropy Benchmark Statistics ===

preorder:
  Mean time (µs): 441.0897
  Median time (µs): 387.0245
  Std dev (µs): 71.9765
  Min time (µs): 358.2600
  Max time (µs): 629.2680

postorder:
  Mean time (µs): 359.9071
  Median time (µs): 359.1125
  Std dev (µs): 3.5262
  Min time (µs): 354.4240
  Max time (µs): 398.5260

levelorder:
  Mean time (µs): 207.7913
  Median time (µs): 207.0220
  Std dev (µs): 2.3586
  Min time (µs): 203.5710
  Max time (µs): 220.8430

is_ultrametric:
  Mean time (µs): 421.3205
  Median time (µs): 421.8295
  Std dev (µs): 9.8495
  Min time (µs): 404.5380
  Max time (µs): 594.3320

Plot saved to 'dendropy_traverse.png'

#### 2.5. cogent3

Generating initial test trees with cogent3...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: False

Test tree structure:
((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);

Test tree tip names: ['A', 'B', 'C', 'D']

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 0.98 seconds

=== cogent3 Benchmark Statistics ===

preorder:
  Mean time (µs): 148.4441
  Median time (µs): 147.9575
  Std dev (µs): 2.2507
  Min time (µs): 145.0720
  Max time (µs): 169.5870

postorder:
  Mean time (µs): 404.0007
  Median time (µs): 403.6060
  Std dev (µs): 3.2119
  Min time (µs): 398.0350
  Max time (µs): 425.3370

levelorder:
  Mean time (µs): 160.7974
  Median time (µs): 160.0200
  Std dev (µs): 2.1355
  Min time (µs): 158.2660
  Max time (µs): 189.1750

is_ultrametric:
  Mean time (µs): 260.5246
  Median time (µs): 259.1950
  Std dev (µs): 4.8786
  Min time (µs): 252.7330
  Max time (µs): 307.4160

Plot saved to 'cogent3_traverse.png'

#### 2.6. biopython

Generating initial test trees with biopython...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: False

Test tree structure:
                                        _____________________________________ A
  _____________________________________|
 |                                     |_____________________________________ B
_|
 |                                      _____________________________________ C
 |_____________________________________|
                                       |_____________________________________ D

Test tree leaf count: 4

Extracting trees from nwk files...
Starting benchmarks...
Benchmark completed in 6.32 seconds

=== biopython Benchmark Statistics ===

preorder:
  Mean time (µs): 2927.1824
  Median time (µs): 2873.5250
  Std dev (µs): 595.5434
  Min time (µs): 2769.8110
  Max time (µs): 15863.9040

postorder:
  Mean time (µs): 2883.0258
  Median time (µs): 2867.8440
  Std dev (µs): 501.9825
  Min time (µs): 2788.8970
  Max time (µs): 14118.6600

levelorder:
  Mean time (µs): 133.5562
  Median time (µs): 132.8690
  Std dev (µs): 2.6847
  Min time (µs): 129.8530
  Max time (µs): 170.9100

is_ultrametric:
  Mean time (µs): 363.8407
  Median time (µs): 363.2345
  Std dev (µs): 5.4136
  Min time (µs): 355.2440
  Max time (µs): 427.1190

Plot saved to 'biopython_traverse.png'

#### 2.7. ete3

Generating initial test trees with ete3...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: False

Test tree structure:

      /-A
   /-|
  |   '\'-B
--|
  |   /-C
   '\'-|
      '\'-D
Test tree leaf count: 4

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 82.25 seconds

=== ete3 Benchmark Statistics ===

preorder:
  Mean time (µs): 2196.3424
  Median time (µs): 2169.5270
  Std dev (µs): 120.5144
  Min time (µs): 2074.4540
  Max time (µs): 3266.6260

postorder:
  Mean time (µs): 7921.4069
  Median time (µs): 7769.9080
  Std dev (µs): 357.4379
  Min time (µs): 7502.4470
  Max time (µs): 11253.3690

levelorder:
  Mean time (µs): 1997.7605
  Median time (µs): 1976.1700
  Std dev (µs): 93.9674
  Min time (µs): 1897.3430
  Max time (µs): 3032.6280

is_ultrametric:
  Mean time (µs): 70055.5502
  Median time (µs): 70035.7010
  Std dev (µs): 789.6765
  Min time (µs): 68188.6480
  Max time (µs): 77554.9050

Plot saved to 'ete3_traverse.png'

#### 2.8. toytree

Generating initial test trees with toytree...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: False

Test tree structure:
<toytree.ToyTree at 0x789ff01cddf0>
Test tree leaf count: 4

Extracting trees from nwk files...
Starting benchmarks...
Benchmark completed in 1.14 seconds

=== toytree Benchmark Statistics ===

preorder:
  Mean time (µs): 206.6023
  Median time (µs): 205.9950
  Std dev (µs): 2.2645
  Min time (µs): 201.4370
  Max time (µs): 226.0130

postorder:
  Mean time (µs): 177.3530
  Median time (µs): 172.0920
  Std dev (µs): 9.7117
  Min time (µs): 166.5810
  Max time (µs): 218.2290

levelorder:
  Mean time (µs): 168.5053
  Median time (µs): 167.9745
  Std dev (µs): 3.6049
  Min time (µs): 164.4070
  Max time (µs): 256.3400

is_ultrametric:
  Mean time (µs): 585.8150
  Median time (µs): 585.3355
  Std dev (µs): 6.6104
  Min time (µs): 572.6420
  Max time (µs): 649.8260

Plot saved to 'toytree_traverse.png'

#### 2.9 Phylo.jl

Generating initial test trees with Phylo.jl...
Ultrametric test tree is ultrametric: true
Non-ultrametric test tree is ultrametric: false

Test tree leaf count: 4

Extracting trees from nwk files...
Starting benchmarks...
Benchmark completed in 1.05 seconds

=== Phylo.jl Benchmark Statistics ===

preorder:
  Mean time (µs): 132.2183
  Median time (µs): 116.4985
  Std dev (µs): 370.7202
  Min time (µs): 110.147
  Max time (µs): 11831.304

postorder:
  Mean time (µs): 161.1209
  Median time (µs): 120.912
  Std dev (µs): 922.1456
  Min time (µs): 111.76
  Max time (µs): 29279.99

levelorder:
  Mean time (µs): 301.4037
  Median time (µs): 121.132
  Std dev (µs): 5351.8917
  Min time (µs): 109.956
  Max time (µs): 169372.424

is_ultrametric:
  Mean time (µs): 259.694
  Median time (µs): 222.1015
  Std dev (µs): 457.4447
  Min time (µs): 206.207
  Max time (µs): 10863.949

Plot saved to 'phylo_jl_traverse.png'

#### 2.10. castor (R)

### (3. phyl2vec?)
