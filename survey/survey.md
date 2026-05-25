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


A similar issue with no current fix has been found in microbenchmarking (https://github.com/joshuaulrich/microbenchmark/issues/53).

> mres$time[2:21] - sres[1:20]
 [1] 239.046326   9.046326   9.046326   8.807907   9.046326  19.046326   9.284744  10.284744  19.284744
[10]   8.569489   8.807907   9.284744   9.046326   9.284744   9.046326   9.046326   8.807907  19.046326
[19]   9.046326   9.284744

For these reasons, the function used from here on out is benchmarkSys.

for sres

> max(sres)
[1] 5.245209

#### 2.3. ggtree

A visualizer based on ape, therefore recycles its methods (no need to benchmark tree traversing)

#### 2.4. dendropy

Generating initial test trees...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: False

Test tree root: <Node object at 0x786d7393e720: 'None' (None)>
Test tree number of leaves: 1000

Starting benchmarks...
Benchmark completed in 17.06 seconds

=== dendropy Benchmark Statistics ===

preorder:
  Mean time (µs): 4465.6845
  Median time (µs): 4450.2000
  Std dev (µs): 296.3414
  Min time (µs): 4033.3840
  Max time (µs): 6660.0240

postorder:
  Mean time (µs): 3835.9992
  Median time (µs): 3824.1830
  Std dev (µs): 93.3559
  Min time (µs): 3725.0170
  Max time (µs): 5760.2990

levelorder:
  Mean time (µs): 3789.8141
  Median time (µs): 3777.7300
  Std dev (µs): 71.3077
  Min time (µs): 3686.2640
  Max time (µs): 4459.5830

is_ultrametric:
  Mean time (µs): 4953.9725
  Median time (µs): 4943.2780
  Std dev (µs): 63.0396
  Min time (µs): 4823.2640
  Max time (µs): 5363.0050

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

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 67.51 seconds

=== biopython Benchmark Statistics ===

preorder:
  Mean time (µs): 28054.4702
  Median time (µs): 28041.5680
  Std dev (µs): 1249.2541
  Min time (µs): 26635.7930
  Max time (µs): 41059.8290

postorder:
  Mean time (µs): 28477.1477
  Median time (µs): 28391.1860
  Std dev (µs): 1083.6928
  Min time (µs): 28063.8840
  Max time (µs): 41767.1810

levelorder:
  Mean time (µs): 7359.9267
  Median time (µs): 7353.7135
  Std dev (µs): 42.9094
  Min time (µs): 7269.4550
  Max time (µs): 7551.5240

is_ultrametric:
  Mean time (µs): 3525.9377
  Median time (µs): 3520.1175
  Std dev (µs): 36.3710
  Min time (µs): 3458.8080
  Max time (µs): 3696.2130

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
<toytree.ToyTree at 0x71cf36b6faa0>
Test tree leaf count: 4

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 11.50 seconds

=== toytree Benchmark Statistics ===

preorder:
  Mean time (µs): 2086.0348
  Median time (µs): 2083.2180
  Std dev (µs): 27.1078
  Min time (µs): 2036.3880
  Max time (µs): 2620.2790

postorder:
  Mean time (µs): 1857.8428
  Median time (µs): 1857.0380
  Std dev (µs): 19.1475
  Min time (µs): 1811.0340
  Max time (µs): 1948.3860

levelorder:
  Mean time (µs): 1790.2446
  Median time (µs): 1786.3290
  Std dev (µs): 19.0562
  Min time (µs): 1755.3910
  Max time (µs): 1886.1220

is_ultrametric:
  Mean time (µs): 5755.9146
  Median time (µs): 5754.3930
  Std dev (µs): 46.4490
  Min time (µs): 5624.8860
  Max time (µs): 6012.2280

Plot saved to 'toytree_traverse.png'

#### 2.9 Phylo.jl

Generating initial test trees with Phylo.jl...
Ultrametric test tree is ultrametric: true
Non-ultrametric test tree is ultrametric: false

Test tree leaf count: 4

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 7.99 seconds

=== Phylo.jl Benchmark Statistics ===

preorder:
  Mean time (µs): 266.1369
  Median time (µs): 48.19
  Std dev (µs): 4988.5093
  Min time (µs): 43.732
  Max time (µs): 154932.34

postorder:
  Mean time (µs): 2891.3937
  Median time (µs): 2605.087
  Std dev (µs): 915.1422
  Min time (µs): 2397.22
  Max time (µs): 13842.358

levelorder:
  Mean time (µs): 1429.3421
  Median time (µs): 1121.1755
  Std dev (µs): 5151.6154
  Min time (µs): 1080.685
  Max time (µs): 161629.768

is_ultrametric:
  Mean time (µs): 3216.5896
  Median time (µs): 2892.338
  Std dev (µs): 1034.2751
  Min time (µs): 2727.295
  Max time (µs): 11723.618

Plot saved to 'phylo\_jl\_traverse.png'

### (3. phyl2vec?)
