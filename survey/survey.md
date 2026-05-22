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

#### 2.3. ggtree

A visualizer based on ape, therefore recycles its methods (no need to benchmark tree traversing)

#### 2.4. dendropy

#### 2.5. cogent3

Generating initial test trees with cogent3...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: True

Test tree structure:
((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);

Test tree tip names: ['A', 'B', 'C', 'D']

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 1.01 seconds

=== cogent3 Benchmark Statistics ===

preorder:
  Mean time (µs): 156.6650
  Median time (µs): 155.2700
  Std dev (µs): 5.1808
  Min time (µs): 151.2630
  Max time (µs): 225.8520

postorder:
  Mean time (µs): 417.4327
  Median time (µs): 415.3070
  Std dev (µs): 8.5776
  Min time (µs): 407.8930
  Max time (µs): 593.0200

levelorder:
  Mean time (µs): 165.7403
  Median time (µs): 164.3870
  Std dev (µs): 3.3601
  Min time (µs): 162.3540
  Max time (µs): 184.5750

is_ultrametric:
  Mean time (µs): 272.8828
  Median time (µs): 269.8000
  Std dev (µs): 10.8101
  Min time (µs): 262.4910
  Max time (µs): 419.0740

Plot saved to 'cogent3_traverse.png'


#### 2.6. biopython

>Terminal output:

Generating initial test trees with biopython...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: True

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
Benchmark completed in 94.18 seconds

=== biopython Benchmark Statistics ===

preorder:
  Mean time (µs): 30317.9912
  Median time (µs): 30182.3950
  Std dev (µs): 1604.8641
  Min time (µs): 28933.1090
  Max time (µs): 51965.5300

postorder:
  Mean time (µs): 31152.5805
  Median time (µs): 31124.3390
  Std dev (µs): 1904.8460
  Min time (µs): 29246.3460
  Max time (µs): 55589.4210

levelorder:
  Mean time (µs): 28840.6855
  Median time (µs): 28883.4305
  Std dev (µs): 1677.7494
  Min time (µs): 24819.7810
  Max time (µs): 53835.1430

is_ultrametric:
  Mean time (µs): 3749.9837
  Median time (µs): 3730.8465
  Std dev (µs): 92.3718
  Min time (µs): 3648.7880
  Max time (µs): 4572.5280

Plot saved to 'biopython_traverse.png'


#### 2.7. ete3

>Terminal output

**ISSUE: RETURNS TRUE FOR BOTH IS.ULTRAMETRIC; fix later**

Generating initial test trees with ete3...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: True

Test tree structure:

      /-A
   /-|
  |   "\"-B
--|
  |   /-C
   "\"-|
      "\""-D
Test tree leaf count: 4

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 85.71 seconds

=== ete3 Benchmark Statistics ===

preorder:
  Mean time (µs): 2165.6035
  Median time (µs): 2178.7040
  Std dev (µs): 47.4163
  Min time (µs): 2045.3240
  Max time (µs): 2490.1380

postorder:
  Mean time (µs): 7949.7280
  Median time (µs): 7921.4980
  Std dev (µs): 106.0309
  Min time (µs): 7696.7870
  Max time (µs): 9315.1710

levelorder:
  Mean time (µs): 2043.6679
  Median time (µs): 2037.7290
  Std dev (µs): 53.9969
  Min time (µs): 1980.8630
  Max time (µs): 2894.1540

is_ultrametric:
  Mean time (µs): 73466.2594
  Median time (µs): 72549.1290
  Std dev (µs): 2082.2927
  Min time (µs): 70841.1350
  Max time (µs): 79277.1680

Plot saved to 'ete3_traverse.png'

#### 2.8. toytree

Generating initial test trees with toytree...
Ultrametric test tree is ultrametric: True
Non-ultrametric test tree is ultrametric: True

Test tree structure:
<toytree.ToyTree at 0x7be0f1a49190>
Test tree leaf count: 4

Generating large ultrametric tree for benchmark...
Starting benchmarks...
Benchmark completed in 9.48 seconds

=== toytree Benchmark Statistics ===

preorder:
  Mean time (µs): 2109.4448
  Median time (µs): 2107.3090
  Std dev (µs): 19.0287
  Min time (µs): 2049.1200
  Max time (µs): 2221.6430

postorder:
  Mean time (µs): 1891.2966
  Median time (µs): 1892.2310
  Std dev (µs): 15.2445
  Min time (µs): 1840.6790
  Max time (µs): 1975.8530

levelorder:
  Mean time (µs): 1849.0401
  Median time (µs): 1847.5975
  Std dev (µs): 17.0936
  Min time (µs): 1798.4110
  Max time (µs): 1978.1570

is_ultrametric:
  Mean time (µs): 3614.3272
  Median time (µs): 3605.8925
  Std dev (µs): 71.7640
  Min time (µs): 3541.3570
  Max time (µs): 5260.7080

Plot saved to 'toytree_traverse.png'



### (3. phyl2vec?)
