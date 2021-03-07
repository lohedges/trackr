# TrackR

Track reconstruction in C++ using a Kalman filter implemented with
[Eigen](https://eigen.tuxfamily.org).

This repository is based on [code](https://github.com/dpohanlon/IPU4HEP)
by Daniel O'Hanlon and is intended to test and benchmark the track
reconstruction implementation described in [this](https://arxiv.org/pdf/2008.09210.pdf)
paper.

## Compiling

First clone the repository and its submodules:

```
git clone --recursive https://github.com/lohedges/trackr.git
```

Now navigate to the repository directory and run `make`:

```
cd trackr
make
```

After this, you should find three executables in the directory: `trackr_test`,
`trackr_benchmark`, and `trackr_benchmark_omp`.

By default the code is compiled using `g++`. To use a different compiler, pass
the `CXX` argument to `make`, e.g.:

```
make CXX=clang++
```

The code requires a `C++11` compliant compiler and is also compiled with `-O3`
optimisation flags by default. To pass additional flags, use the `OPTFLAGS`
argument, e.g.:

```
make OPTFLAGS=-funroll-loops
```

## Tests

TrackR has been numerically tested against two Python reference implementations
[here](https://github.com/dpohanlon/IPU4HEP/blob/master/kalman_filter_python/kf2d.py)
and [here](https://github.com/dpohanlon/IPU4HEP/blob/master/kalman_filter_python/kf2dTF.py)
using 5 tracks (`nGen` in the Python code.) While we use single precision for
all linear algebra operations, there is near perfect numerical agreement with
the output of the Python code.

To run the test:

```
./trackr_test
```

If the output agrees you should then see:

```
Hooray, test passed!
```

## Benchmarks

Benchmarks can be run using `trackr_benchmark`, e.g.:

```
./trackr_benchmark 10 100
```

The first argument is the number of tracks (input size), the second is the
number of repeats. The code performs repeat track reconstructions for the
chosen input size and outputs the mean throughput and variance. When finished
you should see output like the following:

```
      10          1.089306      0.005029            100
```

The columns are `tracks`, `throughput`, `variance`, `repeats`.

A benchmark script, `benchmark.sh`, is provided to test the performance of
the Eigen implementation across and range of input sizes. To run it:

```
./benchmark.sh
```

You should see output like the following:

```
# Tracks        Throughput      Variance        Repeats
       2          0.307449      0.000472         100000
       4          0.567414      0.000972         100000
       8          0.995351      0.003125         100000
      16          1.595489      0.011701         100000
      32          2.364217      0.032905         100000
      64          2.979356      0.040387         100000
     128          3.469605      0.088055         100000
     256          3.888345      0.064635         100000
     512          4.053215      0.073310         100000
    1024          3.720583      0.058653         100000
```

The following figure shows benchmark results run on a ThinkPad P1 Gen 2 laptop
with an [i7-9750H](https://www.intel.co.uk/content/www/uk/en/products/processors/core/i7-processors/i7-9750h.html)
CPU. Error bars show the standard deviation of the throughput.

![Benchmarks.](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_weetabix.png)

Here Clang is found to outperform GCC with a peak throughput of around 5.4
million tracks per second at an input size of roughly 500 tracks. For a larger
number of tracks, efficient parallelisation could likely be achieved by
running in batches of this input size.

As an example, the `trackr_benchmark_omp` script executes the Kalman filter in
parallel using [OpenMP](https://www.openmp.org). For example, we can run 5
batches of the 500 tracks in parallel and average over 100 repeats by running
the following:

```
./trackr_benchmark_omp 500 100 5
```

The arguments are the same as `trackr_benchmark`, with the addition of an extra
`batches` argument. You should see output like:

```
     500         12.323338      0.179636            100       5
```

The columns are `tracks`, `throughput`, `variance`, `repeats`, `batches`.

Another benchmark script, `benchmark_omp.sh`, is provided to test the
parallel performance of the Eigen implementation across range of batch numbers.
To run it:

```
./benchmark_omp.sh
```

You should see output like the following:

```
# Tracks        Throughput      Variance        Repeats Batches
     500          2.554266      0.012224          10000       1
     500          4.826631      0.233636          10000       2
     500          5.853855      1.882993          10000       3
     500          9.105968      0.595891          10000       4
     500          9.306524      3.740205          10000       5
     500          9.972789      3.556862          10000       6
     500         10.024893      0.877617          10000       7
     500         11.270907      0.718888          10000       8
     500         12.549682      0.891290          10000       9
     500         13.740646      1.473547          10000      10
     500         14.917983      1.105789          10000      11
     500         16.175130      1.780425          10000      12
     500         11.250949      0.613828          10000      13
     500         11.687452      0.973873          10000      14
     500         11.670631      1.574130          10000      15
     500         12.149416      1.644515          10000      16
     500         13.382016      2.169050          10000      17
     500         13.234795      1.609296          10000      18
     500         13.418732      0.811662          10000      19
     500         13.884903      1.036523          10000      20
     500         14.490800      1.371911          10000      21
     500         15.059707      1.426638          10000      22
     500         15.618644      1.961555          10000      23
     500         16.201504      1.981564          10000      24
     500         13.450061      1.224356          10000      25
     500         13.795785      1.301039          10000      26
     500         13.603895      1.664917          10000      27
     500         13.509158      6.982504          10000      28
     500         13.670539      1.370342          10000      29
     500         14.281145      1.704331          10000      30
```

If we wish to adjust the number of tracks to the optimum for your CPU, then
pass this as an argument to the script, e.g.:

```
./benchmark_omp.sh
```

The following figure shows benchmark results for parallel track reconstruction
run on the same ThinkPad P1 Gen 2 laptop, where error bars again show the
standard deviation of the throughput.

![Benchmarks.](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_weetabix_omp.png)

By wrapping the Kalman filter execution in a simple `#pramga omp parallel for`
loop, we have increased the peak throughput to around 26 million tracks per
second using 12 batches of 500 tracks running in parallel. Note that the
[i7-9750H](https://www.intel.co.uk/content/www/uk/en/products/processors/core/i7-processors/i7-9750h.html)
CPU has 6 cores and 12 threads. The throughput shows a sawtooth structure
with peaks at multiples of 12 batches, where all threads on the CPU are
fully active.


| **Device**    | **Cores** | **32 bit FLOPS**   | **TDP**  | **Cost**   | **Peak throughput**    |
|---------------|-----------|--------------------|----------|------------|------------------------|
| Graphcore GC2 | 1216      | 31.1 TFLOPS        | 120 W \* | $8112 \*\* | 2.6 x 10^6 tracks / s  |
| i7-9750H      | 6         | 0.7 TFLOPS \*\*\*  | 45 W     | $400       | 26 x 10^6 / tracks / s |

\* Two IPUs per board, so half the board TDP.

\*\* Cost is taken as a quarter of the price of an IPU M2000, containing 4 IPUs.

\*\*\* Single precision estimate at 3.6 GHz, which was sustained multicore frequency under test.

A rough comparison of the throughput per dollar and per watt for the current
implementations gives.

| **Device**    | **Throughput per watt** | **Throughput per dollar** |
|---------------|-------------------------|---------------------------|
| Graphcore GC2 | 2.17 x 10^4 tracks / s  | 321 tracks / s            |
| i7-9750H      | 2.17 x 10^5 tracks / s  | 65000 tracks / s          |

The CPU implementation is roughly an order of magnitude more power efficient
and two-orders of magnitude more cost efficient.
