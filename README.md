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

After this, you should find two executables in the directory: `track_test` and
`track_benchmark`.

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
