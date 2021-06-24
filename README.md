# TrackR

Track reconstruction in C++ using a Kalman filter. Implemented with
[Eigen](https://eigen.tuxfamily.org) on the CPU and using the
[Poplar SDK](https://www.graphcore.ai/products/poplar) on the
[IPU](https://www.graphcore.ai).

This repository is based on [code](https://github.com/dpohanlon/IPU4HEP)
by Daniel O'Hanlon and is intended to test and benchmark the track
reconstruction implementation described in [this](https://link.springer.com/article/10.1007/s41781-021-00057-z)
paper.

## CPU implementation

### Compiling

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

When using Clang you'll need to `libomp` package for OpenMP support. If
compiling on macOS you'll also explicitly need to pass `CXX=clang++` to
use Clang, since Apple use a GCC frontend to an LLVM backend.

The code requires a `C++11` compliant compiler and is also compiled with `-O3`
optimisation flags by default. To pass additional flags, use the `OPTFLAGS`
argument, e.g.:

```
make OPTFLAGS=-funroll-loops
```

### Testing

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

### Benchmarks

(Note that benchmarks are performed over different realisations of the track
data. While it shouldn't affect results in this case, i.e. the speed of the
computation doesn't depend on the data itself, it is nonetheless good practice.)

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

![Benchmarks CPU.](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_weetabix.png)

Here Clang is found to outperform GCC with a peak throughput of around 5.4
million tracks per second at an input size of roughly 500 tracks. For a larger
number of tracks, efficient parallelisation can be achieved by running in batches
of this input size.

As an example, the `trackr_benchmark_omp` binary executes the Kalman filter in
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
./benchmark_omp.sh 100
```

The following figure shows benchmark results for parallel track reconstruction
run on the same ThinkPad P1 Gen 2 laptop, where error bars again show the
standard deviation of the throughput.

![Benchmarks CPU OpenMP.](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_weetabix_omp.png)

By wrapping the Kalman filter execution in a simple `#pramga omp parallel for`
loop the peak throughput is increased to around 26 million tracks per second
when using 12 batches of 500 tracks running in parallel. Note that the
[i7-9750H](https://www.intel.co.uk/content/www/uk/en/products/processors/core/i7-processors/i7-9750h.html)
CPU has 6 cores and 12 threads. The throughput shows a sawtooth structure
with peaks at multiples of 12 batches, where all threads on the CPU are
fully active.

(Note that the single batch performance is lower than that obtained using
`trackr_benchmark`. This is due to the overhead of using OpenMP.)

| **Device**    | **Cores** | **32 bit FLOPS**   | **TDP**  | **Cost**   | **Peak throughput**    |
|---------------|-----------|--------------------|----------|------------|------------------------|
| Graphcore GC2 | 1216      | 31.1 TFLOPS        | 120 W \* | $8112 \*\* | 2.2 x 10^6 tracks / s  |
| i7-9750H      | 6         | 0.3 TFLOPS \*\*\*  | 45 W     | $400       | 26 x 10^6 / tracks / s |

\* Two IPUs per board, so half the board TDP.

\*\* Cost is taken as a quarter of the price of an IPU M2000, containing 4 IPUs.

\*\*\* Taken from published Geekbench 4 SGEMM score.

A rough comparison of the throughput per dollar and per watt for the current
implementations gives:

| **Device**    | **Throughput per watt** | **Throughput per dollar** |
|---------------|-------------------------|---------------------------|
| Graphcore GC2 | 0.18 x 10^5 tracks / s  | 271 tracks / s            |
| i7-9750H      | 5.77 x 10^5 tracks / s  | 65000 tracks / s          |

The CPU implementation is roughly 30 times more power efficient and 240 times
more cost efficient at present.

## IPU: Redux

Starting from scratch, the IPU implementation has now be re-worked to consolidate
all operations into a single [codelet](https://github.com/lohedges/trackr/raw/main/src/KalmanFilterCodelet.cpp),
i.e. the entire projection, filter, and smoothing stages of the Kalman filter.
The codelet can operate on an arbitrary number of tracks in serial within a
single vertex. (Within the memory limits of an IPU tile.) For efficient
parallelisation, multiple batches of tracks can be processed by separate
vertices that are mapped to different tiles, forming a _compute set_. Each
vertex receives its own copy of the required Kalman matrices, with any constant
terms in the equations pre-computed to avoid unnecessary operations. At present,
the codelet uses a set of hand-rolled, and unoptimised matrix operations. There
should be plenty of scope for improving performance by vectorising these. (Note
that the codelet is compiled with `-O3` optimisations, however.)

### Compiling

The IPU code can be compiled using:

```
make ipu
```

(Note that paths to the Poplar SDK are currently hardcoded, so would need to
be adjusted for your particular setup. Our version of the SDK is also compiled
using the old CXX11 ABI.)

### Testing

To test the IPU implementation, run:

```
./track_test_ipu
```

If successful, you should once again see:

```
Hooray, test passed!
```

Note that the test uses a single vertex, i.e. the test tracks are batched and
run in serial. We have also validated that the test passes when parallelising
over tiles using a single track per tile.

### Benchmarks

(Note that benchmarks were run on a Graphcore Collosus Mk1. The newer MK2 IPUs
have triple the memory per tile and around 20% more tiles per IPU, and have
shown up to 9x the performance in other benchmarks. As before, repeats were run
over independent, randomly generated track data.)

The IPU implementation can be benchmarked using using `trackr_benchmark_ipu`, e.g.:

```
./trackr_benchmark_ipu 256 100 100
```

The first argument is the number of tiles to run on, the second is the number
of tracks per tile, and the third is the number of repeats. When finished you
should see output like the following:

```
     256	              100	 18.846300	1.294070	    100
```

The columns are `tiles`, `tracks-per-tile`, `throughput`, `variance`, `repeats`.

Another benchmark script, `benchmark_ipu.sh`, is provided to test the
performance of the IPU implementation using a different number of tiles and
tracks per tile. To run it:

```
./benchmark_ipu.sh 200
```

Here, the value of 200 that is passed to the script is the number of tracks per
tile. The script will then measure performance for this batch size across a
range of tile numbers using 100 repeats.

You should see output like the following:

```
#  Tiles	Tracks (per tile)	Throughput	Variance	Repeats
       2	              200	  0.184652	0.000000	    100
       4	              200	  0.372116	0.000001	    100
       8	              200	  0.743286	0.000020	    100
      19	              200	  1.769149	0.000023	    100
      38	              200	  3.478308	0.003301	    100
      76	              200	  7.054743	0.000795	    100
     152	              200	 14.127182	0.001852	    100
     192	              200	 17.839186	0.003172	    100
     256	              200	 21.787933	0.200980	    100
     304	              200	 25.587167	0.435657	    100
     342	              200	 28.950420	2.432793	    100
     456	              200	 38.786932	4.290949	    100
     608	              200	 54.308519	1.226884	    100
     912	              200	 70.847941	104.611908	    100
    1216	              200	 92.133103	208.197903	    100
```

Note that we only consider data transfer of the smoothed tracks from a single
tile for the purposes of testing. Here the benchmarks are concerned with the
raw processing power and we will consider strategies for efficient data transfer
later. It also should be noted that connecting a data stream to each tile has a
memory overhead that reduces the total number of tracks that it is possible to
batch on each tile, so the throughput numbers should be taken with a grain of
salt.

The following figure shows benchmark results for a range of batch sizes run on
a single IPU, i.e parallelisation across the tiles of a single IPU only. Once
again, error bars show the standard deviation of the throughput.

![Benchmarks IPU.](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_ipu.png)

(The benchmarks code was compiled with Clang 10.0.0. Near identical performance
was found using GCC 10.1.0.)

The throughput for all batch sizes is seen to scale approximately linearly as
the number of tiles is increased, as expected for an embarrassingly parallel
(beautifully serial) workload. Increasing the number of tracks that are processed
on each tile increases the throughput, with a peak of approximately 110 million
tracks per second when processing 400 tracks on all 1216 tiles of the IPU. (In
the current setup, 400 tracks is approximately at the memory limit of each
tile.)

The following tables show updates of the CPU and IPU performance comparison.

| **Device**    | **Cores** | **32 bit FLOPS**   | **TDP**  | **Cost**   | **Peak throughput**    |
|---------------|-----------|--------------------|----------|------------|------------------------|
| Graphcore GC2 | 1216      | 31.1 TFLOPS        | 120 W    | $8112      | 110 x 10^6 tracks / s  |
| i7-9750H      | 6         | 0.3 TFLOPS         | 45 W     | $400       | 26 x 10^6 / tracks / s |

A rough comparison of the throughput per dollar and per watt for the current
implementations gives:

| **Device**    | **Throughput per watt** | **Throughput per dollar** |
|---------------|-------------------------|---------------------------|
| Graphcore GC2 | 9.17 x 10^5 tracks / s  | 13560 tracks / s          |
| i7-9750H      | 5.77 x 10^5 tracks / s  | 65000 tracks / s          |

The IPU implementation is now more energy efficient, although the CPU gives
still gives around 5 times the performance per dollar.

Note that the benchmark code was compiled with `-O3` optimisation flags. Further
tests have shown that the `-O2` optimisation level performs slightly better at
roughly 120 million tracks per scond using batches of 400 tracks on all 1216
tiles of a single IPU. In constrast, when turning off optimisations, i.e. `-O0`,
the performance drops to around 1.6 million tracks per second. Having removed
all additional compiler flags passed to `graph.addCodelets`, it appears that the
default optimisation level is `-O2`, hence the original code _should_ have been
compiled with this level of optimisation enabled.

#### Hardware threads

Each IPU tile has 6 hardware worker threads that run in a time-sliced fashion.
By scheduling multile vertices on the same tile within the same compute set, we
can use task-based parallelism to hide instruction and memory latency and
increase throughput by approximately 6x (theoretically). For the Kalman filter
code, this can be achieved by splitting the tracks assigned to each tile over
six separate vertices, i.e. we first parellelise the tracks over the tiles
of the IPU, then within each tile we further parellelise over the 6 hardware
threads. When the number of tracks that are assigned to each tile doesn't
perfectly divide between the 6 threads, then the remainder are assigned to
thread zero. For example, if 100 tracks were assigned to a tile then 20 tracks
would be assigned to thread 0 and 16 to threads 1 through 5.

The following figure shows benchmark comparison between batch sizes of 100,
200, and 400 tracks run on different numbers of tiles, using both a single
vertex per tile, and 6 vertices per tile, i.e. expoiting all hardware threads.
All data has been averaged over 100 independent repeats.

![Benchmarks IPU (threaded).](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_ipu_threaded.png)

For all batch sizes there is a marked increase in performance when using
6 hardware threads on all tiles of the IPU. For the largest batch size, 400,
the peaks throughput reachs around 430 million tracks per second, a four-fold
increase on the single-threaded performance. The throughput scaling shows
some interesting features that will be investigated further. While the
threaded performance initially increases approximately linearly, there
is then a slight drop-off in performance as the number of tiles increases.
Beyond approximately 400 tiles the performance increases linearly once more.
(In particular, the 6-thread throughput for batches of 100 tracks actually
dips below that of the single threaded code when using between approximately
250 and 550 tiles.)

#### Warm-up

All benchmarks reported previously were run using a single instance of the graph
compute engine, i.e. the IPU was run from cold. (We used 100 repeats, but each
repeat is a single run.) However, a benchmarking paper (available from the
Graphcore website
[here](https://www.graphcore.ai/hubfs/assets/pdf/Citadel%20Securities%20Technical%20Report%20-%20Dissecting%20the%20Graphcore%20IPU%20Architecture%20via%20Microbenchmarking%20Dec%202019.pdf?hsLang=en))
advocates the use of an untimed warm-up iteration of any benchmark prior to
timing. This is used to exlude any warm-up overheads from the results,
ensuring that steady-state measurements are achieved. The following figure
shows a comparison between the _cold_ threaded benchmarks shown previously
and those run once the system was _warm_.

![Benchmarks IPU (threaded and warm).](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_ipu_threaded_warm.png)

Following the warm-up run, the scaling of throughput versus tile number is near
perfectly linear. In addition, run-to-run fluctuations are massively diminished,
as seen in the size of the error bars. For the largest batch size of 400 tracks,
the peak throughput when running on all 1216 tiles is now approximately 590
million tracks per second.

#### Profiling

The Poplar SDK comes with some excellent
[profiling tools](https://github.com/graphcore/examples/tree/master/tutorials/poplar/tut4_profiling).
To enable profiling during benchmarking add `true` at the end of the command-line,
e.g.:

```
./trackr_benchmark_ipu 256 100 100 true
```

This will write a profiling report to `profile.txt` in the current directory.
Sensible execution profiling options have been set, but can be adjusted by
setting appopriate entries in the `POPLAR_ENGINE_OPTIONS` environment variable,
e.g.:

```
 POPLAR_ENGINE_OPTIONS='{"debug.instrumentCompute":"true", "debug.computeInstrumentationLevel":"tile"}'
```

(Here `debug.computeInstrumentationLevel` can be one of `device`, `ipu`, or `tile`.)

#### Vectorisation (sort of)

Following advice [here](https://www.graphcore.ai/hubfs/assets/pdf/Citadel%20Securities%20Technical%20Report%20-%20Dissecting%20the%20Graphcore%20IPU%20Architecture%20via%20Microbenchmarking%20Dec%202019.pdf?hsLang=en)
and [here](https://github.com/thorbenlouw/BabelStream/blob/master/PoplarKernels.cpp)
we have tried several approaches to encourage the `popc` compiler to
vectorise the various matrix operations within the Kalman filter codelet.
These approaches make use of aligned pointers and type conversion, which is used
to access the vertex fields via the vector type `float2`. This enables
(theoretically) the compiler to emit 64-bit instructions for 32-bit values,
i.e. allowing two `float` elements to be operated for a single `float2`, e.g.:

```cpp
class MyVertex : public poplar::Vertex
{
    // Some read/writeable field, request 8-byte alignment.
    poplar::InputOut<Vector<float, poplar::VectorLayout::SPAN, 8>> in;

    ...

    bool compute()
    {
        // Reinterpret bit pattern as float2.
        float2 *f2in = reinterpret_cast<float2 *>(&in[0]);

        ...
```

In addition, the `popc` compiler also supports the use of the `__restrict` type
qualifier, meaning that we can specify that the pointer is not aliased, which
should help the compiler to vectorise any loops over any arrays accessed via
the pointer. The compiler also allows uses of `#pragma unroll` directives to
indicate the appropriate unroll factor for any loop.

Combining the above approaches the codelet helper functions were re-written to
cast all tensor rows (track hits) to `float2 *` and pass these to sub
functions to perform per-row operations (copies, additions, subtractions) in
an unrolled fashion. From the benchmarking paper linked to above, it would
appear that an unroll factor of 8 would be most appropriate. With this in
mind, the bechmark code was updated to require track numbers that are multiples
of 6 and 8 only, i.e. divisible over the hardware threads on each IPU and
possible to be unrolled by a factor of 8.

The multiplication function poses a problem for vectorisation since both
tensors have the same storage order. Since its not possible to easily
transpose one of them on the vertex, we apply a manual optimisation by
completely unrolling the inner-most loop, which is of fixed size, i.e.:

```cpp

// All rows in 'in1' have same number of columns.
int size = in1[0].size();

for (int i=0; i<4; ++i)
{
    for (int j=0; j<size; ++j)
    {
        out[i][j] = in0[i][0] * in1[0][j];
        for (int k=1; k<4; ++k)
        {
            out[i][j] += in0[i][k] * in1[k][j];
        }
    }
}
```

becomes

```cpp

int size = in1[0].size();

for (int i=0; i<4; ++i)
{
    for (int j=0; j<size; ++j)
    {
        out[i][j] = in0[i][0] * in1[0][j]
                  + in0[i][1] * in1[1][j]
                  + in0[i][2] * in1[2][j]
                  + in0[i][3] * in1[3][j];
    }
}
```

The following plot shows a benchmark comparison of 400 tracks per tile following
a warmup run, to 408 tracks per tile run in the _vectorised_ fashion.

![Benchmarks IPU (threaded, warm, vectorised).](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_ipu_threaded_warm_vector.png)

When running on all 1216 tiles of the IPU, the peak throughput now reaches
approximately 1.15 billion tracks per second.

Breaking down the peformance contribution of the different vectorisation tricks
gives some surprising results.

* The use of `__restrict` makes no difference whatsover.
* Loop unrolling generally slows the code (slightly), i.e. the optimum unroll
factor is 1. However, manually unrolling the inner-most multiplication loop
does produce a big performance gain.
* Additional performance gains seem to come from type conversion and performing
the matrix operations on a per-row basis. Reinterpreting as `float2 *` has a
marginal benefit over `float *`. (Or performing no conversion and just operating
on the original tensors row-by row.)
* Using an optimisation level of `-O3` now gives marginally better performance
than `-O2`.

Updated IPU/CPU performance comparison:

| **Device**    | **Throughput per watt** | **Throughput per dollar** |
|---------------|-------------------------|---------------------------|
| Graphcore GC2 | 95.83 x 10^5 tracks / s | 14.18 x 10^4 tracks / s   |
| i7-9750H      |  5.77 x 10^5 tracks / s |  6.50 x 10^4 tracks / s   |

#### Multi-IPU benchmarks

It is trivial to run the benchmark code on multi-IPU devices by simply selecting
the required number of IPUs when connecting to available hardware devices, e.g.:

```cpp
int num_ipus = 2;
auto dm = poplar::DeviceManager::createDeviceManager();
auto hwDevices = dm.getDevices(poplar::TargetType::IPU, num_ipus);
```

This returns a _virtual_ IPU, containing a multiple of the number of tiles
per IPU, e.g. 2432 tiles when selecting 2 IPUs. Since the code is embarassingly
parallel, scaling to more IPUs should produce a near-linear increase in
performance. Benchmarking 408 tracks per tile using all 2432 tiles on 2 IPUs
is found to give a throughput of approximately 2.3 billion tracks per second,
i.e. twice the single IPU throughput.

#### Manual optimisations

As with any piece of numerical software, there is scope for performance gains
by examining the equations that are modelled to remove redundant floating
point operations. In this case, the Kalman matrices contain many zero terms
that lead to wasted performance when performing matrix multiplications, sums,
copies, inverses, etc. For example, in the absence of noise (as considered
here), the multiplicand in any matrix multiplication always has zero entries
in elements (0,2), (0,3), (1,2), (1,3), (2,0), and (2,1). As such, the
multiplication code shown above can be further simplified to:

```cpp
int size = in1[0].size();

// Using two separate loops appears to produce the best optimisation.

for (int i=0; i<2; ++i)
{
    for (int j=0; j<size; ++j)
    {
        out[i][j] = in0[i][0] * in1[0][j]
                  + in0[i][1] * in1[1][j];
    }
}
for (int i=2; i<4; ++i)
{
    for (int j=0; j<size; ++j)
    {
        out[i][j] = in0[i][2] * in1[2][j]
                  + in0[i][3] * in1[3][j];
    }
}
```

(Note that other terms in the multiplicand may be zero. The above is the
easiest generalisation that avoids the need for separate code to handle
different matrices.)

The following plot shows benchmarks after applying these simple optimisations.

![Benchmarks IPU (threaded, warm, vectorised, manual optimisations).](https://github.com/lohedges/trackr/raw/main/benchmarks/benchmark_ipu_threaded_warm_vector_manual.png)

The performance when running batches of 408 tracks on all 1216 tiles has now
reached approximately 1.8 billion tracks per second.

While the manual optimisations described here might not be appropriate for
real world code, where the addition of noise would mean that the number of
zero matrix terms is fewer, it highlights the benefit of thinking about the
equations and reducing redundancy wherever possible. Here we achieved roughly
40% increase in throughput for little effort.
