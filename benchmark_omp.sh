#!/usr/bin/env bash

if [ -z "$1" ]
then
    num_tracks=500
else
    num_tracks=$1
fi

num_batches=$(seq 30)
repeats=10000

echo -e "# Tracks\tThroughput\tVariance\tRepeats\tBatches"
for batches in ${num_batches[@]}; do
    ./trackr_benchmark_omp $num_tracks 100 $batches
done
