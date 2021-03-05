#!/usr/bin/env bash

num_tracks=(2 4 8 16 32 64 128 256 512 1024)
repeats=100000

echo -e "# Tracks\tThroughput\tVariance\tRepeats"
for num in ${num_tracks[@]}; do
    ./trackr_benchmark $num $repeats
done
