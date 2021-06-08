#!/usr/bin/env bash

if [ -z $1 ]; then
    tracks=1
else
    tracks=$1
fi

num_tiles=(2 4 8 19 38 76 152 192 256 304 342 456 608 912 1216)
repeats=10

echo -e "#  Tiles\tTracks (per tile)\tThroughput\tVariance\tRepeats"
for num in ${num_tiles[@]}; do
    ./trackr_benchmark_ipu $num $tracks $repeats
done
