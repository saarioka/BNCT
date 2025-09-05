#!/bin/bash

set -e

source ../../env.sh

for ((i=2; i<=80; i+=2))
do
    export RUN_ID="${i}mm"
    export MODERATOR_THICKNESS="$i"  # cm
    cmake ..
    make -j20
    echo "Running simulation with moderator thickness ${MODERATOR_THICKNESS} mm"
    ./exampleB2b run.mac > "output_${RUN_ID}.log"
done
