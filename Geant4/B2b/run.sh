#!/bin/bash

set -e

source ../../env.sh

for ((i=0; i<=8; i+=2))
do
    export RUN_ID="${i}mm"
    export MODERATOR_THICKNESS="$i"  # cm
    cmake ..
    make -j20
    ./exampleB2b run.mac > "output_${RUN_ID}.log"
done
