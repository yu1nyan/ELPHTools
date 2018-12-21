#!/bin/bash

export run=$1
export sub=$2
export nfile=$3

./run_calib hodo $run $sub $nfile 20 70 830 830 png
./run_calib proto1s $run $sub $nfile 20 70 830 830 png
./run_calib proto2s $run $sub $nfile 20 70 855 840 png
