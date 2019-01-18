#!/bin/bash

export run=$1
export sub=2
export nfile=$2

./run_calib hodo $run $sub $nfile 20 70 850 850 png
./run_calib proto1s $run $sub $nfile 20 70 780 780 png
./run_calib proto2s $run $sub $nfile 20 70 800 780 png
