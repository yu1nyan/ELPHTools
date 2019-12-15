#!/bin/bash
export run=$1
export sub=0
export nfile=$2

./run_calib hodo $run $sub $nfile 20 70 830 830 png
./run_calib proto1s $run $sub $nfile 20 70 825 825 png
./run_calib proto2s $run $sub $nfile 20 70 845 830 png
