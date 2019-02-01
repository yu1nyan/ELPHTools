#!/bin/bash

export run=$1
export sub=6
export nfile=$2

./run_calib hodo $run $sub $nfile 20 70 830 830 png
./run_calib proto2s $run $sub $nfile 20 70 850 850 png
