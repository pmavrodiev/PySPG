#!/bin/bash
# $1 -i
# $2 input_file

VARIABLES=`spg-get_parameters.py -i $2 -p store_dynamics_filename,lnTruth,t,N`
R_SCRIPT='/home/pmavrodiev/run/rank_detailed_investigations3.R'

mkdir -p /local/scratch/pmavrodiev
ctx-rank.p4 -i $2

R --vanilla --slave --args $VARIABLES < $R_SCRIPT


