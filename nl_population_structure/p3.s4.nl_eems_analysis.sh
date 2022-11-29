#!/bin/bash

# move to the EEMS cluster and convert the plink data into the diffs file

# variables
sequence=NL-ancestry

# use bed2diffs to convert
bed2diffs_v1 --bfile $sequence.geog-filtered --nthreads 4
