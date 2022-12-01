#!/bin/bash

###############################################################################
## run IBDNe on the Irish-British IBD-clusters

##################
# variables

# software
browning_dir=beagle/
ibdne=$browning_dir/ibdne.23Apr20.ae9.jar

build=37
genmap_dir=genetic_maps/build$build
cat $genmap_dir/plink.chr{1..22}.GRCh$build.map > $tmp_dir/tmp.map
mincm=4

##################
# run the IBD Ne analysis on the second level Irish-British clusters
# [1:1, 1:2, 1:3, 1:4, 2:1, 2:2, 2:3, 2:4, 3:1, 3:2, 3:3, 3:4]
# or ["CS England", "Orkney", "Wales", "Eng-Sco-Ire", "SW Ireland", "SE Ireland", "NW Ireland"]

# the .ibd file contains only segments shared between individuals placed in the cluster

for clust in "1:1" "1:2" "1:3" "1:4" "2:1" "2:2" "2:3" "2:4" "3:1" "3:2" "3:3" "3:4"
  do
  zcat NL-IreBrit.target-ref-clusters.${clust}.ibd.gz | \
        java -jar $ibdne \
        map=$tmp_dir/tmp.map \
        nthreads=6 \
        out=ne_data/NFL-IreBrit.target-ref-clusters.${clust}.${mincm}cM.ibdne-out \
        mincm=${mincm} \
        nboots=80
  done
  
  %sh 

###############################################################################
## Run IBDNe on the NFL fs-clusters

##################
# run the IBD Ne analysis

for clust in {1..22}
  do
  zcat NL-IreBrit.nfl-fs-clusters.${clust}.ibd.gz | \
        java -jar $ibdne \
        map=$tmp_dir/tmp.map \
        nthreads=6 \
        out=ne_data/NFL-IreBrit.nfl-fs-clusters.${clust}.${mincm}cM.ibdne-out \
        mincm=${mincm} \
        nboots=80
  done
