#!/bin/bash

###############################################################################
# set up the ASCEND parameters

# https://github.com/sunyatin/ASCEND
# https://github.com/sunyatin/ASCEND/blob/master/docs/present_DNA.pdf

for clust in fs_5 fs_5  fs_6 fs_7 fs_8 fs_9 fs_10 fs_11 fs_12 fs_13 fs_14 fs_15 fs_16 fs_17 fs_18 fs_19 fs_20 fs_21 fs_22 IBD_1_1 IBD_1_2 IBD_1_3 IBD_1_4 IBD_2_1 IBD_2_2
  do
  echo -e "\n####\n Cluster: ${clust}"
  
  echo -e "genotypename: $sequence.geno
snpname: $sequence.snp
indivname: $sequence.ascend-label.ind
targetpop: $clust
outpop: YRI
outputprefix: $sequence.${clust}.ascend
binsize: 0.001
mindis: 0.001
maxdis: 0.3
maxpropsharingmissing: 1
minmaf: 0
haploid: NO
dopseudohaploid: NO
morgans: NO
onlyfit: NO
usefft: YES
qbins: 100
seed: 1234
blocksizename: None" > $sequence.${clust}.ascend.par
  
  python3 ASCEND.py -p $sequence.${clust}.ascend.par
  done
