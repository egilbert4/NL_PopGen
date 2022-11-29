#!/bin/bash

###############################################################################
## Convert the plink data into the diffs file
# variables
sequence=NL-ancestry

# use bed2diffs to convert
bed2diffs_v1 --bfile $sequence.geog-filtered --nthreads 4

###############################################################################
# perform the initial EEMS analysis

# variables
ndemes=300

mkdir intial/

# create param files
for chain in {1..10}
  do
  echo -e "datapath = $sequence.geog-filtered
  mcmcpath = initial/$sequence.EEMS-initial-nDemes${ndemes}.chain${chain}
  nIndiv = 739
  nSites = 685221
  nDemes = ${ndemes}
  diploid = true
  numMCMCIter = 5000000
  numBurnIter = 3000000
  numThinIter = 9999" > $sequence.EEMS-initial-nDemes${ndemes}.chain${chain}.param

  # run EEMS
  runeems_snps --params $sequence.EEMS-initial-nDemes${ndemes}.chain${chain}.param
  done

###############################################################################
# perform the extension EEMS analysis, make sure this is using the EEMS cluster
# variables
mkdir extension/

ndemes=300
start_chain=5

# create param file
for chain in {1..10}
  do
  echo -e "datapath = $sequence.geog-filtered
  mcmcpath = extension/$sequence.EEMS-extension-nDemes${ndemes}.chain${chain}
  prevpath = initial/$sequence.EEMS-initial-nDemes${ndemes}.chain${start_chain}
  nIndiv = 739
  nSites = 685221
  nDemes = ${ndemes}
  diploid = true
  numMCMCIter = 3000000
  numBurnIter = 2000000
  numThinIter = 9999" > $sequence.EEMS-extension-nDemes${ndemes}.chain${chain}.param

  # run EEMS
  runeems_snps --params $sequence.EEMS-extension-nDemes${ndemes}.chain${chain}.param
  done
