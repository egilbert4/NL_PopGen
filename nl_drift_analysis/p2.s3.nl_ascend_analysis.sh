#!/bin/bash

###############################################################################
## Cell 4 - clean the merged data
sequence=nfl-gnomad-merged

# update the ids
plink2 --bfile $sequence.merged \
  --keep nl-gnomad-merged.keep.idfile \
  --make-bed --out $sequence.keep

plink2 --bfile $sequence.keep \
  --update-ids nl-gnomad-merged.update.idfile \
  --make-bed --out $sequence.updated

rm $workdir/$sequence.keep.*

# cleaning
plink2 --bfile $sequence.updated --geno 0.05 --maf 0.01 --hwe 1e-9 midp keep-fewhet --make-bed --out $sequence.filtered
plink2 --bfile $sequence.filtered --mind 0.05 --make-bed --out $sequence.filtered-cleaned

