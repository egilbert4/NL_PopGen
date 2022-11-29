#!/bin/bash

######################################################################
## Clean the merged data

sequence=nl-gnomad-merged

# update the ids
plink2 --bfile $sequence.merged \
  --keep nl-gnomad-merged.keep.idfile \
  --make-bed --out $sequence.keep

plink2 --bfile $sequence.keep \
  --update-ids nl-gnomad-merged.update.idfile \
  --make-bed --out $sequence.updated

# cleaning and filtering SNPs and individuals
plink2 --bfile $sequence.updated --geno 0.05 --maf 0.01 --hwe 1e-9 midp keep-fewhet --make-bed --out $sequence.filtered
plink2 --bfile $sequence.filtered --mind 0.05 --make-bed --out $sequence.filtered-cleaned

######################################################################
## LD prune data and run PCA

# LD
plink2 --bfile $sequence.filtered-cleaned --indep-pairwise 1000 50 0.2 --out $sequence..LD

# pca unprojected
plink --bfile $sequence.filtered-cleaned \
    --extract $sequence..LD.prune.in \
    --pca --out $sequence..unproject-pca

# pca projected onto GNOMAD refs
cat $sequence.filtered-cleaned.fam | awk '{print $1, $2, $1}' > $sequence.clust
  
plink --bfile $sequence.filtered-cleaned \
    --extract $sequence.LD.prune.in \
    --mac 1 \
    --within $sequence.clust --pca-cluster-names GNOMAD \
    --pca --out $sequence.projected-pca
  
######################################################################
## Cell 6 - perform supervised admixture with world-wide references

# perform the K=5 analysis of the 1000 Genomes References
cat $sequence.filtered-cleaned.fam | grep "KGP" | awk '{print $1,$2}' > KGP.idfile

plink --bfile $sequence.filtered-cleaned \
    --extract $sequence..LD.prune.in \
    --keep KGP.idfile \
    --make-bed --out KGP3.admixture

k=5
rep=1
admixture -s time -j4 KGP3.admixture.bed $k

mv KGP3.admixture.$k.P KGP3.admixture.K$k.R$rep.P
gzip -f KGP3.admixture.K$k.R$rep.P

cat KGP3.admixture.fam | awk '{print $1,$2}' > tmp
paste -d " " tmp KGP3.admixture.$k.Q > KGP3.admixture.K$k.R$rep.Q

######################################################################
## Supervised analysis on NL using .P data from KGP3

# extract data
cat $sequence.filtered-cleaned.fam |
  grep "fs_" | awk '{print $1,$2}' > NL.idfile
  
plink --bfile $sequence.filtered-cleaned \
  --extract $sequence.LD.prune.in \
  --keep NL.idfile \
  --make-bed --out NL.admixture

# run supervised admixture for NL
zcat KGP3.admixture.K5.R1.P > NL.admixture.5.P.in
admixture -s time -j4 -P NL.admixture.bed 5

rm NL.admixture.5.P.in
mv NFL.admixture.5.P NL.admixture.K5.P
gzip -f NL.admixture.K5.P
cat NL.admixture.fam | awk '{print $1,$2}' > tmp
paste -d " " tmp NL.admixture.5.Q > NL.admixture.K5.Q
