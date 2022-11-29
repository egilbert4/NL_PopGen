#!/bin/bash

###############################################################################
## Set up
sequence=nfl-gnomad-merged

######################################################################
# stage 1: combine datasets together
## grab the raw data form both NL and gnomad v3.1
# set up the locations
nl_data=vcf38-20200918-gold.vcf.bgz
nl_keep=NL.founder-indigenous_ancestry.idfile
nl_remove=NL-ancestry.related.remove

gnomad_data=gnomad_1kg_hgdp_filtered_vars
gnomad_remove=gnomad.related.txt

######################################################################
# find common snps and merge NFL and GNOMAD
plink2 --vcf $nl_data --autosome --snps-only --max-alleles 2 --make-bed --out nl.raw
plink2 --bfile $gnomad_data --autosome --snps-only --max-alleles 2 --make-bed --out gnomad.raw

cat nl.raw.bim | awk '{print $1, $1":"$4, $3,$4,$5,$6}' > tmp-nl.raw.bim
mv tmp-nl.raw.bim nl.raw.bim

cat gnomad.raw.bim | awk '{print $1, $1":"$4, $3,$4,$5,$6}' > tmp-gnomad.raw.bim
mv tmp-gnomad.raw.bim gnomad.raw.bim

# remove the strand ambigous SNPs
cat nl.raw.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > nl.snps
cat gnomad.raw.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > gnomad.snps

join gnomad.snps nl.snps > common.snps

# extract these common snps
plink2 --bfile nl.raw --extract common.snps --make-bed --out nl.common
cat nl.common.fam | awk '{print "NL",$2,$3,$4,$5,$6}' > tmp.fam
mv tmp.fam nl.common.fam
plink2 --bfile nl.common --keep $nl_keep --remove $nl_remove --make-bed --out nl.common-filtered

plink2 --bfile gnomad.raw --extract common.snps --make-bed --out gnomad.common
cat gnomad.common.fam | awk '{print "GNOMAD",$2,$3,$4,$5,$6}' > tmp.fam
mv tmp.fam gnomad.common.fam
plink2 --bfile gnomad.common --remove $gnomad_remove --make-bed --out gnomad.common-filtered

# merge these together
plink --bfile nl.common-filtered --bmerge gnomad.common-filtered --make-bed --out $sequence.merged



