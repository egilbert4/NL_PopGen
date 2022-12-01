#!/bin/bash

###############################################################################
## Extract and clean raw genotype data 
sequence=NL-IreBrit

#######################
## grab the raw data form both NFL and Ire-Brit
# set up the locations
irebrit_data=ire-brit_data/ATP.raw
irebrit_keep=ire-brit_data/ATP.unrelated.txt

nl_data=vcf37-20200918-gold.vcf.bgz
nl_keep=NL.founder-indigenous_ancestry.idfile
nl_remove=NL-ancestry.related.remove

# extract the common SNPs for each dataset
plink2 --vcf $nl_data --autosome --snps-only --max-alleles 2 --sort-vars --make-pgen --out nl.raw
plink2 --pfile nl.raw --make-bed --out nl.raw
cat nl.raw.fam | awk '{print "NFL",$2,$3,$4,$5,$6}' > tmp.fam; mv tmp.fam nl.raw.fam

cat nl.raw.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > nl.snps
cat $irebrit_data.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > ire-brit.snps

join nl.snps ire-brit.snps > common.snps

#######################
# extract common SNPs and
plink --bfile nl.raw --extract common.snps --keep $nfl_keep --remove $nfl_remove --make-bed --out nl.common
plink --bfile $irebrit_data --extract common.snps --keep $irebrit_keep --make-bed --out ire-brit.common

plink --bfile ire-brit.common --bmerge nfl.common --make-bed --out $sequence
plink --bfile ire-brit.common --flip $sequence-merge.missnp --make-bed --out ire-brit.flipped
plink --bfile ire-brit.flipped --bmerge nl.common --make-bed --out $sequence

#######################
# clean
plink --bfile $sequence --geno 0.05 --maf 0.01 --hwe 1e-6 --make-bed --out $sequence.filtered
plink --bfile $sequence.filtered --mind 0.05 --make-bed --out $sequence.filtered-cleaned

# LD
plink2 --bfile $sequence.filtered-cleaned --indep-pairwise 1000 50 0.2 --out $sequence.LD

# perform PCA
plink --bfile $sequence.filtered-cleaned --extract $sequence.LD.prune.in --pca 20 --out $sequence.pruned-pca

# roh
plink --bfile $sequence.filtered-cleaned \
  --homozyg \
  --homozyg-window-snp 50 \
  --homozyg-snp 50 \
  --homozyg-kb 1500 \
  --homozyg-gap 1000 \
  --homozyg-density 50 \
  --homozyg-window-missing 5 \
  --homozyg-window-het 1 \
  --out $sequence.roh

###############################################################################
## Phase data

#######################
# phasing the Ire-Brit and NL
build=37
for chrom in {1..22}
  do
  # the genetic map
  genmap=genetic_maps/chr$chrom.b$build.gmap.gz

  # create fam file
  plink2 --bfile $sequence.filtered-cleaned\
    --chr $chrom \
    --make-just-fam --out phasing/input_data/$sequence.chr$chrom
  
  # extract the genotypes
  plink2 --bfile $sequence.filtered-cleaned \
    --chr $chrom \
    --export vcf bgz id-paste=iid \
    --out phasing/input_data/$sequence.chr$chrom
  
  # create index
  tabix -p vcf phasing/input_data/$sequence.chr$chrom.vcf.gz
  
  # phase
  shapeit4 --input phasing/input_data/$sequence.chr$chrom.vcf.gz \
    --map $genmap \
    --region $chrom \
    --thread 6 \
    --output phasing/output_data/$sequence.chr$chrom.phased.vcf.gz
  done
  
%sh

###############################################################################
## detect IBD segment in NL-Ire-Brit

#######################
# software
browning_dir=beagle/
beagle=$browning_dir/beagle.21Apr21.304.jar
refinedibd=$browning_dir/refined-ibd.17Jan20.102.jar
mergeibd=$browning_dir/merge-ibd-segments.17Jan20.102.jar

#######################
# ibd detection using refinedIBD
for chrom in {1..22}
  do
  # the genetic map
  genmap=genetic_maps/build$build/plink.chr$chrom.GRCh$build.map
    
  # refinedIBD
  java -Xmx4g -jar $refinedibd \
    gt=phasing/output_data/$sequence.chr$chrom.phased.vcf.gz \
    map=$genmap \
    out=ibd_data/$sequence.chr$chrom.refinedIBD \
    chrom=$chrom \
    length=1 nthreads=8
  
  # merge-IBD
  zcat ibd_data/$sequence.chr$chrom.refinedIBD.ibd.gz | \
    java -Xmx4g -jar $mergeibd phasing/output_data/$sequence.chr$chrom.phased.vcf.gz \
    $genmap 0.6 1 > ibd_data/$sequence.chr$chrom.mergeibd.ibd
  gzip -f ibd_data/$sequence.chr$chrom.mergeibd.ibd
  
  done
