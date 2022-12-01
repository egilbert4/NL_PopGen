#!/bin/bash

###############################################################################
## Extract and clean raw genotype data

# variables
sequence=NL-ancestry

raw_data=vcf38-20200918-gold.vcf.bgz
ind_file=NL.founder-indigenous_ancestry.idfile

# extract autosomal bi-allelic SNPs
plink2 --vcf $raw_data --autosome --snps-only --max-alleles 2 --sort-vars --make-pgen --out $sequence.raw

# refomat ID information
plink2 --pfile $sequence.raw --make-bed --out $sequence.raw
cat $sequence.raw.fam | awk '{print "NL", $2,$3,$4,$5,$6}' > tmp.fam
mv tmp.fam $sequence.raw.fam

# exclude strand-ambiguous SNPS
cat $sequence.raw.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > non.at-gc.snps
 
# filter SNPs
plink --bfile $sequence.raw --extract non.at-gc.snps --geno 0.05 --maf 0.01 --hwe 1e-6  --make-bed --out $sequence.filtered

# clean up any poorly genotyped individuals
plink --bfile $sequence.filtered --keep $idfile --mind 0.05 --make-bed --out $sequence.filtered-cleaned

###############################################################################
## Quantify relatedness
king -b $sequence.filtered-cleaned.bed --cpus 8 --related --degree 3 --prefix $sequence.related

%sh 

###############################################################################
## Perform PLINK PCA
removefile=/dbfs/mnt/shared-data/ed_results/NFL-ancestry.related.remove

# prune for LD
plink2 --bfile $sequence.filtered-cleaned --remove $removefile --indep-pairwise 1000 50 0.2 --out $sequence.ld

# run pca
plink2 --bfile $sequence.filtered-cleaned --remove $removefile --extract $sequence.ld.prune.in --make-rel square --pca 20 --out $sequence.pruned-pca

%sh

###############################################################################
## Run PLINK ROH calling
removefile=/dbfs/mnt/shared-data/ed_results/NFL-ancestry.related.remove

plink --bfile $sequence.filtered-cleaned \
  --remove $removefile \
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
## Perform Phasing

# variables
removefile=/dbfs/mnt/shared-data/ed_results/NFL-ancestry.related.remove

mkdir $workdir/phasing
mkdir $workdir/phasing/input_data
mkdir $workdir/phasing/output_data

build=38

# extract the initial genotypes
for chrom in {1..22}
  do
  
  # the genetic map
  genmap=shapeit/genetic_maps/chr$chrom.b$build.gmap.gz

  # create fam file
  plink2 --bfile $sequence.filtered-cleaned --remove $removefile --chr $chrom --make-just-fam --out phasing/input_data/$sequence.chr$chrom
  
  # extract the genotypes in vcf format
  plink2 --bfile $sequence.filtered-cleaned \
    --remove $removefile \
    --chr $chrom \
    --export vcf bgz id-paste=iid \
    --out phasing/input_data/$sequence.chr$chrom
  
  # index the vcf
  tabix -p vcf phasing/input_data/$sequence.chr$chrom.vcf.gz
  
  # phase with SHAPEIT
  shapeit4 --input phasing/input_data/$sequence.chr$chrom.vcf.gz \
    --map $genmap \
    --region $chrom \
    --thread 6 \
    --output phasing/output_data/$sequence.chr$chrom.phased.vcf.gz
  
  done

###############################################################################
## Perform IBD detection using refinedIBD
mkdir $workdir/ibd_data

# software
browning_dir=/usr/local/share/beagle/
beagle=$browning_dir/beagle.21Apr21.304.jar
refinedibd=$browning_dir/refined-ibd.17Jan20.102.jar
mergeibd=$browning_dir/merge-ibd-segments.17Jan20.102.jar

# sequence
build=38

#######################
# ibd detection using refinedIBD
for chrom in {1..22}
  do
  # the genetic map
  genmap=$browning_dir/genetic_maps/build$build/plink.chr$chrom.GRCh$build.map
    
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

###############################################################################
## Conversion of the phased data to ChromoPainter format 

# variables
mkdir $workdir/fineSTRUCTURE
mkdir $workdir/fineSTRUCTURE/data_files

# software
impute2chromopainter=fs_4.1.1/impute2chromopainter.pl
convertrec=fs_4.1.1/convertrecfile.pl

#######################
# covert to haps/sample and then the cp
for chrom in {1..22}
  do
  ## the vcf to haps
  # the genetic map
  genmap=shapeit/genetic_maps/chr$chrom.b$build.gmap.gz
  
  # converting
  input=phasing/output_data/$sequence.chr$chrom.phased.vcf.gz
  output=phasing/output_data/$sequence.chr$chrom.phased
    
  bcftools convert $input --hapsample $output
  mv $output.sample $output.samples
  
  ## the cp
  # the input/output filenames   
  inhap=phasing/output_data/$sequence.chr$chrom.phased.hap
  outhap=fineSTRUCTURE/data_files/$sequence.chr$chrom.phased
  outhaplog=fineSTRUCTURE/data_files/$sequence.chr$chrom.hap.log

  inrec=$outhap.phase
  outrec=fineSTRUCTURE/data_files/$sequence.chr$chrom.phased.recomrates
  outreclog=fineSTRUCTURE/data_files/$sequence.chr$chrom.phased.recomrates.log

  # convert the haplotypes
  gunzip $inhap.gz
  perl $impute2chromopainter $inhap $outhap | tee $outhaplog
  gzip -f $inhap
    
  # convert the recombination files
  zcat $genmap | tail -n +2 | awk -v CHROM=$chrom '{print CHROM,$1,$2,$3}' > fineSTRUCTURE/data_files/temp.txt
  perl $convertrec -M hapmap $inrec fineSTRUCTURE/data_files/temp.txt $outrec | tee $outreclog
  rm fineSTRUCTURE/data_files/temp.txt
  done

###############################################################################
## Running fineSTRUCTURE - stage 1

#######################
# newly created ID files using the format of fs
full_idfile=$fsdir/$sequence.idfile
subset_idfile=$fsdir/$sequence.subset.idfile

# fs parameters
hpcmode=1
indsperproc=1
s1indfrac=1

# change directory
cd fineSTRUCTURE

# set up the stage 1 analysis for the chromosomes
for chrom in 3 9 16 21
  do
  fs ${sequence}_stage1_chr${chrom}.cp \
    -idfile $subset_idfile \
    -phasefiles data_files/$sequence.chr${chrom}.phased.phase \
    -recombfiles data_files/$sequence.chr${chrom}.phased.recomrates \
    -indsperproc $indsperproc \
    -s1indfrac $s1indfrac \
    -exec fs \
    -n \
    -hpc $hpcmode \
    -go
  done
  
# combine the commands
cat ${sequence}_stage1_chr3/commandfiles/commandfile1.txt \
  ${sequence}_stage1_chr9/commandfiles/commandfile1.txt \
  ${sequence}_stage1_chr16/commandfiles/commandfile1.txt \
  ${sequence}_stage1_chr21/commandfiles/commandfile1.txt > stage1.commands.txt

# run the commands
cat stage1.commands.txt | parallel

# then update the cp files
for chrom in 3 9 16 21
  do
  fs ${sequence}_stage1_chr${chrom}.cp -go
  done

# find the information you need for weighted averages
for chrom in 3 9 16 21; do cat ${sequence}_stage1_chr${chrom}.cp | grep "Neinf"; done
echo ""
for chrom in 3 9 16 21; do cat ${sequence}_stage1_chr${chrom}.cp | grep "muinf"; done
echo ""
for chrom in 3 9 16 21; do head -n2 data_files/$sequence.chr${chrom}.phased.phase | tail -n1; done

# then clean up the tmp files
for chrom in 3 9 16 21
  do
  echo $chrom
  rm $fsdir/${sequence}_stage1_chr${chrom}/stage1/*
  done

###############################################################################
## Runnning fineSTRUCTURE - stage 2

#######################
# set up the fs parameters
hpcmode=1
indsperproc=1

s3iters=4000000 # half each assigned to burnin and sampling
s3samples=500
s4iters=200000

 set up the project by starting the stage 1 - at least to make1
fs ${sequence}.cp \
  -s34args:"-T 1 -k 2" \
  -idfile $full_idfile \
  -phasefiles $fsdir/data_files/$sequence.chr{1..22}.phased.phase \
  -recombfiles $fsdir/data_files/$sequence.chr{1..22}.phased.recomrates \
  -indsperproc $indsperproc \
  -s2chunksperregion 50 \
  -s3iters $s3iters \
  -maxretained $s3samples \
  -s4iters $s4iters \
  -exec fs \
  -n \
  -hpc $hpcmode \
  -go

# using the estimations of Ne and mu from stage 1
fs ${sequence}.cp \
    -Neinf 184.4716 \
    -muinf 0.0002307115 \
    -go
    
echo "#### Running stage 2 commands...."
cat ${sequence}/commandfiles/commandfile2.txt | parallel

echo "#### Combining stage 2 results...."
fs ${sequence}.cp -go

echo "#### Cleaning Up stage 2 results...."
for chrom in {1..22}
  do
  echo $chrom
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.log
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.chunkcounts.out
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.regionchunkcounts.out
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.regionsquaredchunkcounts.out
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.chunklengths.out
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.mutationprobs.out
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.EMprobs.out
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.prop.out
  rm ${sequence}/stage2/${sequence}_tmp_mainrun.linked_file${chrom}_ind*.samples.out.gz
  done

###############################################################################
### Process fineSTRUCTURE stage 3
chmod -u+x ${sequence}/commandfiles/commandfile3.txt

./${sequence}/commandfiles/commandfile3.txt

###############################################################################
### CELL 12 - Finish stage 3 and run stage 4
fs ${sequence}.cp -go

cat ${sequence}/commandfiles/commandfile4.txt | parallel

fs ${sequence}.cp -go
