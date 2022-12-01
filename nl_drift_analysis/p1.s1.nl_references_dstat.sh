## script for

###############################################################################
## copy over the id files for nfl, ire-brit, and gnomad
cp ${datadir}/NL-ancestry.K22-clust-membership.csv NL-ancestry.K22-clust-membership.csv
cp ${datadir}/gnomad_idfile.csv gnomad_idfile.csv
cp ${datadir}/NL-IreBrit.IBD-clusters-nf_ref.tsv NL-IreBrit.IBD-clusters-nf_ref.tsv

###############################################################################
## Extract and merge with NL-SB data with the gnomad data
sequence=nfl-gnomad-merged

nl_data=vcf38-20200918-gold.vcf.bgz
nl_keep=NL.founder-indigenous_ancestry.idfile
nl_remove=NL-ancestry.related.remove

gnomad_data=/dbfs/mnt/shared-data/gnomad_1kg_hgdp_filtered_vars
gnomad_remove=$workdir/gnomad.related.txt

echo "gnomad	HG00116
gnomad	HG00475
gnomad	HG01936
gnomad	HGDP00102
gnomad	HGDP00562
gnomad	LP6005441-DNA_E10
" > $gnomad_remove

###############################################################################
# find common snps and merge NL and GNOMAD
plink2 --vcf $nl_data --autosome --snps-only --max-alleles 2 --make-bed --out nl.raw
plink2 --bfile $gnomad_data --autosome --snps-only --max-alleles 2 --make-bed --out gnomad.raw

cat nl.raw.bim | awk '{print $1, $1":"$4, $3,$4,$5,$6}' > tmp-nl.raw.bim
mv tmp-nl.raw.bim nl.raw.bim

cat gnomad.raw.bim | awk '{print $1, $1":"$4, $3,$4,$5,$6}' > tmp-gnomad.raw.bim
mv tmp-gnomad.raw.bim gnomad.raw.bim

cat nl.raw.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > nl.snps
cat gnomad.raw.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > gnomad.snps

join gnomad.snps nl.snps > common.snps

# extract
plink2 --bfile nl.raw --extract common.snps --make-bed --out nl.common
cat nl.common.fam | awk '{print "NL",$2,$3,$4,$5,$6}' > tmp.fam
mv tmp.fam nl.common.fam
plink2 --bfile nl.common --keep $nl_keep --remove $nl_remove --make-bed --out nl.common-filtered

plink2 --bfile gnomad.raw --extract common.snps --make-bed --out gnomad.common
cat gnomad.common.fam | awk '{print "GNOMAD",$2,$3,$4,$5,$6}' > tmp.fam
mv tmp.fam gnomad.common.fam
plink2 --bfile gnomad.common --remove $gnomad_remove --make-bed --out gnomad.common-filtered

# merge
plink --bfile nl.common-filtered --bmerge gnomad.common-filtered --make-bed --out $sequence.merged

###############################################################################
## Liftover the ATP (Irish-British) data from build37 to build38 to merge with
## NL and gnomad

# variables
irebrit_data=ire-brit_data/ATP.raw
irebrit_keep=ire-brit_data/ATP.unrelated.txt

# liftover
cat $irebrit_data.bim | awk '{print $2,$5$6}' | awk '{if ($2 != "AT" && $2 != "TA" && $2 != "GC" && $2 != "CG") print $1}' | sort > atp.snps

plink --bfile $irebrit_data --extract atp.snps --make-bed --out ATP.b37

cat ATP.b37.bim | awk '{print "chr" $1, $4, $4+1, $2}' > ATP.b37.map
sed -i 's/ /\t/g' ATP.b37.map

liftOver ATP.b37.map hg19ToHg38.over.chain.gz ATP.b38.map ATP.b38.unlifted.map
cat ATP.b38.map | awk '{print $4}' > ATP.lifted.snps
plink --bfile ATP.b37 --extract ATP.lifted.snps --make-bed --out ATP.b37.lifted
cat ATP.b38.map | awk '{print $4, $2}' > ATP.lifted.snps
plink --bfile ATP.b37.lifted --update-map ATP.lifted.snps 2 1 --make-bed --out ATP.b38

###############################################################################
## Merge b38 ATP with the NFL-gnomad merged data
sequence1=nfl-gnomad-merged
sequence2=nfl-gnomad-irebrit-merged

# common and merge
cat ATP.b38.bim | awk '{print $0,$5$6}' | awk '{if ($7 != "AT" && $7 != "TA" && $7 != "GC" && $7 != "CG") print $1 ":" $4}' | sort > ATP.snps
cat $sequence1.merged.bim | awk '{print $0,$5$6}' | awk '{if ($7 != "AT" && $7 != "TA" && $7 != "GC" && $7 != "CG") print $1 ":" $4}' | sort > /$sequence1.snps

join $sequence1.snps ATP.snps > common.snps

cat ATP.b38.bim | awk '{print $1, $1 ":" $4, $3,$4,$5,$6}' > tmp; mv tmp ATP.b38.bim
cat $sequence1.merged.bim | awk '{print $1, $1 ":" $4, $3,$4,$5,$6}' > tmp; mv tmp $sequence1.merged.bim

plink2 --bfile ATP.b38 --extract common.snps --make-bed --out ATP.common
cat ATP.common.fam | awk '{print "ATP",$2,$3,$4,$5,$6}' > tmp.fam
mv tmp.fam ATP.common.fam

plink2 --bfile $sequence1.merged --extract common.snps --make-bed --out $sequence1.merged.common

# merge
plink --bfile ATP.common --bmerge $sequence1.merged.common --make-bed --out $sequence2
plink --bfile ATP.common --flip $sequence2-merge.missnp --make-bed --out ATP.common.flipped
plink --bfile ATP.common.flipped --bmerge $sequence1.merged.common --make-bed --out $sequence2

plink --bfile ATP.common.flipped --exclude $sequence2-merge.missnp --make-bed --out ATP.common.flipped.excluded
plink --bfile $sequence1.merged.common --exclude $sequence2-merge.missnp --make-bed --out $sequence1.merged.common.excluded
plink --bfile ATP.common.flipped.excluded --bmerge $sequence1.merged.common.excluded --make-bed --out $sequence2
