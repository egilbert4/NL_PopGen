######################################################################
## Cell 4 - QC the merged data

sequence=nfl-gnomad-irebrit-merged

kgp_idfile=/dbfs/mnt/shared-data/ed_results/KGP3.reference.idfile
nl_idfile=/dbfs/mnt/shared-data/ed_results/NL.reference.idfile
atp_idfile=/dbfs/mnt/shared-data/ed_results/ATP.reference.idfile

# update the ids
plink2 --bfile $sequence \
  --keep nl-gnomad-irebrit-merged.keep.idfile \
  --make-bed --out $sequence.keep

plink2 --bfile $sequence.keep \
  --update-ids nl-gnomad-irebrit-merged.update.idfile \
  --make-bed --out $sequence.updated


# QC markers and inds
plink2 --bfile $sequence.updated --geno 0.05 --maf 0.01 --hwe 1e-9 midp keep-fewhet --make-bed --out $sequence.filtered
plink2 --bfile $sequence.xtrahwe-filtered --mind 0.05 --make-bed --out $sequence.filtered-cleaned

######################################################################
## Convert this merged dataset to eigenstrat format
sequence1=nfl-gnomad-irebrit-merged
sequence2=NFL-YRI-Brit

grep -e IBD $sequence1.updated.fam | awk '{print $1,$2}' > tmp1
grep -e fs $sequence1.updated.fam | awk '{print $1,$2}' > tmp2
grep -e YRI $sequence1.updated.fam | awk '{print $1,$2}' > tmp3
  
cat tmp1 tmp2 tmp3 > ${sequence2}.idfile
rm tmp1 tmp2 tmp3

plink2 --bfile $sequence1.updated \
      --geno 0.05 --maf 0.01 --hwe 1e-9 midp keep-fewhet --mind 0.05 \
      --keep ${sequence2}.idfile \
       --make-bed --out tmp

plink --bfile tmp --recode --out $sequence2.admixtools
rm tmp.*
  
grep "YRI" $sequence2.admixtools.ped | cut -f 1-6 -d ' ' | wc -l
  
# convertf
echo -e "genotypename:    $sequence2.admixtools.ped
snpname:         $sequence2.admixtools.map
indivname:       $sequence2.admixtools.ped
outputformat:    EIGENSTRAT
genotypeoutname: $sequence2.admixtools.geno
snpoutname:      $sequence2.admixtools.snp
indivoutname:    $sequence2.admixtools.ind
familynames:     YES" > $sequence2.plink2eigen.par

convertf -p $workdir/$sequence2.plink2eigen.par > $workdir/$sequence2.plink2eigen.log
  
