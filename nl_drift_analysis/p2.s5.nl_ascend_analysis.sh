!#/bin/bash

###############################################################################
sequence=nfl-gnomad-merged

# subset markers with cM information
echo "Extracting ped/map..."

plink --bfile $sequence.filtered-cleaned \
  --extract temp.cm.extract \
  --recode --out $sequence.cm-updated

# convert with converf
echo -e "genotypename:    $sequence.cm-updated.ped
snpname:         temp.cm.map
indivname:       $sequence.cm-updated.ped
outputformat:    EIGENSTRAT
genotypeoutname: $sequence.geno
snpoutname:      $sequence.snp
indivoutname:    $sequence.ind
familynames:     NO" > $sequence.plink2eigen.par

convertf -p $sequence.plink2eigen.par > $sequence.plink2eigen.log
