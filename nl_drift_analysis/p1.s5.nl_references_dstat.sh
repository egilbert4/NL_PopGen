###############################################################################
## General d stat analysis D(YRI, NL; Ireland, England)
sequence=nl-gnomad-irebrit-merged
sequence1=NL-YRI-Brit

# run the D-stats analysis

echo -e "indivname:    $sequence1.admixtools.dstat-test1.ind  
  snpname:      $sequence1.admixtools.snp
  genotypename: $sequence1.admixtools.geno
  popfilename:  $sequence1.dstat.poplist
  printsd:  YES" > $sequence1.dstat.par

echo -e "YRI fs_5 England Ireland
YRI fs_6 England Ireland
YRI fs_7 England Ireland
YRI fs_8 England Ireland
YRI fs_9 England Ireland
YRI fs_10 England Ireland
YRI fs_11 England Ireland
YRI fs_12 England Ireland
YRI fs_13 England Ireland
YRI fs_14 England Ireland
YRI fs_15 England Ireland
YRI fs_16 England Ireland
YRI fs_17 England Ireland
YRI fs_18 England Ireland
YRI fs_19 England Ireland
YRI fs_20 England Ireland
YRI fs_21 England Ireland
YRI fs_22 England Ireland" > $sequence1.dstat.poplist

# run the dstat analysis
qpDstat -p $sequence1.dstat.par > $sequence1.dstat.log
cat $sequence1.dstat.log | grep "result:" | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > $sequence1.dstat.results

###############################################################################
## dstat analysis for the Irish versus English versus NL drift
sequence=NFL-YRI-Brit

# run the D-stats analysis
for clst in fs_5 fs_6 fs_7 fs_8 fs_9 fs_10 fs_11 fs_12 fs_13 fs_14 fs_15 fs_16 fs_17 fs_18 fs_19 fs_20 fs_21 fs_22
  do
  # set up D-stats
  echo -e "indivname:    $sequence.${clst}.admixtools.dstat-test2.ind
    snpname:      $sequence.admixtools.snp
    genotypename: $sequence.admixtools.geno
    popfilename:  $sequence.${clst}.admixtools.dstat-test2.popfile
    printsd:  YES" > tmp.dstat.par
    
  qpDstat -p tmp.dstat.par > $sequence.${clst}.dstat-nldrift.log
  cat $sequence.${clst}.dstat-nldrift.log | grep "result:" | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > $sequence.${clst}.dstat-nldrift.results

  rm $workdir/tmp.*
  done
