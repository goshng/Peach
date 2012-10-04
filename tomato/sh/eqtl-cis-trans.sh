# IL12-2
NUMBERCHROM=12
  rm -rf output/position4; mkdir output/position4
  rm -rf output/cis2; mkdir output/cis2
  rm -rf output/trans2; mkdir output/trans2
  rm -rf output/links2; mkdir output/links2
  rm -rf output/ilname; mkdir output/ilname
  rm -rf output/cissnpinsideil ; mkdir output/cissnpinsideil
  rm -rf output/cissnpoutsideil ; mkdir output/cissnpoutsideil
  rm -rf output/transsnpil ; mkdir output/transsnpil
for i in `ls output/position3`; do
    ILNUM=$(echo $i | sed 's/\.csv//')
    CHRNUM=$(echo $i | sed 's/IL\([0-9]*\)-.*/\1/')
    if [[ "$i" =~ IL ]]; then
      Rscript R/eqtl-cis-trans.R $NUMBERCHROM \
        data/karyotype/karyotype.sol.txt \
        output/position3/$ILNUM.csv \
        output/IL_SNP_gene/${ILNUM}_SNP \
        output/ilmap2.dat \
        output/cis2 output/trans2 output/links2 output/ilname \
        output/cissnpinsideil output/cissnpoutsideil output/transsnpil
    fi
done
