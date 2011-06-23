# Author: Sang Chul Choi
# Date  : Tue Jun 14 10:28:12 EDT 2011

function sim1-imaa {
  REPETITION=1
  global-variable $SPECIES $REPETITION

  # Read options from sim1 file to run 
  OUTFILE=$(grep OUTFILE $SPECIESFILE | cut -d":" -f2)
  BURNIN=$(grep REPETITION${REPETITION}-BURNIN $SPECIESFILE | cut -d":" -f2)
  TREESAMPLESIZE=$(grep REPETITION${REPETITION}-TREESAMPLESIZE $SPECIESFILE | cut -d":" -f2)
  THIN=$(grep REPETITION${REPETITION}-THIN $SPECIESFILE | cut -d":" -f2)
  QPRIOR=$(grep REPETITION${REPETITION}-QPRIOR $SPECIESFILE | cut -d":" -f2)
  MPRIOR=$(grep REPETITION${REPETITION}-MPRIOR $SPECIESFILE | cut -d":" -f2)
  TPRIOR=$(grep REPETITION${REPETITION}-TPRIOR $SPECIESFILE | cut -d":" -f2)
  PRINTINTERVAL=$(grep REPETITION${REPETITION}-PRINTINTERVAL $SPECIESFILE | cut -d":" -f2)
  echo $IMA2 \
    -z$PRINTINTERVAL \
    -a1 \
    -b$BURNIN \
    -d$THIN \
    -q$QPRIOR \
    -m$MPRIOR \
    -t$MPRIOR \
    -i$DATADIR/$OUTFILE \
    -l$TREESAMPLESIZE \
    -o$IMA2DIR/$OUTFILE-$(date +%m%d%y%H%M%S).ima \
    -r5 \
    -s$(date +%s) > log.txt
    # -hfg -hn30 -ha0.99 -hb0.8 \
}
