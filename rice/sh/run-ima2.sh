# Author: Sang Chul Choi
# Date  : Sat Jun 11 17:00:39 EDT 2011

function run-ima2 {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION

      BURNIN=$(grep REPETITION${REPETITION}-Burnin $SPECIESFILE | cut -d":" -f2)
      CHAINLENGTH=$(grep REPETITION${REPETITION}-ChainLength $SPECIESFILE | cut -d":" -f2)
      THIN=$(grep REPETITION${REPETITION}-Thin $SPECIESFILE | cut -d":" -f2)

      echo $IMA2 \
        -b$BURNIN \
        -c3 \
        -d$THIN \
        -g$DATADIR/prior-file.txt \
        -hfg -hn30 -ha0.98 -hb0.8 \
        -i$DATADIR/summarize-sequence-data.txt \
        -l$CHAINLENGTH \
        -o$IMA2DIR/$FUNCNAME-$(date +%m%d%y%H%M%S).out \
        -r5 \
        -s$(date +%s) > log.txt
      break
    fi
  done
}
