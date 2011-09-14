# Author: Sang Chul Choi
# Date  : Fri Jun 10 15:32:07 EDT 2011

function summarize-sequence-data {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION

      echo perl pl/$FUNCNAME.pl \
        -command polymorphicsite \
        -data $DATADIR \
        $DATADIR/$FUNCNAME.txt
      echo "Check $DATADIR/$FUNCNAME.txt"
      break
      perl pl/$FUNCNAME.pl \
        -assignment \
        -command polymorphicsite \
        -data $DATADIR \
        > $DATADIR/$FUNCNAME-for-assignment.txt
      echo "Check $DATADIR/$FUNCNAME-for-assignment.txt"

      break
    fi
  done
}
