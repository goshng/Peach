function conf {
  CONFFILE=conf/global
  PROJECTNAME=$(grep PROJECTNAME $CONFFILE | cut -d":" -f2)
  CAC_USERNAME=$(grep CAC_USERNAME $CONFFILE | cut -d":" -f2)
  CAC_LOGIN=$(grep CAC_LOGIN $CONFFILE | cut -d":" -f2)
  CAC_ROOT=$(grep CAC_ROOT $CONFFILE | cut -d":" -f2)
  BATCHEMAIL=$(grep BATCHEMAIL $CONFFILE | cut -d":" -f2)
  BATCHACCESS=$(grep BATCHACCESS $CONFFILE | cut -d":" -f2)
  QUEUENAME=$(grep QUEUENAME $CONFFILE | cut -d":" -f2)

  # Binary files
  IMR=$(grep IMR $CONFFILE | cut -d":" -f2)
  FSC=$(grep FSC $CONFFILE | cut -d":" -f2)

  # Other global variables
  CAC_USERHOST=$CAC_USERNAME@$CAC_LOGIN
  CAC_SSHROOTDIR=$CAC_USERHOST:$CAC_ROOT
  CAC_OUTPUTDIR=$CAC_ROOT/output
  CACBASE=$CAC_ROOT/output

  # The main base directory contains all the subdirectories.
  OUTPUTDIR=$ROOTANALYSISDIR/output
}
