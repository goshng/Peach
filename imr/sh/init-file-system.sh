###############################################################################
# Functions: structuring file system (or creating directories)
###############################################################################
function init-file-system {
  echo -n "Creating $ROOTANALYSISDIR/output ..." 
  mkdir $ROOTANALYSISDIR/output 
  echo -e " done"
  echo -n "Creating $CAC_ROOT/output at $CAC_USERHOST ..."
  ssh -x $CAC_USERHOST mkdir -p $CAC_ROOT
  ssh -x $CAC_USERHOST mkdir -p $CAC_ROOT/output
  echo -e " done"
}
