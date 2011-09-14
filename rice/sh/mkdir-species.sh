# Create direcotires for storing analyses and their results.
# ----------------------------------------------------------
# The species directory is created in output subdirectory. The cluster's file
# system is almost the same as the local one. 
# The followings are the directories to create:
# 
# /Users/goshng/Documents/Projects/mauve/output/cornell
# /Users/goshng/Documents/Projects/mauve/output/cornell/1/data
# /Users/goshng/Documents/Projects/mauve/output/cornell/1/run-mauve
# /Users/goshng/Documents/Projects/mauve/output/cornell/1/run-clonalframe
# /Users/goshng/Documents/Projects/mauve/output/cornell/1/run-clonalorigin
# /Users/goshng/Documents/Projects/mauve/output/cornell/1/run-analysis
# 
# if 
# BASEDIR=/Users/goshng/Documents/Projects/mauve/output/cornell
# 
# I use 
function mkdir-species {
  mkdir -p $DATADIR $IMA2DIR
  ssh -x $CAC_USERHOST mkdir -p $CAC_DATADIR $CAC_IMA2DIR
}


