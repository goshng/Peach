#!/bin/bash

function initialize {
NAMEANALYSIS=Tomato
LOCALBASEDIR=$HOME/Documents/Projects/$NAMEANALYSIS-Analysis
STORAGEBASEDIR=$HOME/data/$NAMEANALYSIS-Analysis
echo WARNING! Be careful in calling this menu because it would create
echo local and storage directories. 
echo Local base directory at $LOCALBASEDIR
echo Storage base directory at $STORAGEBASEDIR
echo We also create symbolic links at the local base directory.
echo Choose your local and storage base directoires by editting the shell script
echo sh/make-soft-links.sh.
echo Delete the line of exit in the script of sh/make-soft-links.sh
echo and execute this shell script again.
exit 

mkdir -p $LOCALBASEDIR
mkdir -p $STORAGEBASEDIR
for i in data downloads email log output species; do
mkdir $STORAGEBASEDIR/$i
done

BASEDIR=`pwd`
echo Creating soft links to directories in the source code base directory...
for i in c R doc pl sh src run tex COPYING; do
rm $LOCALBASEDIR/$i 
ln -s $BASEDIR/$i $LOCALBASEDIR/$i 
done

echo Creating soft links to directories in the storage base directory...
for i in data downloads email log output species; do
rm $LOCALBASEDIR/$i 
ln -s $STORAGEBASEDIR/$i $LOCALBASEDIR/$i
done
}
