#!/bin/bash

function initialize {
echo Create symbolic links at the local base directory.
echo Choose your local and storage base directoires.
echo Delete the line of exit in the script, and execute this shell script again.
NAMEANALYSIS=Tomato
LOCALBASEDIR=$HOME/Documents/Projects/$NAMEANALYSIS-Analysis
STORAGEBASEDIR=$HOME/data/$NAMEANALYSIS-Analysis
# exit 

mkdir -p $LOCALBASEDIR
mkdir -p $STORAGEBASEDIR
for i in data downloads email log output species; do
mkdir $STORAGEBASEDIR/$i
done

BASEDIR=`pwd`
echo Creating soft links to directories in the source code base directory...
for i in R doc pl sh src run COPYING; do
rm $LOCALBASEDIR/$i 
ln -s $BASEDIR/$i $LOCALBASEDIR/$i 
done

echo Creating soft links to directories in the storage base directory...
for i in data downloads email log output species; do
rm $LOCALBASEDIR/$i 
ln -s $STORAGEBASEDIR/$i $LOCALBASEDIR/$i
done

echo Creating conf directory in the local base directory...
echo Copying a local conf file to the create conf directory...
echo Edit the conf file available at $LOCALBASEDIR/conf/README
# mkdir $LOCALBASEDIR/routput
# mkdir $LOCALBASEDIR/conf
# cp $BASEDIR/conf/README.local $LOCALBASEDIR/conf/README
}
