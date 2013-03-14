#!/bin/bash

#########################################################################
# http://stackoverflow.com/questions/230266/vi-vim-how-to-pipe-visually-selected-text-to-a-unix-command-and-append-output-to
# http://stackoverflow.com/questions/3752785/what-can-you-execute-on-visually-selected-text-in-vim
# 
# You select parts of vi buffer and execute them.
# 1. Select parts of any text using VI command `V'.
# 2. Push `:!bash' to execute them.
# 3. Push `u' to make the executed parts be reappear.

#########################################################################
# Export these shell variables for easy access. You need to run them in the BASH
# shell directly.
# Create a directory, ucsc, at $HOME/Documents/Projects
# git clone git://genome-source.cse.ucsc.edu/kent.git
# Follow the instructions available at kent/src/product
# Start a Apache Web server in Mac OSX (Mountain Lion)
# $ sudo apachectl start
# $ sudo defaults write /System/Library/LaunchDaemons/org.apache.httpd Disabled -bool false
# See
# http://reviews.cnet.com/8301-13727_7-57481978-263/how-to-enable-web-sharing-in-os-x-mountain-lion/
# for detail of setting up a Web Server in Mac OSX Mountain Lion
# Install MySQL
# $ sudo /usr/local/mysql/support-files/mysql.server start
# $ /usr/local/mysql/bin/mysqladmin -u root password 'yourpasswordhere'
export UCSCDIR=/Users/goshng/Documents/Projects/ucsc
export KENT=/Users/goshng/Documents/Projects/ucsc/kent
export KENTPRODUCT=/Users/goshng/Documents/Projects/ucsc/kent/src/product
export SQL_PASSWORD=$USER

#########################################################################
# Let's add a new genome to streptoccocus genus.
#
# SdeqATCC12394
# REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/cornell_sde1/CP002215.gbk
# SdeqGGS124
# REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_dysgalactiae_equisimilis_GGS_124_uid59103/NC_012891.fna
# SddyATCC27957
# REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/cornell_sdd/sdd.gbk
# SpyMGAS315
# REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_pyogenes_MGAS315_uid57911/NC_004070.fna
# SpyMGAS10750
# REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_pyogenes_MGAS10750_uid58575/NC_008024.fna
# Seeq4047
# REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_equi_4047_uid59259/NC_012471.fna
# 
# 1. Copy the source file: genbank or FASTA format.

KENT=/Users/goshng/Documents/Projects/ucsc/kent
KENTPRODUCT=/Users/goshng/Documents/Projects/ucsc/kent/src/product
SQL_PASSWORD=$USER

###########################################################
# BEGIN of load bacterial genome function 
###########################################################
function loadBacterialGenome {
echo -n "What is the database name for the genome? (e.g., SddyATCC27957) "
read DBNAME
REFGENOME=$DBNAME

echo -n "Will you use a genbank format for creating $DBNAME? (Type in y or n) "
read WISH
if [ "$WISH" == "y" ]; then
  echo -n "What is the file name of a genbank format? (e.g., /path/to/$DBNAME.gbk) "
  read REFGENOMEFASTA
  echo "  copying $REFGENOMEFASTA ..."
  cp $REFGENOMEFASTA $DBNAME.gbk
  echo "  converting the gbk to the FASTA file using Kent's gbToFaRa command ..."
  gbToFaRa /dev/null $REFGENOME.fna $REFGENOME.ra $REFGENOME.ta $REFGENOME.gbk 
  toUpper $REFGENOME.fna $REFGENOME.fna.upper
  faSize $REFGENOME.fna.upper
  rm $REFGENOME.fna $REFGENOME.ra $REFGENOME.ta 
  mv $REFGENOME.fna.upper $REFGENOME.fna
  echo -n "Please, Edit `pwd`/$REFGENOME.fna so that the header is chr1, and enter:"
  read DONE
fi

echo -n "Will you use a FASTA format file for creating $DBNAME? (Type in y or n) "
read WISH
if [ "$WISH" == "y" ]; then
  echo -n "What is the file name in FASTA format? (e.g., /path/to/$DBNAME.fna) "
  read REFGENOMEFASTA
  cp $REFGENOMEFASTA $DBNAME.fna
  echo -n "Please, Edit `pwd`/$DBNAME.fna so that the header is chr1, and enter:"
  read DONE
fi

echo "  creating a 2bit file of the FASTA file ..."
hgFakeAgp -minContigGap=1 $REFGENOME.fna $DBNAME.agp
faToTwoBit $REFGENOME.fna $DBNAME.2bit
mkdir -p /gbdb/$DBNAME/html
cp $DBNAME.2bit /gbdb/$DBNAME

# 5. Check agp and 2bit.
# sort -k1,1 -k2n,2n $DBNAME.agp > $DBNAME.agp.2
# checkAgpAndFa $DBNAME.agp.2 $DBNAME.2bit
echo "  creating a database ..."
twoBitInfo $DBNAME.2bit stdout | sort -k2nr > chrom.sizes
rm -rf bed
mkdir -p bed/chromInfo
awk '{printf "%s\t%d\t/gbdb/DBNAME/DBNAME.2bit\n", $1, $2}' \
  chrom.sizes > bed/chromInfo/chromInfo.tab.tmp
sed s/DBNAME/$DBNAME/g < bed/chromInfo/chromInfo.tab.tmp > bed/chromInfo/chromInfo.tab
hgsql -e "create database $DBNAME;" mysql

echo "  creating grp, chromInfo tables ..."
hgsql $DBNAME < $KENT/src/hg/lib/grp.sql
cp bed/chromInfo/chromInfo.tab /tmp/
hgLoadSqlTab $DBNAME chromInfo $KENT/src/hg/lib/chromInfo.sql \
  /tmp/chromInfo.tab
rm /tmp/chromInfo.tab
hgGoldGapGl $DBNAME $DBNAME.agp

echo "  creating GC5 track ..."
mkdir bed/gc5Base
hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 $DBNAME \
  $DBNAME.2bit | wigEncode stdin bed/gc5Base/gc5Base.{wig,wib}
hgLoadWiggle -pathPrefix=/gbdb/$DBNAME/wib \
  $DBNAME gc5Base bed/gc5Base/gc5Base.wig
mkdir -p /gbdb/$DBNAME/wib/bed/gc5Base
cp bed/gc5Base/gc5Base.wib /gbdb/$DBNAME/wib/bed/gc5Base

echo -n "Please, edit files/dbDbInsert.sql for DB name, genome name, date, and scientific name and enter:"
read DONE
hgsql hgcentral < files/dbDbInsert.sql

echo "  granting permission on the created database ..."
mysql -u root -p${SQL_PASSWORD} -e "GRANT FILE ON *.* to browser@localhost \
	IDENTIFIED BY 'genome';" mysql
# 1. SdeqATCC12394 
# 2. SdeqGGS124 
# 3. SddyATCC27957 
# 4. SpyMGAS315    
# 5. SpyMGAS10750  
DBNAME=SpyMGAS10750
SQL_PASSWORD=$USER
mysql -u root -pgoshng --database=$DBNAME < dbdump/spy2knonwgenes.txt
for DB in $DBNAME # hgcentral hg19 hg18 strMut1 hgFixed # proteins040315
do
    mysql -u root -p${SQL_PASSWORD} -e "GRANT SELECT, INSERT, UPDATE, DELETE, \
	CREATE, DROP, ALTER, CREATE TEMPORARY TABLES on ${DB}.* \
	TO browser@localhost \
	IDENTIFIED BY 'genome';" mysql
done
for DB in $DBNAME # hgcentral hg19 hg18 strMut1 hgFixed # cb1 proteins040315
do
    mysql -u root -p${SQL_PASSWORD} -e "GRANT SELECT, CREATE TEMPORARY TABLES \
	on ${DB}.* TO \
	readonly@localhost IDENTIFIED BY 'access';" mysql
done
for DB in hgcentral
do
    mysql -u root -p${SQL_PASSWORD} -e "GRANT SELECT, INSERT, UPDATE, DELETE, \
	CREATE, DROP, ALTER on ${DB}.* TO readwrite@localhost \
	IDENTIFIED BY 'update';" mysql
done
echo "Please, edit the local tracks and make it"
rm -f $DBNAME.fna $DBNAME.2bit $DBNAME.agp $DBNAME.gbk chrom.size
rm -rf bed
}
###########################################################
# END of load bacterial genome function 
###########################################################

###########################################################
# Load wiggle files of recombination intensity.
# TODO: the other 4 species should be added. The generation of the map for the
# wiggle file takes time.
function loadRI {
  WIGDIR=/Users/goshng/Documents/Projects/Mauve/output/cornellf/3/run-analysis
  WIGFILEBASE=ri1-refgenome4-map
  DBNAME=SpyMGAS315
  wigEncode $WIGDIR/$WIGFILEBASE.wig $WIGDIR/$WIGFILEBASE.temp.wig $WIGDIR/$WIGFILEBASE.wib
  hgLoadWiggle $DBNAME ri $WIGDIR/$WIGFILEBASE.temp.wig
  rm $WIGDIR/$WIGFILEBASE.temp.wig
  mv $WIGDIR/$WIGFILEBASE.wib /gbdb/$DBNAME/wib/
  hgsql $DBNAME -e "update ri set file='/gbdb/$DBNAME/wib/$WIGFILEBASE.wib'"
}

###########################################################
# Load sequence alignments in wigMaf format.
# Use /Users/goshng/Documents/Projects/Mauve/test/wigmaf/run
# And make makeDb/trackDb file.
# All of the mauve alignments were loaded.

###########################################################
# Load recombination probability in BED format.
# Let me make something in /Users/goshng/Documents/Projects/Mauve/test/wigmaf/run

###########################################################
# Load virulence genes in BED format.
# Check /Users/goshng/Documents/Projects/Mauve/test/virulence/run

###########################################################
# Load recombination rate per block in BED format.
#

###########################################################
# Load recombination probability.
# 

