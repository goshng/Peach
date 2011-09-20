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
export KENT=/Users/goshng/Documents/Projects/ucsc/kent
export KENTPRODUCT=/Users/goshng/Documents/Projects/ucsc/kent/src/product
export SQL_PASSWORD=$USER

#########################################################################
# Let's add a new genome to streptoccocus genus.
# Change the DBNAME.
# strMut2 is the genome for RNA-Seq data.
DBNAME=strMut2
REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.fna
cp $REFGENOMEFASTA .
# Edit the genome FASTA-file so that the header is chr1.
DBNAME=strMut2
REFGENOMEFASTA=NC_004350.fna
hgFakeAgp -minContigGap=1 $REFGENOMEFASTA $DBNAME.agp
faToTwoBit $REFGENOMEFASTA $DBNAME.2bit
mkdir -p /gbdb/$DBNAME/html
cp $DBNAME.2bit /gbdb/$DBNAME
# sort -k1,1 -k2n,2n $DBNAME.agp > $DBNAME.agp.2
# checkAgpAndFa $DBNAME.agp.2 $DBNAME.2bit
twoBitInfo $DBNAME.2bit stdout | sort -k2nr > chrom.sizes
rm -rf bed
mkdir -p bed/chromInfo
awk '{printf "%s\t%d\t/gbdb/strMut1/strMut1.2bit\n", $1, $2}' \
  chrom.sizes > bed/chromInfo/chromInfo.tab
hgsql -e "create database $DBNAME;" mysql
hgsql $DBNAME < $KENT/src/hg/lib/grp.sql
cp bed/chromInfo/chromInfo.tab /tmp/
hgLoadSqlTab $DBNAME chromInfo $KENT/src/hg/lib/chromInfo.sql \
  /tmp/chromInfo.tab
rm /tmp/chromInfo.tab
hgGoldGapGl $DBNAME $DBNAME.agp
#
DBNAME=strMut2
mkdir bed/gc5Base
hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 $DBNAME \
  $DBNAME.2bit | wigEncode stdin bed/gc5Base/gc5Base.{wig,wib}
hgLoadWiggle -pathPrefix=/gbdb/$DBNAME/wib \
  $DBNAME gc5Base bed/gc5Base/gc5Base.wig
mkdir -p /gbdb/$DBNAME/wib/bed/gc5Base
cp bed/gc5Base/gc5Base.wib /gbdb/$DBNAME/wib/bed/gc5Base
# Now, edit the database.
hgsql hgcentral < files/dbDbInsert.sql
# Grant permissions.
mysql -u root -p${SQL_PASSWORD} -e "GRANT FILE ON *.* to browser@localhost \
	IDENTIFIED BY 'genome';" mysql
for DB in strMut2 # hgcentral hg19 hg18 strMut1 hgFixed # proteins040315
do
    mysql -u root -p${SQL_PASSWORD} -e "GRANT SELECT, INSERT, UPDATE, DELETE, \
	CREATE, DROP, ALTER, CREATE TEMPORARY TABLES on ${DB}.* \
	TO browser@localhost \
	IDENTIFIED BY 'genome';" mysql
done
for DB in strMut2 # hgcentral hg19 hg18 strMut1 hgFixed # cb1 proteins040315
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
# Edit the local tracks

