# Make your private trackDb with 
#	make update [DBS=dbname]
# Make it for the Cornell browser
#	make cornell [DBS=dbname]

# The cornell browser is configured to look at trackDB_cornell and
# hgFindSpec_cornell tables before the default UCSC ones.  To avoid
# cornell-specific tables/tracks getting mixed up with mirrored tracks,
# use "make cornell" to add tracks to Cornell mirror.  You may have to 
# add your genome to DBS below.

# Browser supports multiple trackDb's so that individual developers
# can change things rapidly without stepping on other people's toes. 
# Usually when updating it is best to update your own trackDb and
# test it to make sure it works and that you have cvs updated all
# of trackDb/ before doing a make cornell. Note that you
# must specify which trackDb you are using in your .hgconf file
# or in the cgi-bin-$(USER)/hg.conf file. Something like: 
# db.trackDb=trackDb_YourUserName

# note:  new group ARCHIVED_DBS created for assemblies still needed
# to support Conservation tracks in other assemblies.  These should
# have only chromInfo table on hgwbeta and RR.  make will not rebuild
# trackDb on hgwdev unless the db is specified explicitly on the 
# command line.



KENT=/Users/goshng/Documents/Projects/ucsc/kent
LOADTRACKS_OPTS=-sqlDir=${KENT}/src/hg/lib
EXEC=${KENT}/src/hg/makeDb/trackDb/loadTracks

DBS = 	hg19 \
        strMut1 

tagTypes:
	rsync -u ${KENT}/src/hg/makeDb/trackDb/tagTypes.tab .

update: tagTypes
	${EXEC} ${LOADTRACKS_OPTS} trackDb_${USER} hgFindSpec_${USER} ${DBS}

cornell: clean tagTypes
	${EXEC} -release=alpha ${LOADTRACKS_OPTS} trackDb_cornell hgFindSpec_cornell ${DBS}

# Get rid of symbolic links (created by lower-level makefiles):
clean:
	find . -type l -exec rm {} \;

