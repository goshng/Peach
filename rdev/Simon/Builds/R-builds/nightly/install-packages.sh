#!/bin/sh
if [ -z "$1" -o x$1 == x-h -o x$1 == x--help ]; then
    echo ""
    echo " $0 <oscode>-<arch> <r-source>"
    echo ""
    echo " default: <directory> native <oscode>-<arch>"
    echo " example: $0 R-2.14.0 snowleopard-i386"
    echo ""
    echo " sources are expected in <directory>"
    exit 0;
fi
RD=$1
TARGET=$2
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/BiocInstaller_1.2.1.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/Biobase_2.14.0.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/akima_0.5-4.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/AnnotationDbi_1.16.10.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/annotate_1.32.1.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/Biobase_2.14.0.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/DBI_0.2-5.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/genefilter_1.36.0.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/geneplotter_1.32.1.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/IRanges_1.12.5.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/limma_3.10.0.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/locfit_1.5-6.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/RColorBrewer_1.0-5.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/RSQLite_0.11.1.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/xtable_1.6-0.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/DESeq_1.6.1.tar.gz
$TARGET/$RD/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads/edgeR_2.4.1.tar.gz
$TARGET/$RD/bin/Rscript R/x2.R $TARGET
echo "Check $TARGET/not-yet-installd.txt"
