###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of Tomato Analysis.
# 
# Tomato Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tomato Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Tomato Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################
function package {
  d=tomato-`date +%m%d%y`
  mkdir $d
  echo "Zipping $d ... to create $d.tar.gz ..."
  cp doc/README $d
  cp output/gene2gene-pval.csv $d
  cp output/go-all.csv $d
  cp -r output/circosout $d
  cp -r output/summary $d
  cp -r output/metabolite $d 
  cp -r output/metabolite-cor $d 
  cp -r output/cor-ab $d 
  cp -r output/cor $d 
  cp -r output/go $d 
  cp -r output/ilmap.pdf $d 
  tar zcf $d.tar.gz $d
  rm -rf $d
}

function batch {
  summary 
  go
  figures
}

