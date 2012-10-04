#!/bin/bash
###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of tux.
# 
# tux is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tux is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Mauve Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

function initialize {
  OUTSRC=output/src
  mkdir $OUTSRC/b
  cp src/ezlogger/* $OUTSRC
  cp src/simpleopt/* $OUTSRC
  cp downloads/FindGSL.cmake $OUTSRC
  cp src/gpl-3.0.txt $OUTSRC
}
