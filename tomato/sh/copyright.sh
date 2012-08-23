###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNAseq Analysis.
# 
# RNAseq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNAseq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNAseq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################
function short-notice {
  cat <<EOF
RNAseq Analysis Copyright (C) 2012 Sang Chul Choi
This program comes with ABSOLUTELY NO WARRANTY; for details select menu warranty.
This is free software, and you are welcome to redistribute it
under certain conditions; select menu copyright for details.
EOF
}

function warranty {
  cat <<EOF
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
EOF
}

function copyright {
  less COPYING
}
