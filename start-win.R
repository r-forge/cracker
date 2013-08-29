#..####...#####....####....####...##..##..######..#####..
#.##..##..##..##..##..##..##..##..##.##...##......##..##.
#.##......#####...######..##......####....####....#####..
#.##..##..##..##..##..##..##..##..##.##...##......##..##.
#..####...##..##..##..##...####...##..##..######..##..##.
#........................................................
#
#   |`.              \`-.
#   |  `-.    -----   \  `.     ---------     ___
#   |  __ `.           \   `-.  ____      ----------
#   |  |R|  :.----------+-----`-.----------------.
#    |  ssssss o  o  o  o  o  o  o  o  o       ooo`---__
#     \                |        .'   cRacker             \
#       --.____________/       /__________________.------'
#                     +      .'		-------
#           -----     /    .'___
#    ----            |   .'			-------
#         ___        /  /
#                   | .'
#                   /'

##########################################################################
#	 "Copyright (C) <2011>  <Henrik Zauber; Max Planck Institute for Molecular Plant Physiology>
#
#   This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>"
##########################################################################

#### Please run this script!!!!
# mac: cmd+e
# windows: drag and drop script in the console and press enter!


#path1 <- "/Users/henno/documents/Skripte/R-Functions/cRacker1.30/"
.paths <- commandArgs(TRUE)
print(.paths)
path1 <- paste(.paths[1],"", sep = "")


setwd(path1)
source("initial-script.R")
