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
#     \                 |       .'   cRacker             \
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




	

#### Options to run this script:
# 1. set your path:
# 		Insert into the quotation marks the path that leads to your cRacker folder. It should end with /. \ are not allowed!		
# 2. save the script
# 3. run the script:
# mac: cmd+e
# windows: drag and drop script in the console and press enter



### path 1 #

rm(list=ls())

### Change here the path to the folder which contains your newest cRacker distribution! NO backslashes allowed! Path must end with /:
#######CRACKERPATH#####################
path1 <- "./cRacker1.470b8/cRacker/"
#######################################

####
# Change nothing here:
####
setwd(path1)
source("initial-script.R")
