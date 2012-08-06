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

hz.start <- function(){
	#.path1 <- "/Users/henno/documents/Skripte/R-Functions/cRacker1.483-function/cRacker/"
.path1 <- paste(path.package("cRacker"),"data",sep = "/")
path1 <- .path1
setwd(.path1)

require("cRacker")

license.text <-  paste( "Copyright (C) <2011>  <Henrik Zauber; Max Planck Institute for Molecular Plant Physiology>
 This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
")
    
 cat(license.text)


# Set cran mirror
repos <- getOption("repos")
repos["CRAN"] <- "http://mirrors.softliste.de/cran"
options(repos = repos)

### Loading required packages:
cat("Started cRacker script suite..\n")
cat("Installing/Loading Packages:\n")
check 	<- library()
check	<- check$results[,1]

if(length(grep("pcaMethods",check)) == 0){
try(source("http://www.bioconductor.org/biocLite.R"))
try(biocLite("pcaMethods"))
}

cat("Loading pcaMethods:\n")
try(library("pcaMethods"))

.packages <- c("tcltk2","gplots"
#,"igraph",
#"plotrix",
,"gtools"
,"gdata"
,"caTools"
)
for(i in 1:length(.packages)){
	temp.i <- .packages[i]
if(length(grep(temp.i,check)) == 0){
temp.error <- try(	install.packages(temp.i,repos = "http://cran.r-project.org", dependencies = FALSE)
)

}}

for(i in c(.packages,"pcaMethods")){
	print(i)
temp.i <- grep(i,library()$results)
if(length(temp.i)== 0){
	require("tcltk")
			tkmessageBox(title="Warning",message=paste("The package",i,"is not available. Please restart cRacker after you connected to the internet. Please run cRacker with administrator rights."),icon="error",type="ok")

}

}


		#library(rgl)
#plotrary("tcltk2")
cat("Loading tcltk:\n")
require("tcltk")
require("tcltk2")
tk2font.set("TkDefaultFont",settings= "-family Tahoma -size 10 -weight normal")   

### misc
show.path.mac <- function(){
input <-getwd()
try(system(paste("open ",input), intern = TRUE, ignore.stderr = TRUE))}
show.path.win <- function(){
input <-getwd()
try(system(paste("explorer ",input), intern = TRUE, ignore.stderr = TRUE))}

# read paths:
ReadAffy 		<- function(celfile.path) return(celfile.path)
print("gerg")


.import.list 	<- read.csv(paste(.path1,"/import-config",sep = ""), stringsAsFactors = FALSE)



if(length(grep("hz.matrix.creator",ls())) == 0&1==0){
	
cat("Loading cRacker functions..")
	.path1 <- print(sub("\\\\","/",.path1))
	setwd(.path1)
	wd.install <- .path1
	if(length(grep("cracker-functions",.path1)) == 0){
	setwd("./cracker-functions")
	}
	
files <- list.files()
files <- files[#c(
				grep("hz",files,fixed = TRUE)#,
				#grep("cracker",files,fixed = TRUE))
				]
if(length(files) > 0){
for(i in 1 : length(files)){
	print(files[i])
	source(files[i])
	}
}

#library("tcltk")
cat("The following functions have been loaded:\n")
print(files)
}


.path2.set 		<- hz.path.set(import.list = .import.list,path1 = .path1)
.path2 			<-	normalizePath(.path2.set$path)
	ratio.prog <- 10000


.path2.test 		<- class(try(setwd(.path2)))
if(.path2.test == "try-error"){
 	.path2.input.file	<- basename(.path2)
	.path2				<- dirname(.path2)
}else{
	.path2.input.file <- ""
}	




value = "yes"

while(value == "yes"){

setwd(.path1)
print(.path2)
test.time  <- system.time(hz.script.return<- hz.script(path2 = .path2,path2.set = .path2.set, path2.input.file=  .path2.input.file, import.list = .import.list,path1 = .path1))	

message <- paste("The cRacker has been eaten!\nOutput has been written into:\n\n",getwd(),"\n\nDo you like to run a new session?")



require(tcltk)
value <- tkmessageBox(message =message,    icon = "question", type = "yesnocancel", default = "yes",title = "Thx for using cRacker..")
 
value <- tclvalue(value)
if(value != "yes"){



try(hz.show.path())
}

}
	
}
#hz.start(.path1 <- "/Users/henno/documents/Skripte/R-Functions/cRacker1.483-function/cRacker/")