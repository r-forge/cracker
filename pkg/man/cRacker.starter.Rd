\name{cRacker.starter}
\alias{cRacker.starter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial starter for cRacker%%  ~~function to do ... ~~
}
\description{
Use this function to start a cRacker analysis of your quantified proteomic peptide lists.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cRacker.starter()
}
% %- maybe also 'usage' for other objects documented here.
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }
% 
% %% ~Make other sections like Warning with \section{Warning }{....} ~
% 
% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% 
% ## The function is currently defined as
% function () 
% {
%     .path1 <- paste(path.package("cRacker"), "data", sep = "/")
%     path1 <- .path1
%     setwd(.path1)
%     require("cRacker")
%     license.text <- paste("Copyright (C) <2011>  <Henrik Zauber; Max Planck Institute for Molecular Plant Physiology>\n This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\n\nThis program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.\n\nYou should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>\n")
%     cat(license.text)
%     repos <- getOption("repos")
%     repos["CRAN"] <- "http://mirrors.softliste.de/cran"
%     options(repos = repos)
%     cat("Started cRacker script suite..\n")
%     cat("Installing/Loading Packages:\n")
%     check <- library()
%     check <- check$results[, 1]
%     if (length(grep("pcaMethods", check)) == 0) {
%         try(source("http://www.bioconductor.org/biocLite.R"))
%         try(biocLite("pcaMethods"))
%     }
%     cat("Loading pcaMethods:\n")
%     try(library("pcaMethods"))
%     .packages <- c("tcltk2", "gplots", "gtools", "gdata", "caTools")
%     for (i in 1:length(.packages)) {
%         temp.i <- .packages[i]
%         if (length(grep(temp.i, check)) == 0) {
%             temp.error <- try(install.packages(temp.i, repos = "http://cran.r-project.org", 
%                 dependencies = FALSE))
%         }
%     }
%     for (i in c(.packages, "pcaMethods")) {
%         print(i)
%         temp.i <- grep(i, library()$results)
%         if (length(temp.i) == 0) {
%             require("tcltk")
%             tkmessageBox(title = "Warning", message = paste("The package", 
%                 i, "is not available. Please restart cRacker after you connected to the internet. Please run cRacker with administrator rights."), 
%                 icon = "error", type = "ok")
%         }
%     }
%     cat("Loading tcltk:\n")
%     require("tcltk")
%     require("tcltk2")
%     show.path.mac <- function() {
%         input <- getwd()
%         try(system(paste("open ", input), intern = TRUE, ignore.stderr = TRUE))
%     }
%     show.path.win <- function() {
%         input <- getwd()
%         try(system(paste("explorer ", input), intern = TRUE, 
%             ignore.stderr = TRUE))
%     }
%     ReadAffy <- function(celfile.path) return(celfile.path)
%     print("gerg")
%     .import.list <- read.csv(paste(.path1, "/import-config", 
%         sep = ""), stringsAsFactors = FALSE)
%     if (length(grep("hz.matrix.creator", ls())) == 0 & 1 == 0) {
%         cat("Loading cRacker functions..")
%         .path1 <- print(sub("\\\\", "/", .path1))
%         setwd(.path1)
%         wd.install <- .path1
%         if (length(grep("cracker-functions", .path1)) == 0) {
%             setwd("./cracker-functions")
%         }
%         files <- list.files()
%         files <- files[grep("hz", files, fixed = TRUE)]
%         if (length(files) > 0) {
%             for (i in 1:length(files)) {
%                 print(files[i])
%                 source(files[i])
%             }
%         }
%         cat("The following functions have been loaded:\n")
%         print(files)
%     }
%     .path2.set <- hz.path.set(import.list = .import.list, path1 = .path1)
%     .path2 <- normalizePath(.path2.set$path)
%     ratio.prog <- 10000
%     .path2.test <- class(try(setwd(.path2)))
%     if (.path2.test == "try-error") {
%         .path2.input.file <- basename(.path2)
%         .path2 <- dirname(.path2)
%     }
%     else {
%         .path2.input.file <- ""
%     }
%     value = "yes"
%     while (value == "yes") {
%         setwd(.path1)
%         print(.path2)
%         test.time <- system.time(hz.script.return <- hz.script(path2 = .path2, 
%             path2.set = .path2.set, path2.input.file = .path2.input.file, 
%             import.list = .import.list, path1 = .path1))
%         message <- paste("The cRacker has been eaten!\nOutput has been written into:\n\n", 
%             getwd(), "\n\nDo you like to run a new session?")
%         require(tcltk)
%         value <- tkmessageBox(message = message, icon = "question", 
%             type = "yesnocancel", default = "yes", title = "Thx for using cRacker..")
%         value <- tclvalue(value)
%         if (value != "yes") {
%             try(hz.show.path())
%         }
%     }
%   }
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
