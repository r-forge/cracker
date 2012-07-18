\name{hz.script}
\alias{hz.script}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main cRacker script function
%%  ~~function to do ... ~~
}
\description{Main Script for performing cRacker analysis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hz.script(path1, path2.set = 
	list(path = /data/list.txt, engine = maxquant, settings = default),  import.list = default)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path1}{
Path to cRacker library data folder. By default the path.package(cRacker) is called. 
%%     ~~Describe \code{path1} here~~
}

  \item{path2.set}{
List containing path to your data containing folder or file, used quantitation engine and the path to a saved settings .Rdata file.
}


  \item{import.list}{
Optional path to saved parameter.Rdata file from of previous analysis.
%%     ~~Describe \code{import.list} here~~
}
}
