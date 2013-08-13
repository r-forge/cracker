hz.inf.replace <-
function(x,rel.replace){
	if(all(is.infinite(as.numeric(x)))){
		.max <- 0
		.min <- 0
	}else{
		.max <- max(x[!is.infinite(x)],na.rm = T)
		.min <- min(x[!is.infinite(x)],na.rm = T)

	}
	inf.vec <- as.character(x)
	inf.vec[inf.vec != "Inf" & inf.vec != "-Inf"] <- ""
	
	
	x[x == "Inf"] <- .max+.max*rel.replace
	x[x == "-Inf"] <- .min-.min*rel.replace
	
	
	return(list(x=x,inf.vec = inf.vec))

	
}
