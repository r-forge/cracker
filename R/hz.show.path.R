hz.show.path <-
function(path = NULL){
	if(!is.null(path)){
		input <- path
	}else{
		input <- getwd()
	}
	if(.Platform$OS.type == "unix"){
		try(system(paste("open ",input), intern = TRUE, ignore.stderr = TRUE))
	}else{
		#try(system(paste("explorer ",input), intern = TRUE, ignore.stderr = TRUE))
	}
}
