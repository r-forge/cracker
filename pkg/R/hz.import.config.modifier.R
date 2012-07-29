hz.import.config.modifier <- 
function(path1= NA){
	if(is.na(path1)){
		path1 <- paste(path.package("cRacker.proteomics"),"/data",sep = "")
	}
	
	import.config <- read.csv(paste(path1,"/import-config"
,sep = ""))
	import.config <- fix(import.config)
	write.csv(import.config,paste(path1,"/import-config"
,sep = ""))	
	
}
