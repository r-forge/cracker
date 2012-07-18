hz.loading.cracker.functions <- function(path1){
	path1 				<- paste(path1,"/R",sep = "")
	print(path1)
	import.functions 	<- list.files(path1)
	print("Start loading cRacker functions")
	for( i in 1:length(import.functions)){
	cat(import.functions[i],"   ")	
	source(paste(path1,"/",import.functions[i],sep = ""))
	}
	print("Finished loading cRacker functions")
	
}
