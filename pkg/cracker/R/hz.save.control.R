hz.save.control <-
function(template,path,object,name){
messagetext <- "ok"

if(!is.na(template)){
	grep.name <- grep(as.character(template),list.files(path))
	if(length(grep.name) > 0){
		print("Warning, filename already exists.")
		messagetext <- tk_messageBox("Warning",message = "Filename already exists.\nOverwrite existent file? ",type ="okcancel")
	}

	if(messagetext == "ok"){
		
	assign(name,object)	
	
	save(list = name,file =paste(path,template,sep = ""), precheck = TRUE)
	print(paste("Saved file in",paste(path,template,sep = "")))
	}
		
}else{
}
}
