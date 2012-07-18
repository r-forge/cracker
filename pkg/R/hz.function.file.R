hz.function.file <-
function(input,shortcut){
	
	try(temp.x <- as.character(input[input$shortcut == shortcut,]))
	if(exists("temp.x")){
	if(length(temp.x)==0){temp.x <- "no help available"}
	return(temp.x)
	}else{return("no help available")}
	
}
