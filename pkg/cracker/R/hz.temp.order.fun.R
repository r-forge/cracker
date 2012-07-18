hz.temp.order.fun <-
function(gui.input,.data2,.design){	
		if(gui.input$group.norm & gui.input$exp.design != ""){
		
		
		if(as.logical(gui.input$raw)){
			temp.order.cols <- tolower(gsub(" ","",colnames(.data2$x)))
			temp.order <- c()
			for(i.ord in 1:length(temp.order.cols)){
				
				grep.i.ord	 <- as.character(temp.order.cols[i.ord])== tolower(as.character(.design[,1]))
				if(length(grep.i.ord) == 0){
					temp.order[i.ord] <- "Error, not mapped"
				}else{
					temp.order[i.ord] <- .design[grep.i.ord,3]
				}
			}
			if(is.null(temp.order)){temp.order <- rep("error in mapping",dim(.data2$x)[2])}	
		}
	
		if(!as.logical(gui.input$raw)){
			temp.order.cols <- tolower(gsub(" ","",colnames(.data2$x)))
			temp.order.cols <- gsub(" $","",temp.order.cols)
			temp.order <- c()
			for(i.ord in 1:length(temp.order.cols)){
				grep.i.ord	 <- as.character(temp.order.cols[i.ord]) == tolower(as.character(.design[,2]))
				if(length(grep.i.ord) == 0){
					temp.order[i.ord] <- "Error, not mapped"
					
				}else{
					
					temp.order[i.ord] <- .design[grep.i.ord,3]
					
				}
			
			}
			if(is.null(temp.order)){temp.order <- rep("error in mapping",dim(.data2$x)[2])}	
		}

	}else{
		temp.order <- rep(1,dim(.data2$x)[2])
	}
	
	return(temp.order)

}
