hz.script.bp <-
function(.data2,gui.input,.col,prog.max,ui,pb){
if(length(.data2$x.sd) != 0 & length(which(unique(is.na(.data2$x.sd) == FALSE))) > 0){
		if(!exists("ratio.prog")){ratio.prog <- 1000}

	print("plot SD distribution")
	boxplot.data <- apply(as.matrix(.data2$x.sd),2,as.numeric)
	length.char <- max(nchar(colnames(boxplot.data)))*0.1
	if(length.char < 1){
		length.char <- 1
	}
box.nchar <- max(nchar(colnames(boxplot.data)))
if(box.nchar < 8){
	box.nchar <- 1
	}else{
	box.nchar <- 1 + (box.nchar-8)*0.12	
	}
	
	if(!all(.data2$prot.n == 0)){
show.sd.data 	<- boxplot.data/(apply(.data2$prot.n,2,as.numeric)^(0.5))
rownames(show.sd.data) <- rownames(.data2$prot.n)

	}
	
	
if(gui.input$graphic.type == "pdf"){
	pdf("boxplot-sd.pdf",height = 6+length.char)
	par(mai = c(box.nchar,0.9,0.5,0.1))

}else{
	dir.create(.wd.set <- "relative-sd-eps")
}

if(gui.input$graphic.type == "eps"){
	postscript(paste(".",.wd.set,"boxplot-sd.eps",sep = "/"),height = 6+length.char,width = 7, paper = "special",onefile = FALSE,horizontal = FALSE)
	par(mai = c(box.nchar,0.9,0.5,0.1))
}

		boxplot(boxplot.data,ylab = "relative sd",main = "Boxplot of relative SDs",las = 2,col = .col)
		
		if(max(boxplot.data,na.rm = T) >3){
			if(gui.input$graphic.type == "eps"){
				postscript(paste(".",.wd.set,"boxplot-sd-zoom.eps",sep = "/"),height = 6+length.char,width = 7, paper = "special",onefile = FALSE,horizontal = FALSE)
				par(mai = c(box.nchar,0.9,0.5,0.1))
			}


			boxplot(boxplot.data,ylab = "relative sd",main = "ZOOM Boxplot of relative SDs",las = 2,ylim = c(0,2), col = .col)
		}
		
		if(exists("sd.error")){
			
			if(gui.input$graphic.type == "eps"){
				postscript(paste(".",.wd.set,"boxplot-se.eps",sep = "/"),height = 6+length.char,width = 7, paper = "special",onefile = FALSE,horizontal = FALSE)
				par(mai = c(box.nchar,0.9,0.5,0.1))
			}	
			
			boxplot(sd.error,ylab = "relative se",main = "Boxplot of relative SEs",las = 2, col = .col)
	
		if(max(sd.error, na.rm = T)> 3){
				if(gui.input$graphic.type == "eps"){
				postscript(paste(".",.wd.set,"boxplot-se-zoom.eps",sep = "/"),height = 6+length.char,width = 7, paper = "special",onefile = FALSE,horizontal = FALSE)
				par(mai = c(box.nchar,0.9,0.5,0.1))
			}
			
			
			boxplot(sd.error,ylab = "relative se",main = "ZOOM Boxplot of relative SEs",las = 2,ylim = c(0,2), col = .col)
		}
		
		
		}
	graphics.off()
}
}
