hz.exp.des.parse <-
function(x,exp.des,raw.type = FALSE,gui.input){
	
	exp.des <- apply(exp.des,2,as.character)
	exp.des <- apply(exp.des,2,tolower)

	
	#x <- sclus$labels
	exp.des = exp.des
	
	.exp 		<- unique(exp.des[,2])
	.exp.col 	<- rainbow(length(.exp),alpha = 0.8,s = 1,v = 0.85)
	.names 		<- x
	.pch 		<- x
	#print(.exp)
	
	if(gui.input$raw == TRUE & !all(c(gui.input$empai.sd,gui.input$calc.empai))){
	#print("exp")
	for(a in 1:length(.exp)){
		print(.exp[a])
		temp.a <- .exp[a]
		temp.a.exp <- exp.des[exp.des[,2] == temp.a,]
		
		if(is.vector(temp.a.exp)){
		temp.a.exp <- t(as.matrix(temp.a.exp))
		}
		temp.a.exp <- temp.a.exp[,1]
		
		for(i in 1:length(temp.a.exp)){
			temp.i <- grep(temp.a.exp[i],.names)
			
			if(length(temp.i) != 0){
			.names[temp.i] <- .exp.col[a]
			.pch[temp.i] <- a




			}else{
				.names[i] <- 1
				.pch[i] = 1
				#print(temp.a.exp[i])
				}
		}
	}
	

	.exp.col <- .names
	.exp.col <- sub(" ","",.exp.col)
}else{.exp.col <- 1;.pch <- 1}

return(cbind(.exp.col,x,.pch))

}
