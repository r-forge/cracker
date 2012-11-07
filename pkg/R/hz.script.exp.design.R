hz.script.exp.design <-
function(exp.design ,gui.input, colorblind.set, color.blind,.data2){
hz.exp.des.parse.data2 <- hz.exp.des.parse(x = colnames(.data2$x) ,exp.des = exp.design,raw.type = .data2$gui.input$raw,gui.input = gui.input)	
if(dim(exp.design)[1] > 2){.design <- exp.design}
.col 			<- as.numeric(hz.exp.des.parse.data2[,3])
	if(length(unique(.col))> 1){
		if(colorblind.set){
			.col.rainbow		<- colorRampPalette(unlist(color.blind)[-1])(length(unique(.col)))
		}else{
			.col.rainbow 		<- rainbow(length(unique(.col)),alpha = 0.8,s = 1,v = 0.85)
		}
		.col.rainbow    <- sample(.col.rainbow)
		.col <- .col.rainbow[.col]
	}else{.col <- "white"}
	
hz.exp.des.parse.data2[,1] <- .col

if(all(gui.input$calc.empai, gui.input$empai.sd) | all(!gui.input$calc.empai, !as.logical(.data2$gui.input$raw)) |all(
#gui.input$plot.only != "" & 
gui.input$exp.design != "")){
	
		if(colorblind.set){
			for(c.i in 1:length(unique(hz.exp.des.parse.data2[,1]))){
				hz.exp.des.parse.data2[hz.exp.des.parse.data2[,1] == unique(hz.exp.des.parse.data2[,1])[c.i],1]		<- colorRampPalette(unlist(color.blind)[-1])(length(unique(hz.exp.des.parse.data2[,1])))[c.i]
			
			}
		}else{
			hz.exp.des.parse.data2[,1] <- sample(rainbow(length((hz.exp.des.parse.data2[,1])),alpha = 0.8,s = 1,v = 0.85))
		}
	
	
	
#"white"
	hz.exp.des.parse.data2 <- hz.exp.des.parse.data2[order(hz.exp.des.parse.data2[,2]),]

	hz.exp.des.parse.data2[,3] <- 1:dim(hz.exp.des.parse.data2)[1]
	
	if(gui.input$exp.design != "" &  gui.input$plot.only != "" ){
	
	.design  <- read.table(gui.input$exp.design,header = TRUE)
	.design[,2] <- tolower(make.names(.design[,2],allow = F))
	.design[,1] <- tolower(make.names(.design[,1],allow = F))
	.design[,5] <- tolower(make.names(.design[,5],allow = F))

	
	}	
	
	if(exists(".design")){


			.design[,2] <- tolower(make.names(.design[,2],allow = F))
			.design[,1] <- tolower(make.names(.design[,1],allow = F))
			.design[,5] <- tolower(make.names(.design[,5],allow = F))

			.design[,2] <- tolower(make.names(.design[,2],allow = F))
			if(length(grep(colnames(.data2$x),as.character(.design[,1])))!=0){
				.design[,1] <- tolower(make.names(.design[,1],allow = F))
			}else{
				.design[,1] <- tolower(make.names(.design[,5],allow = F))

			}
			
			if(any(!.design$Order == "")){
				if(!as.logical( gui.input$raw)){
					order.vec <- c()
					for(f in 1:length(hz.exp.des.parse.data2[,2])){
						print(unique(.design$Experiment)[f])
						
						temp.f <- min(.design[.design$Experiment == gsub(" ","",hz.exp.des.parse.data2[,2])[f],]$Order,na.rm = T)
						order.vec <- c(order.vec,temp.f)
					}
					order.vec[is.infinite(order.vec)] <- max(order.vec[!is.infinite(order.vec)])+1


				}else{
					order.vec <- c()

					for(f in 1:length(hz.exp.des.parse.data2[,2])){
						print(unique(.design$Name)[f])
						
						temp.f <- .design[.design$Alternative.name == gsub(" ","",hz.exp.des.parse.data2[,2])[f],]$Order
						order.vec <- c(order.vec,temp.f)
					}

				}
				
				
				
			}else{
				order.vec <- ""
			}
			
			if(as.logical(.data2$gui.input$raw ) & !all(c(gui.input$empai.sd,gui.input$calc.empai))){
					temp.design.raw <- 5
				}else{
					temp.design.raw <- 2
				}
			
			.design.uni <- unique(.design[,c(temp.design.raw,3)])
			.design.pch <- c()
			
			for(temp.it in 1:length(hz.exp.des.parse.data2[,2])){
				temp.i <- tolower(hz.exp.des.parse.data2[temp.it,2])
				temp.i <- gsub(" $","", temp.i)
				
				temp.i <- unique(.design.uni[temp.i == .design.uni[,1],])
				
				if(length(temp.i[,1]) >1){
					#temp.i[,2] <- sum(temp.i[,2])
					temp.i		<- temp.i[1,]
					print(paste("error in exp.design:",temp.i))
				}
				.design.pch <-c(.design.pch, temp.i[,2])
		
			}
			for(temp.t in 1:length(unique(.design.pch))){
				.design.pch[.design.pch==unique(.design.pch)[temp.t]] <- temp.t
				
			}			
			if(length(.design.pch != 0)){
				if(colorblind.set){
					.design.col <- colorRampPalette(unlist(color.blind)[-1])(max(.design.pch))[as.numeric(.design.pch)]
				}else{
					.design.col <- sample(rainbow(max(.design.pch)))[as.numeric(.design.pch)]
				}
			}else{
				.design.pch <- 1
				.design.col <- rgb(0,133,178,max = 255)
			}
			

			
			try(hz.exp.des.parse.data2[,1] <- .design.col)	
			try(hz.exp.des.parse.data2[,3] <- .design.pch)
			try(hz.exp.des.parse.data2 <-cbind(hz.exp.des.parse.data2,order.vec))
			
	}
	
}

.col 		<- hz.exp.des.parse.data2[,1]

return(list(hz.exp.des.parse.data2= hz.exp.des.parse.data2,.col = .col))
}
