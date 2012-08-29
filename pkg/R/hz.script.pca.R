hz.script.pca <-
function(.data2,gui.input ,hz.exp.des.parse.data2,prog.max,pb,ui,ratio.prog){
require(gplots)
	if(!exists("ratio.prog")){ratio.prog <- 1000}

	########		
	text.out <- paste("*** Starting PCA Calculations ***",collapse = "")		
	cat(paste(rep("*",nchar(text.out)),collapse = "",sep = ""),"\n")
	cat(	text.out		,"\n")
	cat(paste(rep("*",nchar(text.out)),collapse = "",sep = ""),"\n")
	#########
		data.log <- .data2$x		
		data.log <- apply(data.log,2,as.numeric)
		
		
		rownames(data.log) <- rownames(.data2$x)
		data.log[data.log == "NaN"] <- NA

		data.log.na.r <- apply(data.log,1,function(x){!all(is.na(as.numeric(x)))})
		data.log.na.c <- apply(data.log,2,function(x){!all(is.na(as.numeric(x)))})
		data.log <- data.log[data.log.na.r, data.log.na.c]
		
	if(dim(data.log)[2] == 2){
		#data.log[is.na(data.log)] <- 0	
	}
	
	if(length(grep("pcaMethods",library()$results)) > 0){
		error <- class(try(.data3 		<- pca((data.log),method = "bpca",nPcs=2)@completeObs))
		report.pca.bpca <- "Bayesian PCA missing value estimator (package pcaMethods; Stacklies 2007)"
		if(error == "try-error"){
			.data3<- t(apply(data.log,1,function(x){x[is.na(as.numeric(x))]<- median(as.numeric(x),na.rm = TRUE);return(x)}))
			report.pca.bpca <- "error in calling bpca, NA replacement with median."
		}
		t.data3	 	<- t(.data3)
	
	
	}else{
		
		print("Warning! Could not find pcaMethods package, NAs in pcadata are replaced with row median!")
		
		.rows.temp <- rownames(data.log)
		data.log <- apply(data.log,1, function(x){
			x[is.na(x)] <- median(x,na.rm = TRUE)
 			
			return(x)
		})
		report.pca.bpca <- "Missing pcaMethods package, NA replacement with median."

	
	}
	
	
	
	
	.data.pca 	<- prcomp(t.data3)



vec.order <- c()
for(.na in rownames(.data.pca$x)){
	vec.order <- c(vec.order ,c(1:length(hz.exp.des.parse.data2[,2]))[.na== hz.exp.des.parse.data2[,2]])
	
}

.exp 		<- unique(.data2$exp.design[,2])
.exp.col 	<- hz.exp.des.parse.data2[,1][vec.order]

#if(!gui.input$color.plot){
#	.exp.col <- "grey30"
#}

.names 		<- hz.exp.des.parse.data2[,2][vec.order]
.pch 		<- as.numeric(hz.exp.des.parse.data2[,3])[vec.order]

if(all(gui.input$calc.empai, gui.input$empai.sd) | all(!gui.input$calc.empai, !as.logical(gui.input$raw)) ){
	if(!exists(".design")| all(unique(.pch) == 1)){
		#.exp.col <- "grey20"	
		
	}
}



if(gui.input$graphic.type == "pdf"){
	pdf("PCA.pdf")
}
.wd <- getwd()

		if(gui.input$graphic.type == "eps"){
	postscript("PCA.eps",width = 7,height = 7, paper = "special",onefile = FALSE,horizontal = FALSE,pointsize = 8)
		par(mfrow = c(2,2))
		}

		temp.sum <- apply(.data.pca$x,2,function(x){sum(abs(x))})
		temp.sum <- temp.sum/sum(temp.sum)*100

		# correct ranges
		y.width	<- diff(.ylim	<- range(.data.pca$x[,2]))
		x.width <- diff(.xlim 	<- range(.data.pca$x[,1]))
		
		if(y.width > x.width){
			rest <- (y.width - x.width )/2
			.xlim <- .xlim + c(rest*-1,rest)
		}else{
			rest <- (x.width - y.width )/2
			.ylim <- .ylim + c(rest*-1,rest)
		}
		.xlim.1 <- .xlim
		plot(.data.pca$x[,1],.data.pca$x[,2],
		xlim = (.xlim-c(0.5*diff(.xlim),0)),
		ylim = .ylim,
		
		xlab = paste("PC1;",round(temp.sum)[1],"%"),ylab = paste("PC2;",round(temp.sum)[2],"%"),col = .exp.col,pch = as.numeric(.pch),lwd = 2.5, cex = 1.5,type = "n")


		names <- names(.data.pca$x[,1])

		text(.data.pca$x[,1],.data.pca$x[,2],names,pos = 2,col = "grey35")
		points(.data.pca$x[,1],.data.pca$x[,2],	col = .exp.col,pch = as.numeric(.pch),lwd = 2.5, cex = 1.5)
		
			
		if(dim(.data.pca$x)[2] > 2){	
		# correct ranges
		x.width <- diff(.xlim 	<- range(.data.pca$x[,1]))
		y.width	<- diff(.ylim	<- range(.data.pca$x[,3]))
	
		if(y.width > x.width){
			rest <- (y.width - x.width )/2
			.xlim <- .xlim + c(rest*-1,rest)
		}else{
			rest <- (x.width - y.width )/2
			.ylim <- .ylim + c(rest*-1,rest)
		}	
		
		plot(.data.pca$x[,1],.data.pca$x[,3],
		xlim = (.xlim-c(0.5*diff(.xlim),0)),
		ylim = .ylim,
		xlab = paste("PC1;",round(temp.sum)[1],"%"),ylab = paste("PC3;",round(temp.sum)[3],"%"),col = .exp.col,pch = as.numeric(.pch),lwd = 2.5, cex = 1.5,type = "n")

		names <- names(.data.pca$x[,1])

		text(.data.pca$x[,1],.data.pca$x[,3],names,pos = 2,col = "grey35")
		points(.data.pca$x[,1],.data.pca$x[,3],	col = .exp.col,pch = as.numeric(.pch),lwd = 2.5, cex = 1.5)
		
		# correct ranges
		x.width <- diff(.xlim 	<- range(.data.pca$x[,2]))
		y.width	<- diff(.ylim	<- range(.data.pca$x[,3]))
	
		if(y.width > x.width){
			rest <- (y.width - x.width )/2
			.xlim <- .xlim + c(rest*-1,rest)
		}else{
			rest <- (x.width - y.width )/2
			.ylim <- .ylim + c(rest*-1,rest)
		}	
		
		.xlim.3 <- .xlim

		if(dim(.data.pca$x)[2] > 2){
			plot(.data.pca$x[,2],.data.pca$x[,3],
			xlim = (.xlim-c(0.5*diff(.xlim),0)),
			ylim = .ylim,
			xlab = paste("PC2;",round(temp.sum)[2],"%"),ylab = paste("PC3;",round(temp.sum)[3],"%"),col = .exp.col,pch = as.numeric(.pch),lwd = 2.5, cex = 1.5,type = "n")

			names <- names(.data.pca$x[,1])

			text(.data.pca$x[,2],.data.pca$x[,3],names,pos = 2,col = "grey35")		
			points(.data.pca$x[,2],.data.pca$x[,3],	col = .exp.col,pch = as.numeric(.pch),lwd = 2.5, cex = 1.5)

		}
		}
		
		.nchar.factor <- max(nchar(rownames(.data.pca$x)))/2*0.02
		
		if(gui.input$graphic.type == "eps"){
	barplot2(temp.sum,las = 2,ylab = "variance in %")

	graphics.off()

	postscript("PCA-biplot.eps",width = 7,height = 7, paper = "special",onefile = FALSE,horizontal = FALSE,pointsize = 8)
		par(mfrow = c(2,2))
		}
		
plot.biplot <- 		function(){
		
		biplot(.data.pca$x[,1:2],.data.pca$rotation[,1:2]
		,xlim = (.xlim.1+c(-.nchar.factor*diff(.xlim.1),.nchar.factor*diff(.xlim.1))))
		if(dim(.data.pca$x)[2] > 2){
			.xlim.2 <- range(.data.pca$x[,2])
			biplot(.data.pca$x[,2:3],.data.pca$rotation[,2:3]
					,xlim = (.xlim.2+c(-.nchar.factor*diff(.xlim.2),.nchar.factor*diff(.xlim.2))),)

			biplot(.data.pca$x[,c(1,3)],.data.pca$rotation[,c(1,3)]
					,xlim = (.xlim.1+c(-.nchar.factor*diff(.xlim.1),.nchar.factor*diff(.xlim.1))),)

		}

barplot2(temp.sum,las = 2,ylab = "variance in %")
		}
		try(plot.biplot())

dev.off()


	if(length(rownames(.data.pca$rotation) ) == length(rownames(.data2$x))){
		rownames(.data.pca$rotation) <- rownames(.data2$x)
	}else{
		cat("ALERT! Length of rotations is not equal to rownames of matrix.")
	}
	write.csv(.data.pca$x,"PCA-components.csv")
	write.csv(.data.pca$rotation,"PCA-loadings.csv")
	
return(list(report.pca.bpca=report.pca.bpca))
}
