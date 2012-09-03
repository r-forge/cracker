hz.script.heatmap2 <-
function(.data2,gui.input,p.aov, hz.exp.des.parse.data2,.col,colorblind.set,prog.max,pb,ui, plot.type,color.blind,ratio.prog){
	#save(.data2,gui.input,p.aov, hz.exp.des.parse.data2,.col,colorblind.set,prog.max,pb,ui, plot.type,color.blind,ratio.prog,file = "test.heatmap.Rdata")
		if(!exists("ratio.prog")){ratio.prog <- 1000}

	######## heatmap:
	norm.m <- gui.input$norm.method
	if(norm.m == "median"){norm.m <- "mean"}
	

	hm.input 			<- apply(.data2$x,2,as.numeric)
	rownames(hm.input)	<- rownames(.data2$x)
	hm.input1			<- hz.norm(hm.input,1,norm ="z")$x

	hm.input2			<- hz.norm(hm.input,1,norm ="mean")$x
	
	if(!gui.input$n15.log2){
		hm.input2 <- log2(hm.input2)
	}


	hm.input2[is.na(hm.input) != is.na(hm.input2)]<- 2
	if(dim(hm.input)[2] >2){
		hm.input.plot <- hz.norm(hm.input,norm = "z")$x
		report.heatmap.norm <- "z-score"
	}else{
		hm.input.plot <- hz.norm(hm.input,norm = "median")$x
		report.heatmap.norm <- "mean"
	}
	
	

	rownames(hm.input) 	<- rownames(.data2$x)
	hm.input.plot 	= hm.input.plot[rowSums(!is.na(hm.input))!=0, colSums(!is.na(hm.input))!=0]

	hm.input 		= hm.input[rowSums(!is.na(hm.input))!=0, colSums(!is.na(hm.input))!=0]


	
	
	hm.input.NA  <- hm.input.plot 

	hm.input.NA.inf <- hz.inf.replace(hm.input.NA,0.1)
	
	hm.input.NA 	<- hm.input.NA.inf$x
	hm.input.NA[is.na(hm.input.NA)] <- min(hm.input.NA,na.rm = T)
	

	sclus = hclust(dist(t(hm.input.NA)),method = "average")

	
	test.nchar.rownames <- rownames(hm.input.NA)
nchar.rownames 		<- nchar(test.nchar.rownames)
test.nchar.rownames[nchar.rownames > 15] <- substring(test.nchar.rownames,1,15)
test.nchar.rownames[nchar.rownames > 15] <- paste(test.nchar.rownames[nchar.rownames > 15] ,"...",sep = "")
rownames(hm.input.NA) <- test.nchar.rownames
	
	gclus.sd <- apply(hm.input.NA,1, function(x){sd(x,na.rm = T)})
	if(any(as.numeric(gclus.sd) == 0)){
		print("dist")
		gclus <- hclust(dist(hm.input.NA),method = "average")
	}else{
		print("cor")
		gclus <- hclust(dist(hm.input.NA),method = "average")


	}

	error.try <- class(.error<- try(hz.script.hiercl.return <- hz.script.hiercl(sclus,gclus, p.aov,.col, plot.type= plot.type,gui.input = gui.input)))
	print(hz.script.hiercl.return)

	if(error.try == "try-error"){
				print(.error)

		tkmessageBox(title="Message",message=paste("Error in hierarchical clustering!",.error),icon="warning",type="ok")
	}
	
	
	

	height.val = 2+0.02*dim(hm.input)[1]
	width.val = 5 + dim(hm.input)[2]*0.8
	if(height.val < 8){ 
		height.val = 8
	}
	if(width.val < 8){ 
		width.val = 8
	}

	mar.row = max(nchar(rownames(hm.input)))
	mar.col = max(nchar(colnames(hm.input)))

	heatmap.max <- max(nchar(colnames(hm.input)))

if(heatmap.max > 5){
	heatmap.max <- (heatmap.max-5)*0.53+5
	
	}

		if(colorblind.set){
				jet.colors		<- colorRampPalette(c(
					unlist(color.blind)[c(7,6,4)],
					
					unlist(color.blind)[c(2,3)]

					
					
					))
			}else{
				jet.colors <- colorRampPalette(c("#000039","#00007F", "#007FFF",colors()[639],colors()[638],colors()[635],"green", "yellow","red"))
				
		}
	
try(if(!gui.input$color.plot){
	jet.colors <- colorRampPalette(c("black",colors()[276],colors()[338],"white"))
})
	if(gui.input$graphic.type == "pdf"){
		pdf("heatmap.pdf", width = width.val,height = height.val)
	}else{
		postscript("heatmap.eps", width = width.val,height = height.val, paper = "special",onefile = FALSE,horizontal = FALSE)		
	}
	library("grDevices")
	
	try(dendro.gclus <- as.dendrogram(gclus))
	try(dendro.gclus <- dendrapply(dendro.gclus, function(x){hz.change.nodePar(x,gclus,hz.script.hiercl.return$col.aov,hz.script.hiercl.return$temp.lwd, hz.script.hiercl.return$col.aov)}))
	if(!exists("dendro.gclus")){dendro.gclus <- gclus}
	

	Colv.input <- as.dendrogram(sclus)

standard.heatmap <- TRUE
if(gui.input$time.grouped& exists("exp.design")){
error.try.heat <- class(.error<- try(hz.script.heatmap2.time(sclus,gclus,.data2,gui.input,hz.exp.des.parse.data2,prog.max,pb,ui)))



if(error.try.heat == "try-error"){
		tkmessageBox(title="Message",message=paste("Error in heatmap calculation!",.error),icon="warning",type="ok")
	standard.heatmap <- TRUE

}else{standard.heatmap <- TRUE}


	
	}else{standard.heatmap <- TRUE}
	
	if(standard.heatmap){	
	try(
		hm.data	<- heatmap.2(hm.input.plot,
							Rowv= dendro.gclus,
							Colv=as.dendrogram(sclus),
							col = jet.colors,
							cexRow = 0.05+1/dim(hm.input.plot)[1],
							cexCol = 1+1/dim(hm.input.plot)[2],
							#margins = c(10+5*mar.col/height.val,5+mar.row/10),
							tracecol = "darkgrey",
							trace= "none",
							margins = c(heatmap.max,5),
							lwd =1 
							#,scale = "row"
							)
		
	)
	}
	dev.off()
	#cbind(rownames(hm.input)[hm.data$rowInd],hm.data$rowInd)
	return(list(report.heatmap.norm,hz.script.hiercl.return = hz.script.hiercl.return))
}
