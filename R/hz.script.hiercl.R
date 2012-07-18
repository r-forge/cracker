hz.script.hiercl <-
function(sclus,gclus, p.aov,.col,plot.type,gui.input){
		width.temp <- 7
if(length(sclus$labels) > 30){
	width.temp <- width.temp +  (length(sclus$labels)-30)*0.25
	
}

mai.bottom.temp <- 1
if(max(nchar(sclus$labels)) > 10){
	mai.bottom.temp	<- mai.bottom.temp + (max(nchar(sclus$labels))-10)*0.086
	
}
temp.col <- .col
if(all(temp.col == "white")){temp.col <- rep("grey20",length(temp.col))}
		
	if(gui.input$graphic.type == "pdf"){
		pdf("hierarchical clustering samples.pdf",width = width.temp)
	}else{
		postscript("hierarchical clustering samples.eps",width = width.temp,height = 7, paper = "special",onefile = FALSE,horizontal = FALSE)	
	}
		par(mai = c(mai.bottom.temp,1,0,0.1))
	temp.lwd <- 2.5
	col.temp <- "grey35"
	
	try(dendro.sclus <- as.dendrogram(	sclus,
							hang =0
							
							))
	try(dendro.sclus <- dendrapply(dendro.sclus,function(x){hz.change.nodePar(x,sclus,temp.col,temp.lwd,col.temp)}))
	try(plot(dendro.sclus,lwd = temp.lwd, edgePar = list(lwd = temp.lwd,col = col.temp),axes = FALSE))
	try(axis(2,lwd = temp.lwd-0.5,col = col.temp))
		
	dev.off()
	
			width.temp <- 7
if(length(gclus$labels) > 30){
	width.temp <- width.temp +  (length(gclus$labels)-30)*0.25
	
}

mai.bottom.temp <- 1
if(max(nchar(gclus$labels)) > 10){
	mai.bottom.temp	<- mai.bottom.temp + (max(nchar(gclus$labels))-10)*0.086
	
}
temp.col <- .col

aov.merge.gclus <- hz.merge.control(p.aov[,1],gsub(" ","", gclus$labels))

p.v.col <- p.aov[aov.merge.gclus,]
col.aov <- c()

uncor.col <- "grey30"
cor.col <- colors()[639]
col.aov[as.numeric(p.v.col[,2]) <= gui.input$p.value] <- uncor.col
col.aov[as.numeric(p.v.col[,3]) <= gui.input$p.value] <- cor.col
col.aov[as.numeric(p.v.col[,2]) > gui.input$p.value] <- "grey80"

if(plot.type == 2){
	col.aov <- col.aov
}


if(all(temp.col == "white")){temp.col <- rep("grey20",length(temp.col))}
		
	if(gui.input$graphic.type == "pdf"){
		pdf("hierarchical clustering proteins.pdf",height = width.temp,width = 28)
	}else{
		postscript("hierarchical clustering proteins.eps",height = width.temp,width = 27, paper = "special",onefile = FALSE,horizontal = FALSE)	
	}
		par(mai = c(1,1,1,mai.bottom.temp))
	temp.lwd <- 1
	col.temp <- "grey35"
	
	plot.clustering <- gclus
	try(dendro.gclus <- as.dendrogram(	gclus
							#hang =
							
							))
							
	 						
							
	try(dendro.gclus <- dendrapply(dendro.gclus,function(x){hz.change.nodePar(x,sclus,temp.col,temp.lwd,col.temp)}))
	
	try(plot(dendro.gclus,lwd = temp.lwd, edgePar = list(lwd = temp.lwd,col = "grey80"),axes = FALSE,horiz = T,xlab = "height"))
	try(axis(1,lwd = temp.lwd-0.5,col = col.temp))
		legend("topleft",c(paste("p.value <",gui.input$p.value),paste("corrected p.value <",gui.input$p.value)),fill = c(uncor.col , cor.col),cex = 2,border = "transparent",bg = "#FFFFFF99",box.col = "transparent")
	dev.off() 
	return(list(sclus = sclus,temp.col = temp.col,temp.lwd=temp.lwd,col.temp= col.temp))
	}
