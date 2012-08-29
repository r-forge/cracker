hz.script.anova <-
function(.data2,gui.input, plot.type,prog.max,pb,ui){
if(dim(.data2$x)[2] > 1){
	if(plot.type == 1){
		.data2$aov.export.1[.data2$aov.export.1 == "Inf"] <- NA
		.data2$aov.export.1[.data2$aov.export.1 == "-Inf"] <- NA
		
		.data2$aov.export.2[.data2$aov.export.2 == "Inf"] <- NA
		.data2$aov.export.2[.data2$aov.export.2 == "-Inf"] <- NA
		

	#if(gui.input$row.norm == FALSE){
#			aov.export <-  .data2$aov.export.2
#	}else{	
	
#	}

	
	aov.export <-  .data2$aov.export.1
	
	if(gui.input$log2.test & !gui.input$n15.log2){
	aov.export[,3] <- log2(as.numeric(aov.export[,3]))	
	}
	if(gui.input$log2.test|gui.input$n15.log2 ){
	temp.samples.lab <- "tested peptides log2(intensity)"
	}else{
	temp.samples.lab <- "tested intensity"
	
	}
	pdf("density-tested peptides.pdf")
	try(plot(density(as.numeric(aov.export[,3]) ),main = temp.samples.lab))
	dev.off()
	
	pdf("qqplot.pdf")
	try(qqnorm(as.numeric(aov.export[,3])))
	curve(x*1,add=T)
	dev.off()
	p.pt	<- try(hz.aov(as.data.frame(aov.export, stringsAsFactors = FALSE) ,.data2,gui.input,TRUE,as.numeric(prog.max),pb,ui))
	

	p.aov <- p.pt$aov
	p.ttest	<- p.pt$ttest
	ttestlist <- p.pt$ttestlist
	#p.pt <- p.pt$pt
	ttestlist <- cbind(ttestlist,p.adjust(as.numeric(ttestlist[,3])))
	colnames(ttestlist)[5] <-  paste("p.value",gui.input$p.adjust.method,"corrected",sep = ".")
	
	write.csv(p.aov,"anova-p-values.csv")
	
	####
	plot.names 	<- rownames(.data2$x)


	aov.order <- hz.merge.control(gsub(" ","",p.aov[,1]),gsub(" ","",plot.names))
	.aov.new 	<- p.aov[aov.order[!is.na(aov.order)],]

	.aov.new[is.na(.aov.new[,2]),2] <- 1
	.aov.new[is.na(.aov.new[,3]),3] <- 1
	
cRacker.anova.ttest <- list(p.pt = p.pt,p.aov= p.aov, ttestlist= ttestlist,.aov.new= .aov.new,p.ttest= p.ttest)
class(cRacker.anova.ttest) <- "cRacker.anova"
return(cRacker.anova.ttest)	
	
	}else{
	plot.names 	<- rownames(.data2$x)
	
	#.aov.new <- p.all	
	}
}


}
