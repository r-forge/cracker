hz.script.heatmap2.time <-
function(sclus,gclus,.data2,gui.input,hz.exp.des.parse.data2,prog.max,pb,ui){
if(gui.input$time.grouped& exists("exp.design")){
	
	colnames(hm.input.plot) <- gsub(" ","",colnames(hm.input.plot),fixed = TRUE)
	if(gui.input$raw){
		.design.used <- cbind(.design$Alternative.name,as.character(.design$Group),.design$Order,.design$Time)
	}else{
		.design.used <- unique(cbind(.design$Experiment,as.character(.design$Group),.design$Order,.design$Time))
	}
	
	temp <- .design.used[order(.design.used[,4]),]
	temp <- temp[order(as.numeric(temp[,3])),]
	#stop()
	temp.order <- aggregate(as.numeric(as.character(temp[,3])),list(temp[,2]),min)
	for(te.i in unique(temp[,2])){
		temp[temp[,2] == te.i,3] <- as.character(temp.order[temp.order[,1] == te.i,2])
		
	}
	temp <- temp[order(as.numeric(temp[,3])),]
	
	temp.heat.order <- hz.merge.control(colnames(hm.input.plot),as.character(colnames(hm.input.plot)))	
	
	
	hm.input.plot2 <- hm.input.plot[,temp.heat.order]
Colv.input <- FALSE
}
#Colv.input <- as.dendrogram(sclus)

temp.sep <- (table(temp[temp.heat.order,2]))
order.temp.sep <- hz.merge.control(names(temp.sep),unique(temp[,c(2,3)])[,1])
temp.sep <- temp.sep[order.temp.sep]
order.col <- hz.merge.control(gsub(" ","",hz.exp.des.parse.data2[,2],fixed = T),colnames(hm.input.plot2))

.col.heat <- hz.exp.des.parse.data2[order.col,1]
.col.heat[is.na(.col.heat)]<- "white"

temp.sep.2 <- temp.sep
for(i.temp.sep in 2:length(temp.sep)){
	temp.sep.2[i.temp.sep] <- sum(temp.sep[1: i.temp.sep])
}
temp.sep.2 <- temp.sep.2[!temp.sep.2 >= dim(hm.input.plot2)[2]]
temp.sep.2 <- temp.sep.2[!is.na(temp.sep.2)] 
 
hm.data	<- heatmap.2(hm.input.plot2,
							Rowv= dendro.gclus,
							Colv=Colv.input,
							dendrogram = "row",
							col = jet.colors,
							cexRow = 0.05+1/dim(hm.input)[1],
							cexCol = 1+1/dim(hm.input)[2],
							#margins = c(10+5*mar.col/height.val,5+mar.row/10),
							tracecol = "darkgrey",
							trace= "none",
							margins = c(heatmap.max,5),
							lwd =1 ,
							ColSideColors = .col.heat,
							colsep = temp.sep.2
							#,scale = "row"
							)
temp.legend <- unique(cbind(hz.exp.des.parse.data2[order.col,1],temp[,2]))
legend(1,1,temp.legend[,2],fill = temp.legend[,1], xjust = 1,yjust = 0.5,cex = 0.5,xpd = T)
}
