hz.script.row.plot <-
function(.data2,gui.input,y.lab.input,.aov.new, hz.exp.des.parse.data2,colorblind.set,.col,prog.max,ratio.prog,pb,ui){
if(is.vector(.data2$x.sd)){
	
	if(!exists("ratio.prog")){ratio.prog <- 1000}

	
	.data2$x.sd <- t(as.matrix(.data2$x.sd))
}

if(is.vector(.data2$prot.n)){
	
	.data2$prot.n <- t(as.matrix(.data2$prot.n))


}

ratio.prog <- prog.max/8

row.plot.data 	<- .data2$x[order(as.character(rownames(.data2$x))),]
show.sd.data	<- .data2$x.sd[order(as.character(rownames(.data2$x.sd))),]
prot.n			<- .data2$prot.n[order(as.character(rownames(.data2$prot.n))),]


rows.n			<- rownames(prot.n)
rows.sd			<- rownames(show.sd.data)

show.sd.data 	<- apply(show.sd.data,2,as.numeric)
prot.n 			<- apply(prot.n,2,as.numeric)

rownames(prot.n) <- rows.n
rownames(show.sd.data) <- rows.sd




#if(!is.vector(sub.info)){
	print("ordering info")
try(	order.info 	<- hz.merge.control(gsub(" ","",.data2$proteinlist.info[,2]),gsub(" ","",rownames(row.plot.data))))

try(	sub.info 	<- .data2$proteinlist.info[order.info,1:3])


if(!exists("sub.info")){sub.info <- c()}




## correct n if samples are averaged
if(all(gui.input$calc.empai, gui.input$empai.sd) | all(!gui.input$calc.empai, !as.logical(gui.input$raw)) ){
	if(dim(prot.n)[2]!= dim(show.sd.data)[2]){
		rows <- rownames(prot.n)
		prot.n.comb <- c()
		colnames(prot.n)  	<- tolower(gsub("SD ","",colnames(prot.n)))
		colnames(prot.n)  	<- tolower(gsub(" ","",colnames(prot.n)))
		for(i in 1:dim(show.sd.data)[2]){
			temp.templ 	 	<- colnames(show.sd.data)[i]
			temp.templ  	<- tolower(gsub("SD ","",temp.templ))
			temp.templ  	<- tolower(gsub(" ","",temp.templ))
			temp.cols		<- .data2$exp.design[.data2$exp.design[,2] == temp.templ,1]
			temp.m <- c()
			for(find.vec in temp.cols){
				temp.m <- c(temp.m,c(1:(dim(prot.n)[2]))[find.vec == colnames(prot.n)])
			}
	
			temp.sum <- apply(prot.n[,temp.m],1,function(x){sum(x,na.rm = TRUE)})
			prot.n.comb  <- cbind(prot.n.comb,temp.sum)
		
		}
	
		prot.n <- prot.n.comb
		rownames(prot.n) <- rows
		if(gui.input$calc.empai){
			test <- merge(.data2$x,prot.n,by = 0)
		
		}	

	}
}

##
use.se <- T
if(length(show.sd.data) != 0){
	if(!all(prot.n == 0)& use.se ){
protein.plot.se <- TRUE
show.sd.data 	<- show.sd.data/(prot.n^(0.5))
sd.error		<- show.sd.data

}else{"warning! Used sd instead of sr for error bars."; protein.plot.se <- FALSE}


show.sd.data[is.infinite(show.sd.data)] <- 0
show.sd.data[is.na(show.sd.data)] <- 0


}else{show.sd.data <- matrix(0,dim(row.plot.data))}

if(is.matrix(.aov.new)){
.aov.data		<- as.numeric(as.character(.aov.new[hz.merge.control(gsub(" ","",.aov.new[,1]),gsub(" ","",rownames(row.plot.data))),2]))
}else{
	.aov.data <- NA
}




if(is.null(sub.info)|!is.matrix(sub.info)){
	sub.info<- rep("",dim(row.plot.data)[1])
	
	
}else{sub.info <- sub.info[,3]}

if(!exists(".design")){	.design <- NULL}
if(length(unique(.design$Time)) < 2 & gui.input$exp.design != ""){
	#gui.input$barpl <- FALSE
	#tkmessageBox(message = "Time column in ED contains only one entry. Time groups switched off!",icon = "warning")
	gui.input$time.grouped <- FALSE
}

if(dim(.data2$x)[2] == 2 & all(gui.input$time.grouped,gui.input$barpl)){
	gui.input$time.grouped <- FALSE
}
if(gui.input$exp.design== ""){
	gui.input$time.grouped <- FALSE
}

if(!gui.input$time.grouped){
gui.input$lineplot.beside  <- FALSE
}



print("started plotting")
#stop()
	hz.row.plot(	x = row.plot.data,
					show.sd = show.sd.data,
					sd.rel = TRUE,
					
					.aov = .aov.data 
					,p.v = gui.input$p.value,ui = ui
					,sub =sub.info ,
					x.ylab = y.lab.input,
					plot.col =.col,	
					inf.m = prot.n,
					graphic.type = gui.input$graphic.type,
					plot.type = "b",
					barpl = gui.input$barpl,
					time.groups = gui.input$time.grouped,
					group.barplot = F,
					x.xlab = gui.input$x.xlab,
					.design = .design,
					lineplot.beside = gui.input$lineplot.beside,
					gui.input = gui.input,prog.max=prog.max,ratio.prog= ratio.prog,pb=pb,hz.exp.des.parse.data2=hz.exp.des.parse.data2,colorblind.set=colorblind.set,.col=.col
					)
print("finished plotting")
return(list())
}
