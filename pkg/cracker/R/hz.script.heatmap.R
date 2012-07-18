hz.script.heatmap <-
function(.data2,gui.input,prog.max,pb,ui, ratio.prog){
	
	if(!exists("ratio.prog")){ratio.prog <- 1000}
	
	
library("pcaMethods")
.data2$x[sapply(.data2$x, is.na)] <- NA
pca.rows 	<- 	rownames(.data2$x)
data 		<-	apply(.data2$x,2,as.numeric)

data.sum 	<-	apply(data,1,function(x){sum(x,na.rm = TRUE)})
# remove all rows, cols without entry
data 		<- data[data.sum !=0,]
pca.rows	<- pca.rows[data.sum !=0]
if(!gui.input$n15.log2){
data.log 	<- log2(data)}else{data.log <- data}
data.log.inf <- hz.inf.replace(data.log,0.1)
data.log <- data.log.inf$x



data.log[data.log == "NaN"] 		<- NA
#data.log			<- hz.norm(data.log,1,norm = "z")$x


NA.vec 	<- which(apply(as.matrix(data.log),1,function(x){all(is.na(x))}))

if(length(NA.vec) > 0){
	cat("Matrix contains NA rows. Excluded:")
	data.log <- data.log[-NA.vec,]
	pca.rows <- pca.rows[-NA.vec]
}

rownames(data.log) <- pca.rows
NA.vec <- which(apply(as.matrix(data.log),2,function(x){all(is.na(x))}))
if(length(NA.vec) > 0){
	cat("Matrix contains NA cols. Excluded:")
	data.log <- data.log[,-NA.vec]
}


if(dim(data.log)[2] > 1){

	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*2, label=paste( "Heatmap calculation..."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*2, label=paste( "Heatmap calculation..."))))
		}
	##############
		




#

}
return(list(data.log = data.log,data.sum <- data.sum))
}
