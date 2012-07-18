hz.script.kmeans <-
function(.data2,gui.input,.design, y.lab.input,colorblind.set,color.blind,.col,prog.max,pb,ui){
	
		if(!exists("ratio.prog")){ratio.prog <- 1000}

######
# kmeans
######
centers <- gui.input$centers
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4, label=paste( "Kmeans calculation"))))

		while(1==1&pb.check == "try-error"){
				print("Warning: User closed window!")
pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4, label=paste( "Kmeans calculation"))))
		}
	##############
height.set <- 9		
	
if(gui.input$norm.method == "median"){norm.temp <- "mean"}else{norm.temp <- gui.input$norm.method}

	k.data 			<- .data2$x
	#rownames(k.data) <- paste("lakdliawedliawehdliha",rownames(k.data))
#k.data	<- cbind(k.data,k.data)	
	k.data.rows 	<- rownames(k.data)
	k.data 			<- apply(k.data,2,as.numeric)
	rownames(k.data) <- k.data.rows
	if(1==1){
	k.data <- hz.norm(k.data,norm =norm.temp)$x
		
	if(all(is.na(as.numeric(k.data)))){
			k.data <- hz.norm(k.data,1,norm = "median")$x
			print("Warning: Applied median normalization for kmeans.")		
	}
	}else{
	
		if(!as.logical(gui.input$row.norm)){

			k.data <- hz.norm(k.data,norm =norm.temp)$x
		
		if(all(is.na(as.numeric(k.data)))){
			k.data <- hz.norm(k.data,1,norm = "median")$x
			print("Warning: Applied median normalization for kmeans.")		
		}
		}
	}

	
	k.min <- min(k.data,na.rm = T)
	if(k.min < 0 ){
		print("minimal value for negative data")
		k.data <- k.data+ abs(k.min)+abs(k.min)*0.005
		
	}
	


	if(dim(k.data)[1] > 2 & dim(k.data)[2] > 1){
	k.data.NA <- k.data
	#if(!n15.log2){
	#k.data.NA[is.na(k.data.NA)] <- 0}else{
	k.data.NA[is.na(k.data.NA)]	 <- 2*min(k.data.NA,na.rm = TRUE)
#	}
	# kmeans based on Correlation?
	centers.check <- is.na(as.numeric(centers))
	if(centers == "auto" & centers.check& 1==1){
		cutree.cut <- T
		if(is.na(as.numeric(centers))& centers !="auto"){
			centers <- 2
		}
	# hartigan.criteria from Detlef Groth	
		hartigan.criteria = function (data) { 
			maxk=dim(data)[1]-1; 
			nr=nrow(data); 
			for (i in 2:maxk) {
				clu.A=kmeans(data,i) ; 		
				clu.B=kmeans(data,i+1); 
				sum.A=sum(clu.A$withinss); 		
				sum.B=sum(clu.B$withinss); 
				val = (sum.A/sum.B-1)*(nr-i-1) 
				#print(val)
				if (val < 10) { return(i) ;
				}
			}
		}
		
		
	#centers <- hartigan.criteria(t(k.data.NA))
	}else{
		cutree.cut <- F
	}
		if(centers < 2 | dim(k.data)[1] <= centers){centers = 2}
	
	
	
	.min 		<- min(k.data.NA,na.rm = T)
	use.correlation.matrix <- F
	if(use.correlation.matrix){
	dist.object <- k.data.NA

	dist.object[is.na(dist.object)] <- .min - abs(.min)*0.1
	
	dist.object <- cor(t(dist.object),method = gui.input$do.cor,use = "everything")
	
	dist.object[is.na(dist.object)] <- 0
	}else{
	dist.object <- k.data.NA		
	}
	#dist.object <- dist(k.data.NA,method = "manhattan")

	
	error.rep = "try-error"
while(error.rep == "try-error"){
	if(gui.input$hclust.groups){
	
	hclust.test <- hclust(dist(plot.clustering))
		if(cutree.cut){
			error.rep  <- class(try(kmeans.cluster.output <- cutree(plot.clustering,k=centers)))
			
		}else{
			error.rep  <- class(try(kmeans.cluster.output <- cutree(plot.clustering,h=max(plot.clustering$height)*0.666)))
		centers <- max(kmeans.cluster.output)
			
		}
	
	
	kmeans.cluster.output  <- as.matrix(kmeans.cluster.output)
	#centers <- max(kmeans.cluster.output)
	}else{			
	error.rep <- class(try(	kmeans.cluster <- kmeans(dist.object,as.numeric(centers),iter.max = 200)
	))

	kmeans.cluster.output  <-as.matrix(kmeans.cluster$cluster)
	}
	

	
	if(error.rep == "try-error"){
		centers <- centers -1
	}	
	if(centers < 2){
		stop()
	}
}


	if(gui.input$hclust.groups){
		cluster.name <- paste("hclust-",centers,sep = "")
		
	}else{
		cluster.name <- paste("kmeans-",centers,sep = "")
	}
	
	colnames(kmeans.cluster.output) <- "cluster"
	write.csv(kmeans.cluster.output,paste(cluster.name,".csv",sep = ""))
	
	bp.width <- 9
	if(dim(k.data)[2] > 15){
		bp.width <- bp.width+ (dim(k.data)[2]-15)* 0.4
	}
	.wd <- getwd()

	if(gui.input$graphic.type == "pdf"){
		pdf(paste(cluster.name,"-lp.pdf",sep = ""),height = height.set, width = bp.width,pointsize= 18)
	}else{
		dir.create(.wd.set <- paste(cluster.name,"-eps",sep = ""))
		setwd(.wd.set)
	}
	
	
	kmeans.mean <- c()
	kmeans.list <- c()
	kmeans.at	<- c()
	kmeans.col	<- c()
	for(i in 1:centers){
		
	if(gui.input$graphic.type == "eps"){

		postscript(paste(cluster.name,"-number-",i,"-lp.eps",sep = ""), height = height.set,width = bp.width, paper = "special",onefile = FALSE,horizontal = FALSE,pointsize= 18)
	}
		
		temp.i 	<- as.matrix(kmeans.cluster.output[,1][kmeans.cluster.output[,1] == i])
		temp.merge	<-	merge(as.matrix(temp.i),k.data,by = 0)
		temp.merge 	<-  temp.merge[,-c(1,2)]
		norm.k.x <- max(nchar(colnames(temp.merge)))
		kmeans.list[[i]] <- as.matrix(temp.merge)
		kmeans.at[[i]] 	 <- 1:length(colnames(temp.merge))
	
	if(norm.k.x > 8){norm.k.x <- (norm.k.x-8)/2}
	
	par(oma = c(norm.k.x,0.1,0.1,0.1))

	temp.mean <- apply(temp.merge,2,function(x){x<-as.numeric(x);x<- median(x,na.rm = TRUE)})

	sum.test <- apply(temp.merge,1,function(x){sum(abs(x - temp.mean),na.rm = TRUE)})
	sum.test <- round(sum.test*10)
	order.test 	<- rev(order(sum.test))
	sum.test	<- (sum.test)
	range.test 	<- max(sum.test)
#ramp <- colorRampPalette(c("red","orange","yellow","green","cyan","blue","purple"))
#test.ramp <- ramp(range.test)


if(colorblind.set){
				test.ramp		<- rev(colorRampPalette(c(unlist(color.blind)[c(7,6,4,3,2,1)]))(max(unique(sum.test))))
			}else{
test.ramp <- rainbow(max(unique(sum.test)),v = 0.88,alpha = 0.7,end = 0.9)
				
		}


if(!gui.input$color.plot){
	ramp <- colorRampPalette(c(colors()[338],colors()[276]))
	test.ramp <- ramp(length(unique(sum.test)))
}

	

#stop()#test.ramp <- ramp(length(unique(sum.test)))
sum.test <- test.ramp[sum.test]
sum.test[is.na(sum.test)] <- "#808080"
	
plot.data <-  t(temp.merge[order.test,])
if(all(is.na(plot.data))){
temp.merge[is.na(temp.merge)] <- 0	
}	
#gui.input$time.grouped <- T
if(gui.input$time.grouped&gui.input$exp.design != ""){
	main.temp <- paste("cluster",i,"\n", length(temp.i),"proteins")
	hz.script.kmeans.timeplots.return <- hz.script.kmeans.timeplots(kmeans.cluster.output,i,k.data,.design,gui.input, y.lab.input,main.temp,colorblind.set,.col,prog.max,pb,ui)
	kmeans.col  <- hz.script.kmeans.timeplots.return$kmeans.col
	
}else{
	print(dim(t(temp.merge[order.test,])))
	matplot(t(temp.merge[order.test,]) ,type = "n",main = paste("cluster",i,"\n", length(temp.i),"proteins"),axes = FALSE,xlab = "",ylab = y.lab.input,col = sum.test[order.test],pch = 16,lwd = 1.5 ,cex.lab = 0.6#,ylim = range(k.data.NA)
	)
	grid(lwd = 2)
	matlines(t(temp.merge[order.test,]) ,type = "l",col = sum.test[order.test],pch = 16,lwd = 2 ,cex.lab = 0.6#,ylim = range(k.data.NA)
	)

	
	
	axis(2)
		axis(1,at = c(1:dim(temp.merge)[2]),labels = colnames(temp.merge),las = 2)
		
		kmeans.mean <- rbind(kmeans.mean,temp.mean)
		rownames(kmeans.mean)[i] <- paste("Cluster",i)
		points(temp.mean,type = "l",col = "red",lwd = 3)		
		
		
		}
}

	graphics.off()#dev.off()	
	
	
	bp.width <- 7
	if(dim(k.data)[2] > 15){
		bp.width <- bp.width+ (dim(k.data)[2]-15)* 0.2
		
	}
	if(gui.input$graphic.type == "pdf"){
		pdf(paste(cluster.name,"-bp.pdf",sep = ""), height = height.set,width =bp.width,pointsize= 18)
	}
	
	for(i in 1:centers){
	
	if(gui.input$graphic.type == "eps"){
		postscript(paste(cluster.name,"-number-",i,"-bp.eps",sep = ""), height = height.set,width =bp.width,paper = "special",onefile = FALSE,horizontal = FALSE)

	}
	
		if(gui.input$time.grouped&gui.input$exp.design != ""){
	main.temp <- paste("cluster",i,"\n", length(temp.i),"proteins")
	hz.script.kmeans.timeplots.return <- hz.script.kmeans.timeboxplots(kmeans.cluster.output,i,k.data,.design,gui.input, y.lab.input,colorblind.set,.col,prog.max,pb,ui,kmeans.list,kmeans.at,kmeans.col)
	kmeans.col  <- hz.script.kmeans.timeplots.return$kmeans.col
	
	}else{

		
		temp.i 	<- as.matrix(kmeans.cluster.output[,1][kmeans.cluster.output[,1] == i]) 
		temp.merge	<-	merge(as.matrix(temp.i),k.data,by = 0)
		temp.merge 	<-  temp.merge[,-c(1,2)]
		#temp.merge <- temp.merge[,c(1,2,3,4,5,6,1,2,3,4,1,2,3,4)]
		norm.k.x <- max(nchar(colnames(temp.merge)))
	
	if(norm.k.x > 8){norm.k.x <- (norm.k.x-8)/2}

	par(oma = c(norm.k.x,0.1,0.1,0.1))
	
	

	
	
	
#	boxplot(temp.merge ,type = "n",main = paste("cluster",i,"\n", length(temp.i),"proteins"),axes = FALSE,xlab = "",ylab = y.lab.input,col = .col,notch = FALSE,border = TRUE,lwd = 1.5,boxwex = 0.6,cex.lab = 0.6)
	plot(1 ,1,type = "n",main = paste("cluster",i,"\n", length(temp.i),"proteins"),axes = FALSE,xlab = "",ylab = y.lab.input,col = .col,cex.lab = 0.6,xlim = width.2<- range(0.85,(dim(temp.merge)[2]+0.15)),ylim = range(temp.merge,na.rm = T))

	grid(lwd = 2)
	boxplot(temp.merge ,type = "l",col = .col,notch = FALSE,border = TRUE,lwd = 1.5,boxwex = (width.2[2]-width.2[1])/dim(temp.merge)[2]*0.4,cex.lab = 0.6,add = T,axes = F)
	
	
	
		axis(2)
		axis(1,at = c(1:dim(temp.merge)[2]),labels = colnames(temp.merge),las = 2)
		
		temp.mean <- apply(temp.merge,2,function(x){x<-as.numeric(x);x<- median(x,na.rm = TRUE)})
			
		points(temp.mean,type = "l",col = "grey",lwd = 4)		
		points(temp.mean,type = "l",col = "red",lwd = 2)		}
	
	}
graphics.off()
	}
	#============= end kmea
	print("control point")




setwd(.wd)

return(list(kmeans.cluster.output = kmeans.cluster.output, kmeans.col = kmeans.col, kmeans.at = kmeans.at, kmeans.list= kmeans.list))
}
