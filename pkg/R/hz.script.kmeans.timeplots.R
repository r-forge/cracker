hz.script.kmeans.timeplots <-
function(kmeans.cluster.output,i,k.data,.design,gui.input, y.lab.input,main.temp,colorblind.set,.col,prog.max,pb,ui){
	if(!exists("kmeans.mean")){	kmeans.mean <- c()}
if(gui.input$time.grouped){
		print("time")

temp.i 	<- as.matrix(kmeans.cluster.output[,1][kmeans.cluster.output[,1] == i])
		temp.merge	<-	merge(as.matrix(temp.i),k.data,by = 0)
		temp.merge 	<-  temp.merge[,-c(1,2)]
		norm.k.x <- max(nchar(colnames(temp.merge)))
		if(is.vector(temp.merge)){temp.merge <- as.matrix(temp.merge)}
	
	#if(norm.k.x > 8){norm.k.x <- (norm.k.x-8)/2}

	par(oma = c(1.5,0.1,0.1,0.1),mai = c(1.5,1.5,1,1))
	
	if(gui.input$raw){
		.design.used <- cbind(.design$Alternative.name,as.character(.design$Group),.design$Time)
	}else{
		.design.used <- unique(cbind(.design$Experiment,as.character(.design$Group),.design$Time))
	}
	
	
	correct.exp.des <- hz.merge.control(.design.used[,1],gsub(" ","",colnames(temp.merge)))
	correct.exp.des <- correct.exp.des[!is.na(correct.exp.des)]
	.design.used 	<- .design.used[correct.exp.des,]
	
	list.groups.y <- list()
	list.groups.x <- list()
	list.group	 <- c()
	template.groups <- as.numeric(unique(.design.used[,3]))
	for(groups.i in 1:length(unique(.design.used[,2]))){
		groups.temp.i 		<- .design.used[.design.used[,2] == unique(.design.used[,2])[groups.i],]		
		if(is.vector(groups.temp.i)){groups.temp.i <- t(as.matrix(groups.temp.i))}
		temp.merge.merge 	<- hz.merge.control(gsub(" ","",colnames(temp.merge)),groups.temp.i[,1])
		
		
		
		#stop()
		temp.merge.merge 	<- temp.merge.merge[!is.na(temp.merge.merge)]
		
		list.groups.y[[groups.i]]		<- temp.merge[,temp.merge.merge]-min(temp.merge[,temp.merge.merge],na.rm = T)+1
		list.groups.x[[groups.i]]		<- as.numeric(groups.temp.i[,3])
		if(gui.input$barpl){
		temp.list.g <- temp.merge[,temp.merge.merge]
		if(is.vector(temp.list.g)){temp.list.g <- t(as.matrix(temp.list.g))}
		list.groups.x[[groups.i]]	<-  c(1:dim(temp.list.g)[2])
			
		}

		
		list.group[groups.i]			<- unique(.design.used[,2])[groups.i]
	}

if(gui.input$barpl){
	mat.axes <- FALSE
	mat.xlab <- ""

}else{
	mat.axes <- TRUE
	mat.xlab <- gui.input$x.xlab
	
}
matplot(1,type = "n",xlim =range(list.groups.x,na.rm = T),ylim = .range.y<- range(list.groups.y,na.rm = T),ylab = y.lab.input
, main = main.temp
,xlab = mat.xlab,
axes = mat.axes	,pch = 16,lwd = 1.5 ,cex.lab = 0.6,frame = F
)



grid()
temp.mean.collect.x <- c()
temp.mean.collect.y <- c()

for(groups.plot in 1:length(unique(.design.used[,2]))){
	print(groups.plot)
		groups.plot.x <- list.groups.x[[groups.plot]]
	
	groups.plot.y <- list.groups.y[[groups.plot]]
	
	#groups.plot.y <- t(apply(groups.plot.y,1,function(x){x-min(x,na.rm = TRUE)+1})
#)
	
	if(is.vector(groups.plot.y)){groups.plot.y <- t(as.matrix(groups.plot.y))}
	####
	# color
	###
temp.mean <- apply(groups.plot.y,2,function(x){x<-as.numeric(x);x<- median(x,na.rm = TRUE)})


sum.test <- apply(groups.plot.y,1,function(x){sum(abs(x - temp.mean),na.rm = TRUE)})
	sum.test <- round(sum.test*10)
	order.test 	<- rev(order(sum.test))
	sum.test	<- (sum.test)
	range.test 	<- max(sum.test)
#ramp <- colorRampPalette(c("red","orange","yellow","green","cyan","blue","purple"))
#test.ramp <- ramp(range.test)


if(colorblind.set){
				test.ramp		<- rev(colorRampPalette(c("darkgrey" ,unique(.col)[groups.plot]))(max(unique(sum.test))))
			}else{
test.ramp <- rainbow(max(unique(sum.test)),v = 0.88,alpha = 0.7,end = 0.9)
				
		}


if(!gui.input$color.plot){
	ramp <- colorRampPalette(c(unique(.col)[groups.plot],"darkgrey"))
	test.ramp <- ramp(length(unique(sum.test)))
}

if(!gui.input$color.plots){
	ramp <- colorRampPalette(c(colors()[275],"white"))
	test.ramp <- ramp(length(unique(sum.test)))

}

#stop()#test.ramp <- ramp(length(unique(sum.test)))
sum.test <- test.ramp[sum.test]
sum.test[is.na(sum.test)] <- "#808080"




	matlines(groups.plot.x ,t(groups.plot.y),type = "pl",pch = 1 ,col = sum.test,lty = 1,lwd = 2)
temp.control.vec <- hz.merge.control(groups.plot.x ,template.groups)

temp.mean.collect.y <- rbind(temp.mean.collect.y ,temp.mean[temp.control.vec])
temp.mean.collect.x <- rbind(temp.mean.collect.x ,groups.plot.x[temp.control.vec])

#	points(temp.mean,type = "l",col = "red",lwd = 5)	
#	points(temp.mean,type = "l",col = unique(.col)[groups.plot],lwd = 3)	

}

if(gui.input$color.plots){
.col.mean.temp <- unique(.col)
.col.mean.temp.bg <- "white"
}else{
.col.mean.temp <- c("white","black")	
.col.mean.temp.bg <- c("black","white")
	
}


matlines(t(temp.mean.collect.x),t(temp.mean.collect.y),lwd = 8,col = .col.mean.temp.bg,lty = 1)
matlines(t(temp.mean.collect.x),t(temp.mean.collect.y),lwd = 5,col = .col.mean.temp,lty = 1)
print(.col.mean.temp)
if(gui.input$barpl){


	axis(2)
	
		for(test.i in 1: length(unique(.design.used[,2]))){
			
			
			.names.test.i <- .design.used[.design.used[,2] == unique(.design.used[,2])[test.i],1]
			if(test.i == 1){
				padj.test.i <- test.i
			}else{
				padj.test.i <- (test.i+0.2* test.i)
			}
			print(.names.test.i)
			
			axis(1,at = c(1:length(.names.test.i)),labels = .names.test.i,las = 2, col.axis = .col.mean.temp[test.i], padj = padj.test.i)
		}
		kmeans.mean <- rbind(kmeans.mean,temp.mean)
#try(		rownames(kmeans.mean)[i] <- paste("Cluster",i))
		#axis(3,at = c(1:dim(temp.merge)[2]),labels = colnames(temp.merge),las = 2)

}
		legend("topright",unique(.design.used[,2]),fill= unique(.col.mean.temp) ,bg = "#FFFFFF99",cex = 0.5,xpd = T,inset = c(-0.05,-0.15))

}
}
