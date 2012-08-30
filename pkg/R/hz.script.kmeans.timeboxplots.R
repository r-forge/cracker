hz.script.kmeans.timeboxplots <-
function(kmeans.cluster.output,i,k.data,.design,gui.input, y.lab.input,colorblind.set,.col,prog.max,pb,ui,kmeans.list,kmeans.at,kmeans.col){
	
	print("box")
temp.i 	<- as.matrix(kmeans.cluster.output[,1][kmeans.cluster.output[,1] == i])
		temp.merge	<-	merge(as.matrix(temp.i),k.data,by = 0)
		temp.merge 	<-  temp.merge[,-c(1,2)]
		norm.k.x <- max(nchar(colnames(temp.merge)))
		if(is.vector(temp.merge)){temp.merge <- as.matrix(temp.merge)}
	
	if(norm.k.x > 3){norm.k.x <- 2+(norm.k.x-3)*0.1}else{
		norm.k.x <- 2
	}

	par(oma = c(norm.k.x,0.1,0.1,0.1),mai = c(norm.k.x*0.4,1.5,1,1))
if(gui.input$raw){
		.design.used <- cbind(.design$Alternative.name,as.character(.design$Group),.design$Time)
	}else{
		.design.used <- unique(cbind(.design$Experiment,as.character(.design$Group),.design$Time))
	}

order.bp <- c()
.gap <- 1
at.temp <- seq(.gap*1.1,length.out = dim(temp.merge)[2])#1:dim(temp.merge)[2]


sum.it <- 0
for(temp.bp in 1:length(unique(.design.used[,2]))){

	temp.bp.i 		<- .design.used[.design.used[,2] == unique(.design.used[,2])[temp.bp],1]
	temp.bp.order  	<- hz.merge.control(gsub(" $","",colnames(temp.merge)),temp.bp.i)
	order.bp 	 <- c(order.bp,temp.bp.order)
	
	sum.it <- sum(sum.it +length(temp.bp.order)+.gap)
	
	print(sum.it)
		at.temp[at.temp > sum.it] <- at.temp[at.temp > sum.it]+.gap
		
	

}
temp.merge 		<- temp.merge[order.bp]
.col.temp		<- .col[order.bp]
.col.temp[is.na(.col.temp)] <- "grey"


	plot(1 ,1,type = "n",main = paste("cluster",i,"\n", length(temp.i),"proteins"),axes = FALSE,xlab = "",ylab = y.lab.input,col = .col.temp,cex.lab = 0.6,xlim = c(min(at.temp)*0.9,max(at.temp))*1.01,ylim = range(temp.merge,na.rm = T)+c(0,range(temp.merge,na.rm = T)[2]*0.05))

	grid()
	kmeans.list[[i]] <- temp.merge
	kmeans.at[[i]] 	 <- at.temp
	kmeans.col[[i]]  <- .col.temp

temp	<- boxplot(			temp.merge ,	type = "n",
							col = .col.temp,
							notch = FALSE,
							border = TRUE,
							lwd = 1.5,
							boxwex = 0.6,
							cex.lab = 0.6,
							add = T,
							axes = F,
							at = at.temp)
	
		axis(2)
	temp.mean <- apply(temp.merge,2,function(x){x<-as.numeric(x);x<- median(x,na.rm = TRUE)})

		for(.col.temp.i in unique(.col.temp)){
			print(.col.temp.i)
			.col.temp.temp.i <- .col.temp == .col.temp.i
			axis(1,at = at.temp[.col.temp.temp.i],labels = colnames(temp.merge)[.col.temp.temp.i],las = 2,col.axis = .col.temp.i)
		
		points(at.temp[.col.temp.temp.i],temp.mean[.col.temp.temp.i],type = "l",col = "grey",lwd = 4)		
		points(at.temp[.col.temp.temp.i],temp.mean[.col.temp.temp.i],type = "l",col = .col.temp.i,lwd = 2)
		
		}
		
		
		legend("topright",unique(.design.used[,2]),fill= unique(.col.temp) ,bg = "#FFFFFF99",cex = 0.7,xpd = T,inset = c(-0.12,-0.15), xjust = 1,box.col = "transparent",title = "legend")
return(list(kmeans.col= kmeans.col,kmeans.at = kmeans.at))
}
