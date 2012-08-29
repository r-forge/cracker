hz.script.graph <-
function(.data2,gui.input,prog.max,pb,ui){
cor.threshold 	<- 0.9
	

.rows <- rownames(.data2$x)
.data2$x <- apply(.data2$x,2,as.numeric)
rownames(.data2$x) <- .rows


zero.ex <- apply(.data2$x,1,function(x){all(x == 0)})

.min <- min(t(.data2$x[!zero.ex,]),na.rm = T)

.data2$x[is.na(.data2$x)] <- .min + 0.1*.min

print("Create correlation matrix/list")
data.graph2 <- cor(t(.data2$x[!zero.ex,]),method = gui.input$do.cor,use = "everything")



matrix.to.list <- function(x){
	.init <- c()
	for(i in 2:dim(x)[1]){
		temp.i <- x[i,1:(i-1)]
		if(length(temp.i)== 1){
			temp.i <- cbind(rownames(x)[i],colnames(x)[i-1],temp.i)
		}else{
			temp.i <- cbind(rep(rownames(x)[i],length(temp.i)),names	(temp.i),temp.i)
		}
		.init <- rbind(.init,temp.i)
		}
	return(.init)
	}


#library(igraph)

.graph.matrix 	<- 	matrix.to.list(data.graph2)
.graph.matrix <- .graph.matrix[,c(1,3,2)]
colnames(.graph.matrix) <- c("protein 1","correlation.coefficient","protein 2")
rownames(.graph.matrix) <- NULL
.graph.matrix <- as.data.frame(.graph.matrix)
write.table(.graph.matrix,paste("network",gui.input$do.cor,"-cor.txt"),sep = "\t",row.names = F)

#.graph.matrix  	<- apply(.graph.matrix,2,as.numeric)


#NA.exclude	<- apply(.graph.matrix,1,function(x){
#	x <- grep("NA id",as.character(x))
#	if(length(x) >0){ return(FALSE)
#	}else{return(TRUE)}	
#	})

#.graph.matrix <- .graph.matrix[NA.exclude,]

.graph.matrix.signi <- .graph.matrix[!is.na(.graph.matrix[,2]),]

if(is.vector(.graph.matrix.signi)){
	.graph.matrix.signi <- t(as.matrix(.graph.matrix.signi))
	
}


.graph.matrix.signi <- .graph.matrix.signi[!as.character(.graph.matrix.signi[,1] )== as.character(.graph.matrix.signi[,3]),]
if(is.vector(.graph.matrix.signi)){
	.graph.matrix.signi <- t(as.matrix(.graph.matrix.signi))
	
}
.graph.matrix.signi <- .graph.matrix[abs(as.numeric(.graph.matrix.signi[,2])) > cor.threshold,]

if(length(dim(.graph.matrix.signi)) > 0){
colnames(.graph.matrix.signi) <- c("protein 1","correlation.coefficient","protein 2")
rownames(.graph.matrix.signi) <- NULL

write.table(.graph.matrix.signi,paste("network",gui.input$do.cor,"-cor-threshold.txt"),sep = "\t",row.names = F)


#stop("end of script")
if(gui.input$do.network&1==0){

if(dim(.graph.matrix.signi)[1] > 1){

.color <- as.numeric(.graph.matrix.signi[,3])
.color[.color > cor.threshold] <- 1
.color[.color < cor.threshold] <- -1
.color[.color ==-1] <- 2
.color[is.na(.color)] <- 9
  
 print("Creating Network")                    
g <- graph.data.frame(as.data.frame(.graph.matrix.signi), directed=FALSE)   

# sets names for label
V(g)$label <- V(g)$name
# sets layout
g$layout <- layout.spring
# sets color
E(g)$color <- .color
# 
 print("Plotting Network")                    

	.size = 10+0.01*dim(data.graph2)[1]
	if(gui.input$graphic.type == "pdf"){
		pdf("networkdata.graph.pdf",width = .size,height =.size)
	}else{
		postscript("networkdata.graph.eps",width = .size,height =.size, paper = "special",onefile = FALSE,horizontal = FALSE)	
	}	
	plot(g	#,vertex.shape = "none"
			,vertex.label.cex = 0.4
			,vertex.label.dist = 0
			,vertex.color = "transparent"
			,edge.width = c(0.2)
			,edge.arrow.size = 0	
	)
	dev.off()
	
#
d <- degree(g,mode = "in")
#dd <- degree.distribution(g, mode="in", cumulative=TRUE)

degree.matrix <- cbind(V(g)$name,d)
write.csv(degree.matrix,"degree-network.csv")
}
#alpha <- power.law.fit(d, xmin=20)

#plot(dd)
}
}else{
	
	cat("No entries left after setting threshold.\nAborted Creating Network.")
}
}
