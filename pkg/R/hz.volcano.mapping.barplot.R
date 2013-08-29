hz.volcano.mapping.barplot <- 
function(x,gui.input){

	temp <- getwd()
test 		<- list()
.wd <- getwd()

test$path	<- .wd #gui.input$path.data
dir.create(temp.path <- paste(test$path,"volcano-barplot-mappings",sep = "/"))

if(1==1){

for(mapping.type in  grep("mapping.",colnames(x))){


for(i in 1:length(unique(x$samples))){
temp.test <- x[x$samples== unique(x$samples)[i],]


test$x1 	<- temp.test[as.numeric(as.character(temp.test$ratio)) > 0,mapping.type] 
test$x2 	<- temp.test[as.numeric(as.character(temp.test$ratio)) < 0,mapping.type]
test$sample 	<- gsub(" ","",unlist(strsplit(unique(as.character(x$samples))[i]," /",fixed = T)))[1]
test$reference 	<- gsub(" ","",unlist(strsplit(unique(as.character(x$samples))[i]," /",fixed = T)))[2]
test$exclude	<- 1
temp.name <- paste(test$sample,"-vs-",test$reference,sep = "")

exclude.count 	<- as.numeric(test$exclude)

objects.plot	<-"full"
pdf.width		<- 5

p.value.setting <- 0.05
#setwd(test$path)
#wd.dir <- paste("fishy-results","-ex",exclude.count,"-", make.names(test$sample),Sys.Date(),sep = "")
#dir.create(wd.dir)
#setwd(wd.dir)

hz.merge.control <- function(input,ref){
#input 	<- sub.info[,2]
#ref 	<- rownames(row.plot.data)

	input <- as.character(input)
	input <- cbind(input,c(1:length(input)))
	
	ref <- as.character(ref)
	ref <- cbind(ref,c(1:length(ref)))

	.merge	<- merge(as.matrix(ref),as.matrix(input),by = 1,all.x = T)
	.merge 	<- .merge[order(as.numeric(as.character(.merge[,2]))),]
	return(as.numeric(as.character(.merge[,3])))
	
	
}


try(x1 <- unlist(strsplit(as.character(unlist(test$x1)),"\n|\r")))
try(x2 <- unlist(strsplit(as.character(unlist(test$x2)),"\n|\r")))
try(x1 <- unlist(strsplit(unlist(x1),"|",fixed = TRUE)))
try(x2 <- unlist(strsplit(unlist(x2),"|",fixed = TRUE)))
if(!exists("x1")){x1 <- ""}
if(!exists("x2")){x2 <- ""}


x2 <- x2[!nchar(x2)== 0]
x1 <- x1[!nchar(x1)== 0]

if(length(x2) == 0){
	x2 <- ""
}

x2 <- aggregate(x2,list(x2),length)
x1 <- aggregate(x1,list(x1),length)

rm(merge.data)
merge.data <- merge(as.matrix(x1),as.matrix(x2),all = TRUE,by = 1)

merge.data <- apply(merge.data,2,as.character)
merge.data[is.na(merge.data)] <- "0"
print("hui")
if(length(merge.data) > 0){
exclude.count.vec <- as.numeric(merge.data[,2]) <= exclude.count&as.numeric(merge.data[,3]) <= exclude.count

excluded.strings <- merge.data[exclude.count.vec,1]
merge.data <- merge.data[!exclude.count.vec,]
}
if(!is.vector(merge.data)& length(merge.data) > 0){

if(!is.matrix(merge.data)){
	merge.data <- t(as.matrix((merge.data)))
}
merge.data <- apply(merge.data,2,as.character)
if(!is.matrix(merge.data)){
	merge.data <- t(as.matrix((merge.data)))
}
merge.data <- merge.data[order(as.character(merge.data[,1])),]
if(!is.matrix(merge.data)){
	merge.data <- t(as.matrix((merge.data)))
}
#stop()


if(objects.plot == "full"){
	objects.plot <- dim(merge.data)[1]
	as.numeric(objects.plot)
	
}
if(is.na(objects.plot)){objects.plot <- 20}

if(as.numeric(objects.plot)> 20){
	height.pdf <- 6 + (objects.plot-20)*0.1 
	pdf.cor <- (objects.plot-20)*0.2
}else{
	height.pdf <- 8
	pdf.cor <- 0
}


#  stop()


if(length(excluded.strings)> 0){
try(x1.control <- hz.merge.control(x1[,1],excluded.strings))
try(x1 <- x1[-x1.control[!is.na(x1.control)],])
try(x2.control <- hz.merge.control(x2[,1],excluded.strings))
try(x2 <- x2[-x2.control[!is.na(x2.control)],])
}

fisher.results <- c()


#pb <- tkProgressBar("test progress bar", "Some information in %",0, length(unique(x1[,1])), 0)

	d <- 1

is.reference<- FALSE

#stop()


if(is.reference){
	templates<- unique(x1[,1])	
}else{
	templates <- unique(c(x1[,1],x2[,1]))
} 

# applying fisher exact test
for(i in templates){
		grep.i	<- merge.data[,1] == i
	
		grep.i[is.na(grep.i)] <- FALSE

	temp.i 	<- as.vector(as.matrix(merge.data[grep.i,]))
		d1		<- as.numeric(as.character(temp.i[2]))
		d2		<- as.numeric(sum(as.numeric(merge.data[,2]),na.rm = TRUE))	
		d3		<- as.numeric(as.character(temp.i[3]))
		d4		<- as.numeric(sum(as.numeric(merge.data[,3]),na.rm = TRUE))	
		
		    #setTkProgressBar(pb, d)
		d = d+1
		
		d.ratio <- (d1/d2)/(d3/d4)
		
		d.con.table 	<- rbind(c(d1,abs(d1-d2)),c(d3,abs(d3-d4)))
		d.con.table[is.na(d.con.table)] <- 0
		d.fisher		<- fisher.test(d.con.table)$p.value
		data.barplot		<- c(i,d1,d2,d3,d4,d.fisher)
		fisher.results 	<- rbind(fisher.results, data.barplot)
	
}

fisher.results 				<- cbind(fisher.results,p.adjust(as.numeric(fisher.results[,6]),"BH") )
colnames(fisher.results) 	<- c(
							test$sample,
							paste("n(",test$sample,")"),
							paste("n(",test$sample,") - n(",test$sample,"data.barplot)"),
							paste("n(",test$reference,")"),
							paste("n(",test$reference,")-n(",test$reference,"data.barplot)"),
							"p-value fisher",
							"corrected p-value BH")
print("owiefj")

data.barplot.init <- fisher.results[,c(2,4)]
if(!is.matrix(data.barplot.init)){data.barplot.init <- t(as.matrix(data.barplot.init))}
data.barplot <- apply(data.barplot.init,2,as.numeric)#print(test)
if(!is.matrix(data.barplot)){data.barplot <- t(as.matrix(data.barplot))}

rownames(data.barplot) <- fisher.results[,1]

data.barplot[is.na(data.barplot)] <- 0
.nchar.cor <- nchar(merge.data[,1])
.nchar.cor.val <- 1


if(max(.nchar.cor) > 8){
	.val <- max(.nchar.cor)
	if(.val> 100){
		.val <- 100
	}
	.nchar.cor.val <- 1 + (.val-8)*0.07
	
}	





.cols <- c(rgb(213,94,0,maxColorValue=255),rgb(0,114,178,maxColorValue=255))

if(gui.input$graphic.type == "pdf"){
pdf(paste(temp.path,"/",colnames(temp.test)[mapping.type],".",temp.name,".pdf",sep = ""),width= pdf.width*2,height = height.pdf)
	#par(mai = c(1,1,0.2,0.2),mfrow= c(1,1))

}else{
postscript(paste(temp.path,"/",colnames(temp.test)[mapping.type],".",temp.name,".eps",sep = ""),width= pdf.width*2,height = height.pdf)
}
par(mai = c(0.9,.nchar.cor.val,0.5,0.1))
print("owiefj")
if(!is.matrix(fisher.results)){fisher.results <- t(as.matrix(fisher.results))}


fisher.results[,2] <- as.numeric(fisher.results[,2])#/d2*100
if(!is.matrix(fisher.results)){fisher.results <- t(as.matrix(fisher.results))}

fisher.results[,4] <- as.numeric(fisher.results[,4])#/d4*100
if(!is.matrix(fisher.results)){fisher.results <- t(as.matrix(fisher.results))}

fisher.results <- fisher.results[order(as.numeric(fisher.results[,2]),decreasing = TRUE),]
if(!is.matrix(fisher.results)){fisher.results <- t(as.matrix(fisher.results))}

fisher.results <- fisher.results[order(as.numeric(fisher.results[,4]),decreasing = TRUE),]
if(!is.matrix(fisher.results)){fisher.results <- t(as.matrix(fisher.results))}

#fisher.results <- fisher.results[order(as.numeric(fisher.results[,4])+as.numeric(fisher.results[,2]),decreasing = TRUE),]
fisher.results.init <- fisher.results[,c(2,4)]
if(!is.matrix(fisher.results.init)){fisher.results.init <- t(as.matrix(fisher.results.init))}
true.vec <- apply(fisher.results.init,1,function(x){all(x!=0)})
fisher.results <- rbind(fisher.results[true.vec,],fisher.results[!true.vec,])


backup <-data.barplot
merge.data.init <- merge.data[,2:3]
if(!is.matrix(merge.data.init)){merge.data.init <- t(as.matrix(merge.data[,2:3]))}
max.v <- max(max(apply(merge.data.init,2,function(x){as.numeric(x)/sum(as.numeric(x))}),na.rm = T)*100,na.rm = T)

if(dim(data.barplot)[1]>1){
	counter <- ceiling(dim(merge.data)[1]/as.numeric(objects.plot))
	start 	<- 1
	
	for(i in 1:counter){
		end <- start+as.numeric(objects.plot)-1
		if(end > dim(backup)[1]){
			end <- dim(backup)[1]
		}
		#stop()
		data.barplot <- fisher.results[,c(2,4)]
		
		data.barplot <- apply(data.barplot,2,as.numeric)
		data2 <- fisher.results[,c(6,7)]
		
		data.barplot <- data.barplot[start:end,]
		data2 <- data2[start:end,]
		if(is.vector(data2)){
			data2 <- t(as.matrix(data2))
			
		}
		.names <- fisher.results[start:end,1]

		start <- start+as.numeric(objects.plot)	

data.barplot[is.na(data.barplot)]<- 0	
	
.max <- max(data.barplot,na.rm = T)/100*10+max(data.barplot,na.rm = T)
.p	<- .max-max(data.barplot)/100*10/2
#if(max(nchar(.names)) > 100){stop("stopped")}
.names[nchar(.names) >100 ]<- paste(substr(.names[nchar(.names)> 100],0,90),"string cutted",sep = " # ")
if(is.vector(data.barplot)){
	data.barplot <- t(as.matrix(data.barplot))
}
data.barplot <- apply((data.barplot),2,rev)
.names <- rev(.names)

if(is.vector(data.barplot)){
	data.barplot <- t(as.matrix(data.barplot[c(2,1)]))	
}else{
data.barplot <- data.barplot[,c(2,1)]
	
}

p.value = p.value.setting
p.v <- rep("",length= length(dim(data2)[1]))
p.v[as.numeric(data2[,1]) < p.value] <- "*"
p.v[as.numeric(data2[,2]) < p.value] <- "**"
p.v[is.na(p.v)] <- ""

.names <- paste(.names,rev(p.v))
print("owiejf")
try.error <- class(try(bp.points <- barplot2(	t(data.barplot),
							beside = TRUE,
							xlim = c(0,max(data.barplot,na.rm = T))
							,las = 2,
							col = .cols,
							xlab = "n", 
							horiz = T,
							cex.names = 0.9,
							names.arg = .names,plot.grid = T))
)
if(try.error == "try-error"){stop()}
#stop()
print("owiejf")

#if(exists("bp.points")){
try(bp.points.add <- apply(bp.points,2,mean))
#text(max(data.barplot,na.rm = TRUE)+0.05*max(data.barplot,na.rm = TRUE),bp.points.add,labels =rev(p.v),xpd = T)

if(dim(data.barplot)[1] < as.numeric(objects.plot)){
	if(dim(data.barplot)[1] > 20){
	pdf.cor <- (dim(data.barplot)[1]-20)*0.1 
}else{
	pdf.cor <- -(dim(data.barplot)[1]-20)*0.1
}

	
}

legend(0-0.01,0,rev(c(test$sample,test$reference)[c(2,1)]),fill = rev(.cols),bg = "#ffffff99",xpd = NA,bty = "n",cex = 1 ,xjust = 1,yjust = 0)
#}


}

}
# stop()
 
graphics.off()
write.csv(fisher.results,paste(temp.path,"/",colnames(temp.test)[mapping.type],".",temp.name,"fisher.results-",".csv",sep = ""))

if(is.vector(merge.data)){merge.data <- t(as.matrix(merge.data))}

.nchar.cor <- nchar(merge.data[,1])
.nchar.cor.val <- 1
merge.data.init <- merge.data[,2:3]
if(!is.matrix(merge.data.init)){merge.data.init <- t(as.matrix(merge.data[,2:3]))}

merge.data[,2:3] <- apply(merge.data.init,2,function(x){as.numeric(x) / sum(as.numeric(x),na.rm = T)})
order.vec <- order(as.numeric(merge.data[,2])/as.numeric(merge.data[,3]))
merge.data <- merge.data[rev(order.vec),]
if(is.vector(merge.data)){merge.data <- t(as.matrix(merge.data))}

fisher.results <- fisher.results[hz.merge.control(fisher.results[,1],merge.data[,1]),]
fisher.results.init <- fisher.results[,c(2,4)]
if(is.vector(fisher.results.init)){fisher.results.init <- t(as.matrix(fisher.results.init))}
all.not.zero <- apply(fisher.results.init,1,function(x){any(as.numeric(x)==0)})
fisher.results 	<- fisher.results[!all.not.zero,]
merge.data		<- merge.data[!all.not.zero,]

#stop()
#merge.data <- merge.data[merge.data[,1] < 0.05,]
merge.data.backup <- merge.data


if(max(.nchar.cor) > 8){
	.nchar.cor.val <- 1 + (max(.nchar.cor)-8)*0.355
	
	
}	

max.data <- 40
if(is.vector(merge.data)){
	merge.data <- t(as.matrix(merge.data))
}
if(dim(merge.data)[1]>1){
	counter <- ceiling(dim(fisher.results)[1]/max.data)
	start 	<- 1

	#stop()
	for(i in 1:counter){
		merge.data <-		merge.data.backup

		end <- start+ max.data -1
		if(end > dim(fisher.results)[1]){
			end <- dim(merge.data)[1]
		}
		data.barplot <- merge.data[,c(2,3)]
		data.barplot <- apply(data.barplot,2,as.numeric)
		data2 <- fisher.results[,c(6,7,1)]
		
		data.barplot <- data.barplot[start:end,]
		data2 <- data2[start:end,]
		.names <- merge.data[start:end,1]
		merge.data <- merge.data[start:end,]

		start <- start+ max.data	
#data2 <- data2[data2[,1] < 0.05,]
#data.barplot <- data.barplot[data2[,1] < 0.05,]

data.barplot[is.na(data.barplot)]<- 0	

data.log <- log2(as.numeric(data.barplot[,1])/as.numeric(data.barplot[,2]))

data.log.vec <- data.log

data.log.vec[!is.infinite(data.log.vec) ] <- ""
#stop()
#data.log[is.infinite(data.log)] <- max(data.log[!is.infinite(data.log)])*0.1
data.log <- data.log[inf.vec <- !is.infinite(data.log)]
data.barplot <- data.barplot[inf.vec,]
data2 <- data2[inf.vec,]
.names <- .names[inf.vec]


.max <- max(data.log[!is.infinite(data.log)],na.rm = T)/100*10+max(data.log[!is.infinite(data.log)],na.rm = T)
.p	<- .max-max(data.log)/100*10/2
.min <- min(data.log,na.rm = TRUE)
if(.min > 0){.min <- 0}
if(.max < 0){.max <- 0}
.col.vec <- c()
.col.vec[data.log > 0] <- .cols[1]
.col.vec[data.log < 0] <- .cols[2]

order.names <- order(.names)
data.log <- data.log[order.names]
data.log.vec <- data.log.vec[order.names]
.col.vec <- .col.vec[order.names]
.names <- .names[order.names]
data2 <- data2[order.names,]

.names <- .names[order.log <- order(data.log)]
.col.vec <- .col.vec[order.log]
data.log.vec <- data.log.vec[order.log]

data.log <- data.log[order.log]
data2 <- data2[order.log,]



if(objects.plot == "full"){
	objects.plot <- length(data.log)
	as.numeric(objects.plot)
	
}
if(is.na(objects.plot)){objects.plot <- 20}

if(as.numeric(objects.plot)> 20){
	height.pdf <- 4 + (objects.plot-20)*0.1 
	pdf.cor <- (objects.plot-20)*0.2
}else{
	height.pdf <- 8
	pdf.cor <- 0
}

.nchar.cor <- nchar(.names)
.nchar.cor.val <- 1


if(max(.nchar.cor) > 8){
	.val <- max(.nchar.cor)
	if(.val> 100){
		.val <- 100
	}
	.nchar.cor.val <- 1 + (.val-8)*0.065
	
}	
if(1==0){

if(gui.input$graphic.type == "pdf"){
pdf(paste(temp.path,"/",colnames(temp.test)[mapping.type],".",temp.name,"log2-barplot-mapping.pdf",sep = ""),width= 7,height = height.pdf)
	#par(mai = c(1,1,0.2,0.2),mfrow= c(1,1))

}else{

postscript(paste(temp.path,"/",colnames(temp.test)[mapping.type],".",temp.name,"barplot-mapping.pdf",sep = ""),width= 7,height = height.pdf)
}
par(mai = c(0.9,.nchar.cor.val,0.6,0.1))

#text(0,bp.points,data.log.vec,col = "grey",pos = 4)

#if(exists("bp.points")){
p.value = p.value.setting
p.v <- rep("",dim(data2)[1])
p.v[as.numeric(data2[,1]) < p.value] <- "*"
p.v[as.numeric(data2[,2]) < p.value] <- "**"
#p.v[as.numeric(data2) > p.value] <- ""
try(bp.points <- barplot(t(data.log),beside = TRUE,xlim = c(.min,.max),las = 2,col = .col.vec,xlab = "", horiz = T,cex.names = 0.9,names.arg = paste(.names,(p.v)),cex.lab = 0.5))
#legend(max(data.log),0,legend = paste("log2( n",test$sample,"/ n",test$reference,")"),xpd = T,xjust = 1)
mtext(paste("log2( n",test$sample,"/ n",test$reference,")"),side = 1,line = 3,adj = 1,xpd = T)
bp.points.add <- bp.points


signi.pos <- c()
signi.pos[data.log >0] <- 4
signi.pos[data.log <0] <- 2
signi.pos[is.na(signi.pos)] <- 1

#text(0,rev(bp.points.add+(bp.points.add[1]-bp.points.add[2])/3*0),labels = rev(p.v),xpd = T,col = "grey",cex = 2,pos =rev(signi.pos))
#legend("topright",c("Sample","Reference"),fill = c("black","white"),bg = "#ffffff99")
#}

legend(min(t(data.log),na.rm = T),max(bp.points)*1.3-pdf.cor,(c(test$sample,test$reference)),fill =(.cols),bg = "#ffffff99",xpd = NA,bty = "n")

#legend(0,max(bp.points)*1.17,rev(c(test$sample,test$reference)),fill = rev(.cols),bg = "#ffffff99",xpd = NA)



}

}
try(graphics.off())

}
#close(pb)
#tkmessageBox(title="Message",message=paste("Wrote export into", getwd()),icon="info",type="ok")
#hz.show.path(temp.path)

t("../")

#stop("end of script")

if(1==0){
####
####
# test 
fisher.results.points <- fisher.results[,c(2,4)]
fisher.results.points <- apply(fisher.results.points,2,function(x){as.numeric(x)/sum(as.numeric(x))})
rownames(fisher.results.points) <- fisher.results[,1]

fisher.results.signi <- fisher.results.points[fisher.results[,7] < 0.05,]


fisher.results.signi <- fisher.results.signi[order(fisher.results.signi[,2],decreasing = TRUE),]


fisher.results.signi <- fisher.results.signi[order(fisher.results.signi[,1],decreasing = TRUE),]






.max <- max(fisher.results.signi)
pdf("signi-plot.pdf")
par(mai = c(1,0.5+(max(nchar(rownames(fisher.results.signi)))-6)*0.08,0.1,0.1))

try(bp.points <- barplot(t(fisher.results.signi),beside = TRUE,las = 1,col = c(1,0),xlab = "n in %", horiz = TRUE))

try(graphics.off())
}

}#if exists merge.data
}#mapping type loop
}# i loop
}
print("finished volcano mapping barplot")
}#function
