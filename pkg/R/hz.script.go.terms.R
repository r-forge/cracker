hz.script.go.terms <-
function(.data2, kmeans.cluster.output, info.add,gui.input,kmeans.col,kmeans.at,kmeans.list,prog.max,pb,ui,plot.type,.col,colorblind.set,color.blind, backup.go.input.agg){
		if(!exists("ratio.prog")){ratio.prog <- 1000}

	########## GUI 
	label.string	<- paste( "Loading mapping table")
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4.5, label= label.string)))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4, label= label.string)))
		}
		
	##############
	
	
#if(!exists(".go")){
print("Start Loading mapping library")
#.go <- read.delim(paste(path1,"mapping-library/",GO.library,sep = ""), stringsAsFactors = FALSE,header = FALSE)
temp.list.mapping <- list.files(paste(path.package("cRacker"),"/data/",sep = ""))

temp.grep.mapping <- grep(gui.input$go.library,temp.list.mapping)

if(length(temp.grep.mapping)==0){
try.error <- class(try(load(paste(path.package("cRacker"),"/data/cRackerMapping-",gui.input$go.library,sep = ""))
))
}else{

	
try.error <- class(try(load(paste(gui.input$cracker,"/",gui.input$go.library,sep = ""))))
}



print("Finished Loading mapping library")

########## GUI 
	label.string	<- paste( "Preparing mapping")
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4.5, label= label.string)))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4, label= label.string)))
		}
		
##############

#}
#.data.cluster <- read.csv("/Users/cRacker-DEMO/henrik/sterol-stuff/DRM/20100925-DRM-Mg/ms-analysis--raw-csv-2011-02-01/signi-kmeans-cluster.csv", stringsAsFactors = FALSE)

.go[,1] <- tolower(.go[,1])

.data.cluster <- as.matrix(kmeans.cluster.output)
.data.cluster <- cbind(rownames(.data.cluster),.data.cluster)
.data.cluster[,1] <- gsub(" $","",.data.cluster[,1])

# grep data proteins in .go.cur 

cluster.split <- strsplit(.data.cluster[,1],".",fixed = TRUE)


iso.test  <- grep(".",.go[,1],fixed = TRUE)
if(length(iso.test) == 0){
isoforms = FALSE}else{isoforms = TRUE}

isoforms = FALSE
if(isoforms == FALSE){
	for(ab in 1: dim(.data.cluster)[1]){
		temp.ab  <- cluster.split[[ab]][1]
		.data.cluster[ab,1] <- tolower(temp.ab)	
	}
}else{
	
.data.cluster[,1] <- tolower(.data.cluster[,1])
}


########## GUI 
	label.string	<- paste( "Starting mapping")
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4.5, label= label.string)))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*4, label= label.string)))
		}
##############

grep.prot <- hz.merge.control(tolower(.go[,1]),tolower(.data.cluster[,1]))
grep.prot <- grep.prot[!is.na(grep.prot)]




if(length(grep.prot) != 0){
	
compare.to.all <- TRUE
	
.go.data.cluster	<- .go[grep.prot,]



a.input <- unique(.go[,3])
a.p.value <- c()
a.p.value.test <- c()
if(plot.type == 1){backup.go.input.agg 	<- list()}

for(a in 1:length(a.input)){
#		print(a)

print(paste("Started mapping search for",a.input[a]))	
	
go.input 		<- unique(.go.data.cluster[.go.data.cluster[,3] == a.input[a], ]) # type 

if(dim(go.input)[1] !=0){
	
	go.input.agg <- aggregate(go.input[,2],list(go.input[,2]),length)

	if(plot.type == 1){
		backup.go.input.agg[[a]]		<- go.input.agg		
	}
	
	
		# go through clusters
	nchar.input <- max(nchar(go.input.agg[,1]))
	if(nchar.input > 30){nchar.input <- (nchar.input-30)*0.08 
	}else{nchar.input <- 0}	
	
	if(gui.input$graphic.type == "pdf"){
		pdf(paste("mappings.",a.input[a],".pdf",sep = ""),width = 7 + nchar.input,height = 9)
	}
	.wd <- getwd()
	
	if(dim(go.input.agg)[1] > 25){
		dotchart.cex.subtract <- 0.01*(dim(go.input.agg)[1]-25)
	}else{
		dotchart.cex.subtract <- 0
	}
	dotchart.cex <- 1 - dotchart.cex.subtract
	if(dotchart.cex < 0.1){dotchart.cex = 0.1}

	go.input.agg <- go.input.agg[order(as.numeric(go.input.agg[,2])),]

	#dotchart(go.input.agg[,2],paste(go.input.agg[,1]), main = paste("All; ",length(unique(go.input[,1]))," Proteins",sep = ""),cex = dotchart.cex)
	
	for(b in 1:length((unique(.data.cluster[,2])))){
		
		
		
		########## GUI 
	label.string	<- paste( "Map Cluster",b,"type",a.input[a])
	
	
	.init 	<- (ratio.prog*0.5)
	.init2 <- length((unique(.data.cluster[,2])))*length(a.input)
	pr.input		<- ratio.prog*4.5 + .init/.init2*(b + (a-1)*length((unique(.data.cluster[,2]))))
	print(pr.input)
	pb.check <- class(try(ui$setProgressBar(pb, pr.input, label= label.string)))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, pr.input, label= label.string)))
		}
		
	##############		
		print(paste("Going through cluster",b))	 
		temp.b 		<-	sort(unique(.data.cluster[,2]))[b]
		temp.b.grep	<- 	hz.merge.control(tolower(go.input[,1]),as.character(tolower(.data.cluster[.data.cluster[,2]== temp.b,1])))
		temp.b.grep <- temp.b.grep[!is.na(temp.b.grep)]
		
		
		temp.b.setdiff	<- 	setdiff(.data.cluster[.data.cluster[,2]== temp.b,1],go.input[,1])

		
		if(length(temp.b.grep)!= 0){
			temp.b.go 	<- unique(go.input[temp.b.grep,])
			if(is.vector(temp.b.go)){
			temp.b.go <- t(as.matrix(temp.b.go))
			}
			temp.b.agg 	<- aggregate(temp.b.go[,2],list(temp.b.go[,2]),length)
			
			temp.b.prot	<- aggregate(temp.b.go[,1],list(temp.b.go[,2]),function(x){paste(x,collapse = "|")})
			
			temp.b.mapp <- aggregate(temp.b.go[,2],list(temp.b.go[,1]),function(x){paste(x,collapse = "|")})
			temp.b.mapp <- cbind(temp.b.mapp, a.input[a])
			print(temp.b.mapp[1,])

			temp.b.agg <- temp.b.agg[order(temp.b.agg[,2]),]
			
			plot.go <- TRUE
			
		}else{
			plot.go <- FALSE
		}
		
		if(dim(go.input.agg)[1] > 25){
			dotchart.cex.subtract 		<- 0.01*(dim(temp.b.agg)[1]-25)
		}else{dotchart.cex.subtract <- 0}
			dotchart.cex <- 1 - dotchart.cex.subtract
			
			if(dotchart.cex < 0.1){dotchart.cex = 0.1}
	
		#go through GO entries
		if(plot.go){
			b.p.value <- c()
			for(d in 1:dim(temp.b.agg)[1]){
				d1		<- temp.b.agg[d,]
				d2		<- sum(as.numeric(temp.b.agg[,2]))	
				
				if(compare.to.all){
					#if(plot.type == 2){
				#		control.input <- go.input.agg
				#	}else{
						control.input <- backup.go.input.agg[[a]]
					#}
				}else{control.input <- go.input.agg}
				
				d3		<- control.input[control.input[,1]==d1[,1],]
				d4		<- sum(as.numeric(control.input[,2]))
				
				
				d.ratio <- (d1[,2]/d2)/(d3[,2]/d4)
				
				d.con.table 	<- rbind(	c(d1[,2],abs(d1[,2]-d2)),c(d3[,2],abs(d3[,2]-d4)))
				if(all(dim(d.con.table) != c(2,2))){d.con.table <- matrix(1,nrow = 2,ncol = 2)}
				
				d.fisher		<- try(fisher.test(d.con.table)$p.value)
				b.p.value 		<- rbind(b.p.value,c(as.character(d1[1,]),as.numeric(d.fisher),as.numeric(d.ratio)))
				
			}
			#if(a.input[a] == "NA"){stop()}
			
			b.p.value  <- apply(b.p.value,2,as.character)
			if(is.vector(b.p.value)){b.p.value <- t(as.matrix(b.p.value))}
			temp.b.prot 	<- apply(temp.b.prot,2,as.character)
						if(is.vector(temp.b.prot)){temp.b.prot <- t(as.matrix(temp.b.prot))}

			
			
			data.merge <- merge(as.matrix(b.p.value), temp.b.prot,by = 1,all = TRUE)
	
			b.p.value <- cbind(data.merge[,1:4],b,a.input[a],data.merge[5])
			a.p.value <- rbind(a.p.value, b.p.value)
	
			a.p.value.test <- rbind(a.p.value.test, cbind(temp.b.mapp,b))
			
			
			colnames(b.p.value) <- c("mapping","n","p","ratio","cluster","GO-type","proteins")
			
			b.p.value <- b.p.value[order(b.p.value[,1]),]
			b.p.value <- b.p.value[order(as.numeric(as.character(b.p.value[,4]))),]
			
			
			#define vectors for dotchart
			temp.p.val <-  as.numeric(as.character(unlist(b.p.value[,3])))
			p.val.string	<- rep("", length(temp.p.val))
			p.val.color		<- rep(1, length(temp.p.val))
			p.val.color[temp.p.val < 0.05] <- 2
			p.val.color[temp.p.val < gui.input$p.value] <- 2
			
			p.val.string[temp.p.val < 0.05] <- "*"
			p.val.string[temp.p.val < gui.input$p.value] <- "**"
	
		

			
			temp.p.rat				<- unlist(b.p.value[,4])
			temp.p.rat				<- log2(as.numeric(as.character(temp.p.rat)))
			ratio.data.cluster		<- rep("",length(temp.p.rat))
			ratio.data.cluster[temp.p.rat >0.5] <- "+"
		
			ratio.data.cluster[temp.p.rat < -0.5] <- "-"



			if(dim(temp.b.agg)[1] > 25){
				dotchart.cex.subtract <- 0.01*(dim(temp.b.agg)[1]-25)
			}else{
				dotchart.cex.subtract <- 0
			}
			
			if(length(temp.b.setdiff) > 0){
				temp.b.agg <- rbind(c("unknown",length(temp.b.setdiff)), temp.b.agg)
			}


			dotchart.cex <- 1 - dotchart.cex.subtract
			if(dotchart.cex < 0.1){
				dotchart.cex = 0.1
			}

			if(is.data.frame(temp.b.agg)){
		
			if(gui.input$graphic.type == "eps"){
				dir.create(.wd.set <- "mappings-eps")
				
				postscript(paste(.wd,.wd.set,paste("mappings.",b,a.input[a],".eps",sep = ""),sep = "/"),width = 7 + nchar.input,height = 9, paper = "special",onefile = FALSE,horizontal = FALSE)
			}
		
		
				layout(matrix(c(1,2),ncol = 1,nrow=2),heights = c(2,1))
	
				#par(oma = c(0,0,0,0))
				par(mai = c(0.6,0,0.5,0.2),new = F)
				#stop()
				
				
			dotchart(	as.numeric(log2(as.numeric(as.character(b.p.value[,4])))),	
							paste(	b.p.value[,1],
									unlist(b.p.value[,2]),													ratio.data.cluster ,p.val.string), 
									color = p.val.color, 
									main = paste("Cluster ",b,"; ",length(unique(.data.cluster	[.data.cluster[,2]== temp.b,1]))," Proteins",sep = ""),
									cex = dotchart.cex,
									xlab = paste("log2( n in cluster in % ) - log2( n in all clusters in % )"),
									mgp = c(2,1,0),
									cex.lab = 0.8,
									lwd = 2	
									
									)
									
				par(mai = c(1,0.5,0.2,0.2),new = F)
				
				if(is.list(kmeans.col)& length( kmeans.col) >0){
					box.col <- kmeans.col[[b]]
				}else{box.col <- .col}
				
				plot.error <- class(try(plot(1,1,mgp = c(0,0.5,0),type = "n",labels = TRUE,xlab = "",ylab = "",axes = FALSE,main = paste("Median curve of cluster",b),cex.main = 0.6,col = box.col,xlim = c(min(kmeans.at[[b]])*0.7,max(kmeans.at[[b]]))*1.05,add = TRUE,ylim = range(kmeans.list[[b]],na.rm = TRUE))
))				
			if(plot.error!= "try-error"){
				grid()
				try(boxplot(kmeans.list[[b]],at = kmeans.at[[b]],xlab = "",ylab = "",axes = FALSE,cex.main = 0.6,col = box.col,add = TRUE)
)
	try(				axis(2,cex.axis = 0.4))		
	}else{print("failed in plotting kmeans.at");plot(1,type = "n",axis = F)}
#				axis(1,at = kmeans.at[[b]],labels = colnames(kmeans.list[[b]]),las = 2,cex.axis = 0.3,col = box.col)
		
#				temp.mean <- apply(as.matrix(kmeans.list[[b]]),2,function(x){x<-as.numeric(x);x<- median(x,na.rm = TRUE)})
#				points(kmeans.at[[b]],temp.mean,type = "l",col = "red",lwd = 3)

				if(length(box.col) == 1){box.col <- 1}
				for(box.col.i in unique(box.col)){
					
					nchar.template <- 8
					temp.label <- colnames(kmeans.list[[b]])
					temp.label[nchar(temp.label) > nchar.template] <- paste(substr(temp.label[nchar(temp.label) > nchar.template],0,nchar.template),"...",sep = "")

					print(box.col.i)
					box.col.temp.i <- box.col == box.col.i
try(					axis(1,at = kmeans.at[[b]][box.col.temp.i],labels = temp.label[box.col.temp.i],las = 2,col.axis = box.col.i,cex = 0.1,lwd = 2)
)		
					temp.mean <- apply(as.matrix(kmeans.list[[b]]),2,function(x){x<-as.numeric(x);x<- median(x,na.rm = TRUE)})


					if(gui.input$time.grouped){
					points(kmeans.at[[b]][box.col.temp.i],temp.mean[box.col.temp.i],type = "l",col = "grey",lwd = 4)		
					points(kmeans.at[[b]][box.col.temp.i],temp.mean[box.col.temp.i],type = "l",col = box.col.i,lwd = 2)
					}
		
				}
				
				if(!gui.input$time.grouped){
					points(kmeans.at[[b]],temp.mean,type = "l",col = "grey",lwd = 4)		
					points(kmeans.at[[b]],temp.mean,type = "l",col = 2,lwd = 2)
					}
			}else{
				print(paste("Go object is not a dataframe."))
			}
		}
	}
graphics.off()
}
}

a.p.value <- cbind(a.p.value[,c(1,2,3)],p.adjust(as.numeric(a.p.value[,3]),method = gui.input$p.adjust),log2(as.numeric(as.character(a.p.value[,4]))),a.p.value[,c(4,5,6,7)])
colnames(a.p.value) <- c("GO","n","p .unadjusted",paste("p.adjusted",gui.input$p.adjust,sep ="."),"occurence.ratio.log2","occurence.ratio","cluster","type","proteins")


write.csv(a.p.value,"mapped-data.csv")
colnames(a.p.value.test) <- c("accession","mapping","type","cluster")
a.p.value.test <- a.p.value.test[order(a.p.value.test[,1]),]
info.add <- c()

multi.mapping.wrapper <- function(temp.mappings,bin.classes){

if(is.vector(temp.mappings)){temp.mappings <- t(as.matrix(temp.mappings))}

type.i.matrix <- c()
for(type.i in bin.classes){
	temp.type.i<- temp.mappings[temp.mappings$type == type.i,]
	

	if(dim(temp.type.i)[1] > 0){
		type.i.matrix <- c(type.i.matrix,paste(temp.type.i[,2],collapse = "|"))
		}else{
		type.i.matrix <- c(type.i.matrix,"")
	}	
}
names(type.i.matrix) <- unique(a.p.value.test[,3])
return(type.i.matrix)
}

multi.map.matrix <- c()
for(i in 1:length(unique(a.p.value.test[,1]))){
	
	
		########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*5+(ratio.prog/dim(a.p.value.test)[1]*i), label=paste( "Preparing mapping table"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*5+(ratio.prog/dim(a.p.value.test)[1]*i), label=paste( "Preparing mapping table"))))
		}
	##############
	
	
	
	temp.i <- unique(a.p.value.test[,1])[i]
	temp.mappings  <- a.p.value.test[a.p.value.test[,1] == temp.i,]
	multi.map <-	multi.mapping.wrapper(temp.mappings, unique(a.p.value.test[,3]))
			
	grep.i <- grep(temp.i, .data2$proteinlist.info[,2])
	temp.info <- .data2$proteinlist.info[grep.i,]	

	if( is.vector(temp.info)){temp.info <- t(as.matrix( temp.info))}
	if(is.null(temp.info)){temp.info <- matrix()}
	if(dim(temp.info)[1] == 0){
		#stop()
		temp.info	 <- t(as.matrix(rep("no data",dim(temp.info)[2])))
			colnames(temp.info) <- colnames(.data2$proteinlist.info)
	}
	if( is.vector(temp.info)){temp.info <- t(as.matrix( temp.info))}

	if(dim(temp.info)[1] >1){
	temp.info <- temp.info[1,]
	}
	if( is.vector(temp.info)){temp.info <- t(as.matrix( temp.info))}

	if(dim(temp.info)[1] == 0){
		
	temp.info <- rep("NA",dim(.data2$proteinlist.info)[2])	
		
	}
	if( is.vector(temp.info)){temp.info <- t(as.matrix( temp.info))}

	 try(info.add <- rbind(info.add,temp.info))
	 try(multi.map.matrix <- rbind(multi.map.matrix ,c(temp.i,multi.map)))


}


if(dim(multi.map.matrix)[1]== dim(info.add)[1]){
	extended.info <- cbind(multi.map.matrix,info.add)
}else{
	extended.info <- multi.map.matrix
}


write.csv(extended.info,"protein.mapping.csv")

#write.csv(a.p.value.test,"protein.mapping.csv")
try(order.kmeans <- (hz.merge.control(rownames(kmeans.cluster.output),extended.info[,4])))
try(extended.info <- cbind(extended.info,kmeans.cluster.output[order.kmeans,])
)
colnames(extended.info)[dim(extended.info)[2]] <- "kmeans.cluster"
write.csv(extended.info,"protein.mapping.csv")



}else{
print("Search for GO Terms failed.\nCheck if your Accessions are compatible with loaded GO library!")
}
print("Done")
return(list(extended.info = extended.info,backup.go.input.agg=backup.go.input.agg))
}
