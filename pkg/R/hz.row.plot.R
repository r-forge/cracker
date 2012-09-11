hz.row.plot <-
function(
	x = NA,
	path 			= NA,
	name 			= NA,
	y.axis 			= NA,
	plot.type 		= "b",
	plot.col		= 1,
	barpl			= TRUE,
	sample.names 	= NA,
	sub 			= NA,
	h.abline 		= NA,
	v.abline 		= NA,
	x.ylab			= "intensity",
	x.xlab			= "samples",
	ui 				= NULL,
	show.sd			= NULL,# matrix containing the correpondend sd
	sd.rel  		= FALSE,
	.descr			= NULL, # if sd is relativ sd!
	.aov			= NULL,
	.ttest			= NULL,
	p.v				= 0.01,
	sub.cex 		= 0.5,
	inf.m			= NA,
	graphic.type 	= "pdf",
	.design			= NULL,
	time.groups		= T,
	group.barplot 	= FALSE,
	lineplot.beside = F,
	gui.input,prog.max,ratio.prog,pb,hz.exp.des.parse.data2,colorblind.set,.col

){
	if(!is.null(show.sd)){
		show.sd[is.na(show.sd)] <- 0
		
	}
	
	if(!exists("prog.max")){prog.max <- 10000}
sd.po <- 0
	.aov.cor <- p.adjust(.aov,gui.input$p.adjust.method)
	
	p.v <- as.numeric(p.v)
	tempmean <- x
	total = dim(tempmean)[1]
	if(is.null(ui)== FALSE){
	#ui$setProgressBar(pb, 0, label=paste( "0",  "% done"))
	}

	if(is.na(sample.names)) {
		col.x <- colnames(tempmean)
	} else { 
		col.x <- sample.names
	}
	row.x <- rownames(tempmean)
	x <- apply(x,2,as.numeric)
	
	# 
	
	
	if(is.na(path) == FALSE) {
		#wd <- getwd();setwd(path)
	}
	if(is.na(name)) {
		pdf.name <- paste("single-protein-plot-",Sys.Date(),".",graphic.type,sep = "")
	}else{
		pdf.name <- name
		}
		
	if(max(nchar(col.x)) > 40){
		col.x <- substr(test,nchar(test)-40,nchar(test))
		col.x <- paste("...",col.x)
	}	
		
	if(max(nchar(col.x)) > 4){
	oma.val <- 0.1+(max(nchar(col.x)-4))*0.45
	print(oma.val)
	}else{
	oma.val <- 0.1
	}
	if(time.groups== T & barpl == F & group.barplot == F){oma.val<- 1}
	
	if(max(nchar(col.x)) > 4){
	height <- 5+(max(nchar(col.x)-4))*0.1
	}else{
	height <- 5
	}
#print(col.x)
init.width <- 8
if(dim(as.matrix(x))[2] > 30){
	width <- init.width + (dim(as.matrix(x))[2]-30)*0.13
}else{width = init.width}

if(lineplot.beside){
	width <- width+0.8* length(unique(.design$Group))
	
}

if(graphic.type == "pdf"){
	pdf(pdf.name,pointsize = 13,width = width,height = height)
	try(par(oma = c(oma.val,0.1,0.1,0.1),mai = c(1,1.2,0.5,0)))

}else{
	dir.create(.wd.set <- pdf.name)
}






	for(i in 1:dim(tempmean)[1]
	) {
		
##############	GUI
	ratio.prog2 <- (prog.max/8)/total

	pb.check	<- class(try(ui$setProgressBar(pb, i*ratio.prog2, label=paste(round(i/dim(tempmean)[1]*100),  "% protein barplots"))))
	pb.check	<- class(try(ui$setProgressBar(pb, i*ratio.prog2, label=paste(round(i/dim(tempmean)[1]*100),  "% protein barplots"))))

while(pb.check == "try-error1"){
		print("Warning: User closed window!")

		try(pb 			<- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300))
		pb.check	<- class(try(ui$setProgressBar(pb, i*ratio.prog2, label=paste(round(i/dim(tempmean)[1]*100),  "% done"))))
}
##############	
		
		if(graphic.type == "eps"){
		
	postscript(paste(".",.wd.set, paste(row.x[i],".eps",sep = ""),sep = "/"), paper = "special",onefile = FALSE,horizontal = FALSE,pointsize = 13,width = width,height = height)
	try(par(oma = c(oma.val,0.1,0.1,0.1),mai = c(1,1.2,0.5,0)))

		}
		
		print(oma.val)
		
		if(length(dim(tempmean)[1]) > 1000) {
			limit = 1000
		} else {
			limit = 100
		}
		if(i%%limit==0) {
			cat(paste("Went through",i,"proteins out of",dim(tempmean)[1]),"\n")
		}
		
		temp.y <- as.numeric(tempmean[i,])

		if(is.na(y.axis)) {
			temp.x <- c(1:dim(tempmean)[2])
		}
		if(unique(is.na(temp.y)) == TRUE &length(unique(is.na(temp.y))) == 1) {
			temp.y[is.na(temp.y)] <- 0
		}
		
		if(is.na(sub) == FALSE) {
			sub.x = sub[i]
		} else {
			sub.x = NA
		}
		
		range.y 	<- range(temp.y,na.rm = TRUE)
		
		range.y[1] 	<- min(temp.y,na.rm = TRUE)
		if(length(show.sd) != 0){
			if(length(show.sd) != 0 & unique(dim(x) == dim(show.sd)) == 1 & unique(dim(x) == dim(show.sd)) == TRUE){
			library(gplots)
			
			sd.po <- 0
			sd.po <- as.numeric(show.sd[i,])
			if(sd.rel == TRUE){
				if(range(as.numeric(sd.rel))[2] <= 50){
				sd.po <- as.numeric(temp.y) * as.numeric(sd.po)}else{
					sd.po <- as.numeric(temp.y) * as.numeric(sd.po) / 100
					}
				}
				sd.po[is.na(sd.po)]  <- 0
				range.y <- range(as.numeric(temp.y)+as.numeric(sd.po),na.rm = TRUE)
				range.y[1] <- min(temp.y,na.rm = TRUE)	
				#print(range.y)
				#print(temp.y)
				}
				range.y[2] <-  range.y[2]*1.03
#print(sd.po)
#print(temp.y)
		}
		
		if(is.null(.aov) == FALSE){
			#print(p.v)
			#print(i)
			#print(.aov[i])
			if(is.na(as.numeric(.aov[i]))){.aov[i] <- 1; .aov.cor[i] <- 1}
			
			
			if(as.numeric(.aov.cor[i]) < as.numeric(p.v)){	
				row.x[i] <- paste(row.x[i],"*")
			}
			if(as.numeric(.aov[i]) < as.numeric(p.v)){	
				row.x[i] <- paste(row.x[i],"*",sep = "")
			}
			#print(.aov[i])
			#if(i == 2){stop()}
		}
		
		if(is.null(.ttest) == FALSE){
				if(is.na(sample.names)) {
					col.x <- colnames(tempmean)
				} else { 
					col.x <- sample.names
		}
				
				
				p.v <- 0.01
				.ttest.i	<- as.numeric(.ttest[i,])
				.ttest.i.star	<- .ttest.i
				.ttest.i.star[.ttest.i <= p.v] <- "**"
				.ttest.i.star[.ttest.i <= 0.05 & .ttest.i > p.v ] <- " *"
				.ttest.i.star[.ttest.i > 0.05]  <- "  "
				.ttest.i.star[is.na(.ttest.i)]  <- "  "
				
				col.x 		<- paste(col.x,.ttest.i.star)
		
		}
		
		
		#####
 		plot.timeline <- TRUE

		if(!is.null(.design) & time.groups){
				.design.plot <-(.design[,c(5,3,7)])				



				if(!as.logical(gui.input$raw)|all(as.logical(gui.input$raw),gui.input$calc.empai,gui.input$empai.sd)){		
					.design.plot <- unique(.design[,c(2,3,7)])				

					for(test.design in unique(.design.plot$Group)){
					 	if(length(unique(.design.plot$Experiment[.design.plot$Group == test.design])) != length(.design.plot$Time[.design.plot$Group == test.design])){
					 		plot.timeline <- FALSE
					 	}else{plot.timeline <- TRUE}
					}				
				}	
				
				#plot.timeline <- FALSE
				
				temp.x.m <- list()
				col.vec <- list()
				col.legend <- c()
				temp.x.all <- c()
					.design.plot.backup	<- .design.plot
					plot.col.backup		<- plot.col
	
	
				.names <- c()
				for(set.time.groups in unique(.design.plot$Group)){
					#print(set.time.groups)
					if(set.time.groups ==1){
						.design.plot.backup	<- .design.plot
					}else{
						.design.plot		<- .design.plot.backup
					}
					.design.plot.order <- hz.merge.control(.design.plot[,1],gsub(" ","",col.x))		
					
					.design.plot 		<- .design.plot[.design.plot.order[!is.na(.design.plot.order)],]
					plot.col			<- plot.col.backup[.design.plot.order[!is.na(.design.plot.order)]]
					

					 if(!exists("sd.po")){
					 	sd.po <- rep(0, length(temp.y))
					 }
#if(i == 2){stop()}
					 temp.time.x <- cbind(
					 	.design.plot$Time[.design.plot$Group == set.time.groups],
					 	temp.y[.design.plot$Group == set.time.groups],	
					 	sd.po[.design.plot$Group == set.time.groups],
					 	inf.m[i,] [.design.plot$Group == set.time.groups],		
					 	col.x[.design.plot$Group == set.time.groups]				 						)
			 			.names <- c(.names,set.time.groups)	
				 				
					sd.po <- as.numeric(sd.po) 
					
					 
					temp.x.all <- c(temp.x.all,
					 	.design.plot$Time[.design.plot$Group == set.time.groups]					 	)
					 
					 colnames(temp.time.x) <- c("group","intensity","sd","n","name")
					temp.x.m[[set.time.groups]] <- temp.time.x
				
					
					
					col.vec[[set.time.groups]]	<- plot.col[.design.plot$Group == set.time.groups]	
					col.legend <- c(col.legend,plot.col[.design.plot$Group == set.time.groups][1])	
					
				}

				names(temp.x.m) <- .names
								if(!is.numeric(.design.plot$Time)){plot.timeline <- FALSE;group.barplot <- TRUE}

			}
			
		#####
		if(barpl == FALSE){
			#stop()
			### map .design for setting timegroups
					#print(plot.timeline)	
					
			if(plot.timeline & time.groups ){
				plot.data.all <- c()
				
temp.lim.fun <- function(x){temp.lim <- c()
for(k.t in 1:length(names(x))){
	temp.k.t <- x[[k.t]]
	temp.lim <- c(temp.lim, as.numeric(temp.k.t[,2])+as.numeric(temp.k.t[,3]))
}
return(temp.lim)
}
error.try <- class(try(temp.lim <- temp.lim.fun(temp.x.m))				
)
				
for(plot.matrix in hz.merge.control(names(temp.x.m),as.character(unique(.design.plot$Group)))
){
print(plot.matrix)
					plot.data <- temp.x.m[[plot.matrix]][,1:3]
					assign("temp.x.m",temp.x.m,envir = .GlobalEnv)
					print(plot.data)
					if(is.vector(plot.data)){
						plot.data <- t(as.matrix(plot.data))

					}
					#print(plot.data)
					plot.data <- apply(plot.data,2,as.numeric)
					#plot.data[,2] <- plot.data[,2]-min(plot.data[,2],na.rm = TRUE) 
					plot.data.all <- c((as.numeric(plot.data[,2])-plot.data[,3]),(as.numeric(plot.data[,2])+as.numeric(plot.data[,3]))*1.05)
					}
					
										assign("plot.data.all",plot.data.all,envir = .GlobalEnv)

					
			if(lineplot.beside){
				layout(matrix(c(rep(1,length(unique(.design$Group))+1),2:(length(unique(.design$Group))+2)),2,length(unique(.design$Group))+1, byrow =T), heights = c(0.2,1),widths = c(0.4,rep(1,length(unique(.design$Group)))))
main.temp <- ""
				par(mar = c(0,0,0,0),mai = c(0,0,0,0),oma = c(0,0,0,0))
				plot.new()
				legend("topleft",legend = c(row.x[i],sub.x),box.col = "transparent",cex = 1.5, xjust = 1)
				
								par(mai = c(0.1,0.7,0.1,0))
plot(
					0,
					0,
					type = "n",
					main = main.temp,
					col = plot.col,
					xlab = "",
					ylab = x.ylab,
					#names.arg = col.x,
					ylim = range(plot.data.all,na.rm = T),
					xlim = 	range(.design.plot$Time,na.rm = T) ,
					frame = FALSE,
					axes = F,
					lwd = 3,
					cex = 1.1,
					mgp = c(0,0,0)
					
				)
			}else{

				main.temp <- row.x[i]
				try(par(oma = c(oma.val,0.1,0.1,0.1),mai = c(1,1,0.5,0)))

				plot(
					0,
					0,
					type = "n",
					main = main.temp,
					col = plot.col,
					xlab = x.xlab,
					ylab = x.ylab,
					#names.arg = col.x,
					ylim = range(plot.data.all,na.rm = T),
					xlim = 	range(.design.plot$Time,na.rm = T) ,
					frame = FALSE,
					lwd = 3,
					cex = 1.1
					
				)
				grid()
				}	
								print("test")

				
							#	par(mai = c(0.6,0,0.1,0))

				
				
				

				temp.temp.y <- c()
				temp.temp.x <- c()
				temp.temp.sd.po <- c()
				for(plot.matrix in hz.merge.control(names(temp.x.m),as.character(unique(.design.plot$Group)))){
					

				if(lineplot.beside ){
					
						print("tr")
	#try(par(oma = c(oma.val,0.1,0.1,0.1),mai = c(1,0,0.5,0)))
				par(mai = c(0.6,0.3,0.1,0))
				plot(
					0,
					0,
					type = "n",
					#main = row.x[i],
					col = plot.col,
					xlab = x.xlab,
					ylab = x.ylab,
					#names.arg = col.x,
					ylim = range(plot.data.all,na.rm = T)*c(1,1.2),
					xlim = 	range(.design.plot$Time,na.rm = T) ,
					frame = FALSE,
					lwd = 3,
					cex = 1.1
					
					
				)						
				grid()

					}
					
					if(any(duplicated(temp.x.m[[plot.matrix]][,1]))& i == 1){
						tkmessageBox(message = "Duplicated values in time column in experimental design file. Settings might not be optimal for visual output!")
						
					}
					plot.data <- temp.x.m[[plot.matrix]][,1:3]
					#.names <- names(temp.x.m[plot.matrix])
					
					if(is.vector(plot.data)){
						plot.data <- t(as.matrix(plot.data))

					}
					#print(plot.data)
					plot.data <- apply(plot.data,2,as.numeric)
					#plot.data[,2] <- plot.data[,2]-min(plot.data[,2],na.rm = TRUE) + 0.1
					
					
					print(plot.data[,1:2])
					points(plot.data[,1:2]	,
							type = "b",
							col = col.vec[[plot.matrix]],
							lwd = 3,
							cex = 1.1
							
					)
					#print( col.vec[[plot.matrix]])
					#graphics.off()
					#stop()
					if(lineplot.beside){
			
					plotCI(plot.data[,1],as.numeric(plot.data[,2]),ui =as.numeric(plot.data[,2])+as.numeric(plot.data[,3]),li =as.numeric(plot.data[,2])-as.numeric(plot.data[,3])
						,type = "p",add = TRUE,col =col.vec[[plot.matrix]], gap = 0,lwd = 3)
			
			
			legend(	"top", 
					.names[plot.matrix],
					fill = col.vec[[plot.matrix]],
					border = "transparent", 
					bg = "#FFFFFF99",cex =1.2,box.col = "transparent")
					}
					temp.plot.data <- cbind(plot.data,col.vec[[plot.matrix]])
					temp.temp.y <- rbind(temp.temp.y, temp.plot.data)
				}
				plot.col.backup <- plot.col
				temp.x 	<- as.numeric(temp.temp.y[,1])
				temp.y 	<- as.numeric(temp.temp.y[,2])
				sd.po 	<- as.numeric(temp.temp.y[,3])
				plot.col <- temp.temp.y[,4]

			
			temp.legend.input<- unique(cbind(as.character(.design.plot[,2]), plot.col.backup))
			
			
			#graphics.off()
		#		stop()

			}else{
		
			plot.col<- hz.exp.des.parse.data2[hz.merge.control(hz.exp.des.parse.data2[,2],col.x),1]
			
			plot(	temp.x,
					temp.y,
					main = row.x[i],
					col = plot.col,
					xlab = "",
					ylab = x.ylab,
					#names.arg = col.x,
					ylim = range(c(as.numeric(sd.po)+(as.numeric(temp.y)),as.numeric(temp.y) -as.numeric(sd.po)),na.rm = T),
					axes = FALSE,
					lwd = 3,
					cex = 1.1,
					type = "n"
			)
			
			grid()
			points(temp.x,temp.y,col = "grey",type = "b",lwd = 3,cex = 1.1)

			points(	temp.x,
					temp.y,
					#type = plot.type,
					main = row.x[i],
					col = plot.col,
					xlab = "",
					ylab = x.ylab,
					#names.arg = col.x,
					ylim = range(c(sd.po+temp.y, temp.y -sd.po),na.rm = T),
					axes = FALSE,
					lwd = 3,
					cex = 1.1
			)
			
			
					axis(2)
		axis(1,temp.x,col.x,las = 2)

			}
			
			
			
	
			
			
			
			
		} else {
			
			if(length(show.sd) != 0){
				if(length(show.sd) != 0 & unique(dim(x) == dim(show.sd)) == 1 & unique(dim(x) == dim(show.sd)) == TRUE){
					plot.ci <- TRUE
				}else{
					plot.ci <- FALSE
					sd.po <- rep(0,length(temp.y))
				}
			}else{
				plot.ci <- FALSE
				sd.po <- rep(0,length(temp.y))
			}


			if(time.groups& any(plot.timeline,group.barplot)){
				for(time.groups.i in 1: length( names(temp.x.m))){
					if(time.groups.i == 1){
						time.groups.temp <- as.numeric(temp.x.m[[time.groups.i]][,c(2)] )
						time.groups.n	 <- temp.x.m[[time.groups.i]][,c(4)] 
						time.groups.sd	 <- as.numeric(temp.x.m[[time.groups.i]][,c(3)] )
						time.groups.names<- temp.x.m[[time.groups.i]][,c(5)] 
						
						#print((temp.x.m[[time.groups.i]][,4]))


					}else{
				#		 				temp.x.m.order <- hz.merge.control(temp.x.m		[[time.		groups.i]][,5],time.groups.temp[,5])
						temp.x.m.temp.time  <- temp.x.m[[time.groups.i]][,2]
														 
						time.groups.temp <- cbind(time.groups.temp,temp.x.m.temp.time)		

		 				
						temp.x.m.vec <- temp.x.m[[time.groups.i]][,4]
												
		 				time.groups.n <- cbind(time.groups.n, temp.x.m.vec)



						#print((temp.x.m[[time.groups.i]][,4]))
		 				time.groups.sd <- cbind(time.groups.sd,as.numeric(temp.x.m[[time.groups.i]][,3]))
		 				#print(temp.x.m[[time.groups.i]][,5])
temp.x.m.temp.name  <- temp.x.m[[time.groups.i]][,5]
						if(is.vector(time.groups.names)){
							time.groups.names <- as.matrix(time.groups.names)
						}
						substraction.test	<- dim(time.groups.names )[1]-length(temp.x.m.temp.name)						
						
						if(substraction.test > 0){
							temp.x.m.temp.name 			<- c(temp.x.m.temp.name,rep(NA, substraction.test))
							names(temp.x.m.temp.name)[(length(temp.x.m.temp.name)-substraction.test+1):length(temp.x.m.temp.name)] 	<- 	NA
						}
						if(substraction.test < 0){
							temp.x.m.temp.name 			<- rbind(time.groups.names,matrix(NA, (substraction.test),dim(time.groups.names)[2]))
							names(temp.x.m.temp.name)[(length(temp.x.m.temp.name)-substraction.test+1):length(temp.x.m.temp.name)] 	<- 	NA
						}
								 
						time.groups.names <- cbind(time.groups.names, temp.x.m.temp.name)					
					
					}
				}
				
				colnames(time.groups.temp) <- names(temp.x.m)
				rownames(time.groups.temp) <- NULL
				time.groups.temp[is.na(time.groups.names)] <- NA
					
				colnames(time.groups.sd) <- names(temp.x.m)
				rownames(time.groups.sd) <- NULL
				time.groups.sd[is.na(time.groups.names)] <- NA

				colnames(time.groups.n) <- names(temp.x.m)
				rownames(time.groups.n) <- NULL
				time.groups.n[is.na(time.groups.names)] <- NA

				colnames(time.groups.n) <- names(temp.x.m)
				rownames(time.groups.names) <- NULL
				
				
				

								temp.y 	<- (time.groups.temp)
				temp.y <- apply(temp.y,2,as.numeric)
				sd.po 	<- (time.groups.sd) 
				col.x2 <- as.vector((time.groups.names))
				#hz.merge.control(colnames(x), col.x2)
				#stop()
				
				if(barpl&time.groups){
					temp.y <- t(temp.y)
					col.x2 <- t(time.groups.names)
					sd.po <- t(sd.po)
					color <- col.x2
					

				for(col.grep in 1:length(colnames(tempmean))){
					temp.col.grep <- grep(colnames(tempmean)[col.grep],color,fixed = T)
					#print(temp.col.grep)
					color[temp.col.grep] <- .col[col.grep]
					
				}
				#print(color)
				}
				
							plot.col <- color

				
				#stop()
		try(				inf.m.temp <- as.vector(time.groups.n)[!is.na(as.vector(time.groups.n))]
		)
		
		if(length(inf.m.temp)==0){inf.m.temp <- rep(0,dim(inf.m)[2])}
				inf.m[i,] <- inf.m.temp			
			}
if(barpl&!time.groups){
					col.x2 <- col.x
				}
					
col.temp <- hz.merge.control(hz.exp.des.parse.data2[,2],col.x2)
plot.col <- hz.exp.des.parse.data2[col.temp,1]

library(gplots)
			try.error <- class(try(temp.min <- min(temp.y -sd.po,na.rm = T)*0.9))
			if(try.error == "try-error"){
			try.error <- class(try(temp.min <- 0))

			}
			
			  par(lwd = 2)

				test <- barplot2(
				temp.y-temp.min,
				main = row.x[i],
				col = plot.col,
				names.arg = col.x2,
				las = 2,
				ylab = x.ylab,
				plot.ci = TRUE,
				ci.l = temp.y+sd.po,
				ci.u = temp.y-sd.po,
				ci.color = "black",
				ci.lwd = 2,
				lwd = 2,
				beside = TRUE,
				xpd=F, 
				yaxp=c(0,max(as.numeric(tempmean),na.rm = TRUE), 4) 
				,
				ylim = c(temp.min,max(as.numeric(sd.po)+as.numeric(temp.y),na.rm = TRUE)), 
				offset = temp.min,
				mgp = c(3.9,1,0),
				plot.grid = T,
				grid.col = "darkgrey"
			)
			temp.x <- as.numeric(test)


		}
										 			#	stop()

		if(length(.descr) != 0){
		mtext(.descr[i],3,adj = 0)
		}
	
		if(is.matrix(inf.m)){
			#print("tes")
			l.pos <- rep(3,length(temp.x))
			l.pos[is.infinite(as.numeric(inf.m[i,])) & as.numeric(inf.m[i,] ) < 0 ] <- 1
			#print(l.pos)
			#print(temp.x)
			
			if(barpl&time.groups){n.input <- as.vector(t(time.groups.n))}else{n.input <- inf.m[i,]}
			if(exists("temp.min")){			text(temp.x,y = temp.min,labels=n.input,col = "white",pos = l.pos,cex = 0.8)
}else{
	
}
			
		}
		if(length(show.sd) != 0& barpl == FALSE){
			if(length(show.sd) != 0 & unique(dim(x) == dim(show.sd)) == 1 & unique(dim(x) == dim(show.sd)) == TRUE){
			
			
			#plotCI(temp.x,as.numeric(temp.y),ui =temp.y+sd.po,li = temp.y-sd.po,type = "n",add = TRUE,col = "darkgrey", gap = 0.3,lwd = 2.5)
			if(!lineplot.beside){
						plotCI(temp.x,as.numeric(temp.y),ui =as.numeric(temp.y)+as.numeric(sd.po),li =as.numeric(temp.y)-as.numeric(sd.po)
						,type = "p",add = TRUE,col = plot.col, gap = 0,lwd = 3)
			
			
			if(gui.input$time.grouped){
					legend("topright", temp.legend.input[,1],fill = temp.legend.input[,2],border = "transparent", bg = "#FFFFFF99",cex = 0.8,xpd = T)
			#stop()
			}
			
			
			}
			}}
			
		
		
		if(is.na(sub) == FALSE&!lineplot.beside) {
			
			
			mtext(paste(substring(sub.x,0,125),"..."),adj = 0,cex = sub.cex)
		}
		if(is.null(ui)== FALSE){
			
			

		
		
		}

		if(is.na(h.abline[1]) == FALSE) {
			abline(h= h.abline)
		}
		if(is.na(v.abline[1]) == FALSE) {
			abline(v= v.abline)
		}

	if(graphic.type == "eps"){
		graphics.off()	
	}
	

#	dev.off()
#	stop("er")
				
	}
	

	graphics.off()
	print(paste("Printed row.plot in ",getwd()))

}
