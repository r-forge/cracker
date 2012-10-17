hz.matrix.creator <-
function(	
	x, 						# Output from database
	Raw 	= FALSE, 		# treats different raw-files as different samples! Makes only sense with IonIntensities of non merged empais! 
	type 	= "ion", 		# or ion
	merge.method = "mean", 	# can be median or mean
	outlier = "all.below", 	# can be NA, row, all.below
	score 	= 20,			# exclusion treshold for peptides
	norm.tog.pep = FALSE,   # if normalization of raw files of one resultfile is based on common peptides
	row.norm = TRUE,		# if normalisation over row of peptides over conditions....
	shape = 0.5,			# shape vector?
	add.data = FALSE, 		# extract additional data per proteins, takes the best peptide
	ui = NULL,				# the class to use for user interaction (such as progressBars, messageBoxes)
	re.pep.ex= FALSE,
	maxq.exclu = FALSE,
	cbn.prot = NULL,
	ratios = FALSE,
	phospho = FALSE,
	script.shape = FALSE,
	action.dupli = "exclude", # options: c("mean","sum","max","min","exclude")
	
	prot.norm.shape = 1,  # shape vector of reference protein
	exclude.raw = "blank",
	use.raw = FALSE,
	N15 = FALSE,
	row.target.norm = TRUE,
	path.design = "",
	zero.treat 	= NA,
	n15.correct.method 	= "median",
	n15.correct.expect	= 1,
	#n15.correct			= FALSE,
	n15.log2			= TRUE,
	build.matrix		= TRUE, #loads old matrix from binary
	
	calc.empai		= "FALSE",
	length.pep		= 6,
	empai.sd		= FALSE,
	empai.from.msms = TRUE,
	empai.norm 		= "sum",
	
	n.correction	= FALSE,
	group.norm 		= FALSE,
	group.v			= NULL ,
	norm.method		= "mean",
	group.shape		= NA,
	group.filter	= FALSE,
	
	sum.of.total	= "sum",
	graphic.type 	= "pdf",
	gui.input,
	prog.max,pb,
	.length.matrix
){
	
# Init progress bar if not done yet
print("start matrix.creator function")
library(grid)
if (!exists("prog.max")) {
	prog.max <- 10000;
}

if (!exists("pb")) {
	pb = ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300);
}

##############	GUI
	.label <- "Preparing data"
	pb.check	<- class(try(ui$setProgressBar(pb,0, label=.label)))

	while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))
	}
##############




	cut.path <-	function(x){
	x.uni <- unique(x)
	
			.col		<- strsplit(as.character(x.uni),"\\",fixed = TRUE)
			.cols	<- c()

			ratio.prog <- prog.max/length(.col)
			for(i in 1: length(.col)) {
				temp.f		<- .col[[i]]
				temp.f		<- temp.f[length(temp.f)]
				.cols[i] 	<- temp.f
			#ui$setProgressBar(pb, i*ratio.prog, label=paste("Cutting paths: ","Rf:",i))
			}

			#.cols   	<- sub("\\","", .cols,fixed = TRUE)
			#.cols 	<- sub(".mgf.htm","", .cols,fixed = TRUE)
			#.cols 	<- sub(".raw","", .cols,fixed = TRUE)
			#.cols 	<- sub(".msm.htm","", .cols,fixed = TRUE)

			for( i in 1:length(x.uni)){
				x[x==x.uni[i]] = .cols[i]
				
			}
			
			
			return(x)
		}
		
if(length(cbn.prot) != 0){
	
		
		if(length(grep(cbn.prot,x$code,fixed = TRUE)) == 0) {
try(x$Proteins			<- gsub("__","_",x$Proteins,fixed = TRUE))

try(x$Leading.Proteins	<- gsub("__","_",x$Leading.Proteins,fixed = TRUE))
try(x$code				<- gsub("__","_",x$code,fixed = TRUE))
try(x$Description		<- gsub("__","_",x$Description,fixed = TRUE))
cbn.prot <- tolower(cbn.prot)
			if(length(grep(tolower(cbn.prot),x$code,fixed = TRUE)) == 0) {

			 ui$messageBox(title="",message="Reference Protein is not represented in data.\nPlease chech your settings and rerun cRacker!",icon="error",type="ok") ;stop()}
		}
}


if(calc.empai){
	
exclude.data <- cbind(x$code,x$rawfilename)
exclude.data.uni <- unique(exclude.data)
exclude.data.agg<- aggregate(exclude.data.uni[,2],list(exclude.data.uni[,1]),length)

exclude.prot  <-exclude.data.agg[exclude.data.agg[,2 ]> length(unique(x$rawfilename))/1*gui.input$shape.prot.norm,]

include.vec <- grep(paste(exclude.prot[,1],collapse = "|"),x$code)
x <- x[include.vec,]	
	
}
#stop()
####
# cut paths 
####
#if(is.null(cbn.prot) == FALSE){
#x 		<- x[x$sam_id != exclude.raw,]
#}	

#stop()

result 	<- as.character(unique(x$sam_id))
if(length(result) == 0 & path.design == ""){
	alarm()
	cat("\nALERT: No experimental design available!\n")
	Raw 		<- TRUE
	x$sam_id 	<- rep("NA",dim(x)[1])
}else{
	

	
	# Setting exp.set
	## read own experimental design
if(path.design != ""){

		exp.set.2 		<- read.csv(path.design,sep = "\t", stringsAsFactors = F)
		exp.set.backup 	<- exp.set.2
		#
	
		
		exp.set.2[,1] <- tolower(cut.path(exp.set.2[,1]))
		exp.set.2[,1] <- make.names(tolower(exp.set.2[,1]),allow = FALSE)
		exp.set.2[,2] <- make.names(tolower(exp.set.2[,2]),allow = FALSE)
		exp.set.2[,1] <- make.names(tolower(exp.set.2[,1]),allow = FALSE)
		exp.set.2[,5] <- make.names(tolower(exp.set.2[,5]),allow = FALSE)
		
		if(any(!(exp.set.2[,1] == exp.set.2[,5]))){
			print("replacing raw file name")
			
			for(r in 1:dim(exp.set.2)[1]){
				x$rawfilename[tolower(x$rawfilename) == exp.set.2[r,1]] <- exp.set.2[r,5]
				
			}
			
			exp.group.shape	<- exp.set.2[,c(5,2,4)]			
			exp.set.2 <- exp.set.2[,c(5,2)]
	
		}else{
			exp.group.shape	<- exp.set.2[,c(1,2,4)]
			exp.set.2 <- exp.set.2[,1:2]

		}
		
		
		
		if(dim(exp.set.2)[2] == 1){
			exp.set.2 <- read.table(path.design,sep = "\t",stringsAsFactors = FALSE)
		}

		if(dim(exp.set.2)[2] == 1){
			exp.set.2 <- read.csv2(path.design,stringsAsFactors = FALSE)
		}
		
		.grep.exp <- grep("name|experiment",tolower(colnames(exp.set.2)))
		if(length(.grep.exp) <2){

		
			.grep.exp <- grep("Name",exp.set.2[,1])
			if(is.null(.grep.exp) == FALSE){
				colnames(exp.set.2) 	<- as.vector(exp.set.2[.grep.exp,])
				exp.set.2 			<- exp.set.2[-.grep.exp,]	
				if(dim(exp.set.2)[2] == 3){
					.grep.exp <- grep("Slice",colnames(exp.set.2))
					exp.set.2	<- exp.set.2[,-.grep.exp]
				}
			}
		
			.grep.exp <- grep("name|experiment",tolower(colnames(exp.set.2)))
			exp.set.2 <- exp.set.2[,.grep.exp]

		
		}else{
		
			exp.set.2 <- exp.set.2[,.grep.exp]
		}
		exp.set.2 <- apply(exp.set.2,2,as.character)
	

	##### change result file
	x$rawfilename <- cut.path(as.character(x$rawfilename))
	.result 	<- as.character(x$rawfilename)
	
	
	for( .r in 1:dim(exp.set.2)[1]){
		.result[tolower(.result) == exp.set.2[.r,1]]  <- exp.set.2[.r,2]
		
		
	}
	#stop()
	x$sam_id <- .result
		
	#stop()
	}else{	#		
	
	cat("\n cutting paths...\n")
	
	
	
	.col		<- strsplit(result,"\\",fixed = TRUE)
	.cols	<- c()

	ratio.prog <- prog.max/length(.col)
	for(i in 1: length(.col)) {
		temp.f		<- .col[[i]]
		temp.f		<- temp.f[length(temp.f)]
		.cols[i] 	<- temp.f
		x$sam_id[x$sam_id == result[i]] <- temp.f
		
		####### GUI
		
		pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Cutting paths: ","Rf:",i))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Cutting paths: ","Rf:",i))))
		}
		
		#########
	}

	#x$sam_id   	<- sub("\\","", x$sam_id,fixed = TRUE)
	#x$sam_id 	<- sub(".mgf.htm","", x$sam_id,fixed = TRUE)
	#x$sam_id 	<- sub(".raw","", x$sam_id,fixed = TRUE)
	}
	#x$sam_id 	<- sub(".msm.htm","", x$sam_id,fixed = TRUE)
	
	
	
	
	
}
#stop()
raw.f	<- unique(as.character(x$rawfilename))

.col		<- strsplit(raw.f,"\\",fixed = TRUE)
.cols	<- c()
ratio.prog <- prog.max/length(.col)
for(i in 1: length(.col)) {
	temp.f	<- .col[[i]]
	temp.f	<- temp.f[length(temp.f)]
	.cols[i] <- temp.f
	x$rawfilename[x$rawfilename == raw.f[i]] <- temp.f
	
			####### GUI
		
		pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Cutting paths: ","Raw:",i))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Cutting paths: ","Raw:",i))))
		}
		
		#########
	
}
		
#x$rawfilename   <- sub("\\","", x$rawfilename,fixed = TRUE)
#x$rawfilename 	<- sub(".mgf.htm","", x$rawfilename,fixed = TRUE)
#x$rawfilename 	<- sub(".msm.htm","", x$rawfilename,fixed = TRUE)
#x$rawfilename 	<- sub(".raw","", x$rawfilename,fixed = TRUE)

ratio.prog <- prog.max/3

	####### GUI
		
		pb.check <- class(try(ui$setProgressBar(pb, 0*ratio.prog, label=paste("Sorting of Proteins."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, 0*ratio.prog, label=paste("Sorting of Proteins."))))
		}
		
	#########
x 		<- 	x[order(x$code),]
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, 1*ratio.prog, label=paste("Sorting of rawfilename."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, 1*ratio.prog, label=paste("Sorting of rawfilename."))))
		}
	##############
x 		<- 	x[order(x$rawfilename),]
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, 2*ratio.prog, label=paste("Sorting of sam_id."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, 2*ratio.prog, label=paste("Sorting of sam_id."))))
		}
	##############
x 		<- 	x[order(x$sam_id),]
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, 3*ratio.prog, label=paste("Sorting of sam_id."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, 3*ratio.prog, label=paste("Sorting of sam_id."))))
		}
	##############
		



#temp.my$code <- gsub(" ","",temp.my$code)
	

rm(.col,.cols,result, raw.f) # cleaning workspace

######
# start calculation
######
		
box.ex <- NA

	#build.matrix <- TRUE
#***************************************************************************************************************
#***************************************************************************************************************
if(1==1){
	
	print("start calculation of IonIntensities!")
	temp.my <-as.data.frame(x, stringsAsFactors = FALSE)
	
	if(score > 0){
		temp.my[is.na(as.numeric(temp.my$Score)),] <- 0
		temp.my	<- 	temp.my[as.numeric(temp.my$Score) > score,]
			if(dim(temp.my)[1] == 0){
				ui$messageBox(title="",message="No peptides left, after setting score threshold. Please rerun cRacker with a lowered score threshold.",icon="error",type="ok")
				stop()	
			}
	}
	backup <- temp.my
	repe = 1






	if(build.matrix == FALSE){
			build.matrix.test <- grep("matrix-binary.Rdata",list.files(path2))	
			if(length(build.matrix.test) > 0){
				load(paste(path2,list.files(path2)[build.matrix.test[1]],sep ="/"))
			}else{
				build.matrix == TRUE; print("Error in loading matrix-binary.Rdata. \nRecalculating Matrix...")
			}
	}

	sam.id <- as.character(temp.my$sam_id)
	
		
	if(build.matrix == TRUE){	
		
		if(gui.input$graphic.type == "pdf"){
		pdf("LC.pdf")
		}
		for(p in 1:repe){	
		if(gui.input$graphic.type == "eps"){
		postscript("LC.eps")
		}
		
		
				temp.my <- backup
			
		##
		# exclude non identified Peptides
		###
			quantified.identifier <- as.numeric(temp.my$intensity.1)
			quantified.identifier[is.na(quantified.identifier)] <- 0
			if(calc.empai == FALSE){
				temp.my <- temp.my[quantified.identifier != -777.777,]
			}
			
			if(dim(temp.my)[1] == 0& calc.empai == FALSE){
				ui$messageBox(title="",message="Could not find ion intensities. \nPlease ensure that your data was quantitated!",icon="error",type="ok")
				stop()
			}
			cat("> Using IonIntensities for relative quantitation! \n")
			temp.my.raw			<- tolower(temp.my$rawfilename)
			temp.my..col.uni	<- unique(tolower(temp.my$rawfilename))
			temp.my.pepid		<- as.matrix(unique(paste(
									temp.my$code,
									temp.my$sequence,
									round(as.numeric(temp.my$mcr),digits = 0),
									temp.my$charge,sep = "#"
							)))
		
			pep.all.mean		<-  as.matrix(unique(paste(
									temp.my$code,
									temp.my$sequence,
									round(as.numeric(temp.my$mcr),digits = 0),
									temp.my$charge,sep = "#"
							)))
		
			raw.peptides	<-  as.matrix(unique(paste(
									temp.my$code,
									temp.my$sequence,
									round(as.numeric(temp.my$mcr),digits = 0),
									temp.my$charge,sep = "#"
							)))
		
		
			if(calc.empai == TRUE){
				pep.all.mean <- as.matrix(unique(temp.my$code))
				raw.peptides <- pep.all.mean
				temp.my.pepid <- pep.all.mean
			}

			#modified <- temp.my$Modifications[temp.my$Modifications!= ""]
			###
			# 1. create peptide matrix, and merge same peptides in sample
			###
			cat("\n")
			cat(			"****************************************************************************************\n")
			cat(paste(		"************** 1. start creating peptide matrix*****************************************\n"))
			cat(			"****************************************************************************************\n")
			cat("\n")
	
			double.pep <- c()
			
		
				text.out.add<-""
			
			iterator <- 0
			for(s in 1:length(unique(sam.id))){
				print("s-loop")
				temp.s	<- temp.my[tolower(temp.my$sam_id) == tolower(unique(sam.id))[s],]
				total2	<- length(unique(tolower(temp.s$rawfilename)))
		
				text.out <- paste("*** Resultfile",s,"of",length(unique(sam.id)), "containing" , length(unique(tolower(temp.s$rawfilename))),"raw files", "***",.collapse = "")
				cat(paste(rep("*",nchar(text.out)),.collapse = "",sep = ""),"\n")
				cat(text.out,"\n")
				cat(paste(rep("*",nchar(text.out)),.collapse = "",sep = ""),"\n")
				
				if(calc.empai == FALSE){
				
						temp.my.pepid <- unique(paste(
						temp.my$code,
						temp.my$sequence,
						round(as.numeric(temp.my$mcr),digits = 0),
					temp.my$charge,sep = "#"
					))
					
					
				}else{
					temp.my.pepid <- unique(temp.my$code)
				}
				#############################################################################################################################
				# -- rawfiles: --
				if(calc.empai == FALSE | Raw == TRUE){
					for(i in 1 :length(unique(tolower(temp.s$rawfilename)))){
						temp.i	 <- temp.s[tolower(temp.s$rawfilename) == unique(tolower(temp.s$rawfilename))[i],]
					if(calc.empai == FALSE){
						###### LC
						lc 		<- 	as.numeric(temp.i$Retention.Time)
						lc.max	<-	as.numeric(temp.i$intensity.1)
						mcr		<-	round(as.numeric(temp.i$mcr),digit = 0)
						
						## .colorramp start
						lc.max <- lc.max/max(lc.max,na.rm = TRUE)*1000
						lc.size <- log10(lc.max/4)
						lc.size[lc.size < 0.5] 	<- 0.5
						lc.size[is.na(lc.size)] <- 0.1
		
						
						test <- as.matrix(cbind(lc.max,c(1:length(lc.max))))
						test <- test[order(test[,1]),] 
						if(is.vector(test)){
							test <- rbind(test,rep(NA,length(test)))
						}
	
						pal.crp	<- colorRampPalette( c("#0000CD99","green","orange","red"),bias = 0.5) 
					
						if(all(is.na(unique(test[,1]))) == FALSE){
							.col.test 		<- pal.crp(max(test[,1],na.rm = TRUE))
							test[,1] 		<- abs(round(test[,1],digit = 0)+1)
							test[test[,1]==0,] <- 1
							test[,1]		<- .col.test[(test[,1])]
							
							
							test 			<- test[order(as.numeric(test[,2])),] 
							.col.ramp <- test[,1]
							## .colorramp end
							lc.max <- lc.max / mean(lc.max,na.rm = TRUE)
							
							## plot start 
							library(grDevices)
									layout(matrix(c(4,1,2,0,3,0),2,3,byrow=TRUE) ,c(0.12,1,0.1), c(1,0.1), FALSE)	
							
							add.main = ""
							
							
							if(is.null(mcr)== FALSE){
								par(mar = c( 2.5, 3,4,1))
		
								if(length(lc) !=0 ){
							mcr[is.na(mcr)] <- 0

									plot(lc,mcr,pch=20,cex = 1,col = .col.ramp,ylab = "MCR",xlab = "time in s",main = paste(add.main,unique(tolower(temp.s$rawfilename))[i]),mgp = c(1.5,0.5,0),xaxs = "r",	type	 = "n")
						
									points(lc,mcr,pch = 20,cex =  lc.size,col = .col.ramp,)
									axis(3,mgp = c(2,0.5,0))
									axis(4,mgp = c(2,0.5,0))
							
									par(mar = c( 2.5,0,4,0))
									boxplot(mcr,bty = "n",axes = "FALSE")
									par(mar = c( 0,3,0,0.5))
						
									boxplot(lc, horizontal = TRUE, bty = "n",axes = "FALSE",lwd =1)
									


									x<- c(0.04,0.052,0.057,0.07)
									y<- c(0.9,0.15,0.15,0.9)
									.col.ramp2 <- pal.crp(20)
									yy <- seq(range(y)[1],range(y)[2], length=length(.col.ramp2))
									dx <- diff(yy)
									for(ii in 1:length(yy)){
										#print(ii)
										grid.clip(x=0, y=yy[ii],
										width= 1, 
										height=2*dx[ii],
										just="bottom"
										)
										grid.polygon(x=x, y=y, gp=gpar(fill=((.col.ramp2))[ii]))
									}
									frame()
									mtext("intensity color code",2,outer = F,padj = -1.9,cex = 0.75)
									
																		mtext(paste("low",paste(rep(" ",130),collapse = ""),"high","       ",collapse = ""),4,outer = F,padj = -1.1,cex = 0.75)

									grid.segments(	c(0.01,0.1,0.01,0.01),
												c(.l <- 0.925,.l,.l, .l2<- 0.137),
												c(0.01,0.1,0.1,0.1),
												c(.l2,.l2,.l,.l2))

								}else{
									
									
									layout(matrix(c(4,1,2,0,3,0),2,3,byrow=TRUE) ,c(0.12,1,0.1), c(1,0.1), FALSE)	
									par(mar = c( 2.5, 3,4,1))
			
									plot(mcr,pch=20,cex = 1,col = .col.ramp,ylab = "MCR",main = paste(add.main,unique(tolower(temp.s$rawfilename))[i]),mgp = c(1.5,0.5,0),xaxs = "r",type = "p")			
									points(mcr,pch = 20,cex =  lc.size,col = .col.ramp)
									mtext("Retention time not available!")
									axis(4,mgp = c(2,0.5,0))
									par(mar = c( 2.5,0,4,0))
									boxplot(mcr,bty = "n",axes = "FALSE")
									
									x<- c(0.04,0.052,0.057,0.07)
									y<- c(0.9,0.15,0.15,0.9)
									.col.ramp2 <- pal.crp(20)
									yy <- seq(range(y)[1],range(y)[2], length=length(.col.ramp2))
									dx <- diff(yy)
									for(ii in 1:length(yy)){
										#print(ii)
										grid.clip(x=0, y=yy[ii],
										width= 1, 
										height=2*dx[ii],
										just="bottom"
										)
										grid.polygon(x=x, y=y, gp=gpar(fill=((.col.ramp2))[ii]))
									}
									frame()
									mtext("intensity color code",2,outer = F,padj = -1.9,cex = 0.75)
									
																		mtext(paste("low",paste(rep(" ",130),collapse = ""),"high","       ",collapse = ""),4,outer = F,padj = -1.1,cex = 0.75)

									grid.segments(	c(0.01,0.1,0.01,0.01),
												c(.l <- 0.925,.l,.l, .l2<- 0.137),
												c(0.01,0.1,0.1,0.1),
												c(.l2,.l2,.l,.l2))

									}
							}
						}
						
				## plot end 
				#LC end
				######
						text.out <-paste("Starting with rawfile",unique(temp.my.raw)[i],.collapse = "")
						cat(paste(rep("*",nchar(text.out)),.collapse = "",sep = ""),"\n")
						cat(text.out,"\n")
						cat(paste(rep("*",nchar(text.out)),.collapse = "",sep = ""),"\n")
						
					
							temp.i.pepid <- paste(
							temp.i$code,
							temp.i$sequence,
							round(as.numeric(temp.i$mcr),digits = 0),
							temp.i$charge,sep = "#"
							)
						
# duplication, parser
						dupli.vec 	<- duplicated(temp.i.pepid) # find duplicated entries
						dupli 		<- temp.i.pepid[dupli.vec] # string of duplicates
						uni.dupli 	<- unique(dupli)	# unique string
		
						temp.i.intensities				<- temp.i$intensity.1[dupli.vec == FALSE]
						names(temp.i.intensities)		<- temp.i.pepid[dupli.vec == FALSE]
		
					
		
		
						
						temp.d.uni <- c()
						print(paste("going through",length(uni.dupli),"peptides"))
dupli.vec 	<- duplicated(temp.i.pepid) # find duplicated entries
						dupli 		<- temp.i.pepid[dupli.vec] # string of duplicates
						uni.dupli 	<- unique(dupli)	# unique string
		
						temp.i.intensities	<- temp.i$intensity.1[!is.element( temp.i.pepid, uni.dupli)]
						
						names(temp.i.intensities)		<- temp.i.pepid[!is.element( temp.i.pepid, uni.dupli)]
		
								
				
						temp.d.uni <- c()

						temp.aggregate.peptides <- cbind(temp.i.pepid,temp.i$intensity.1,temp.i	$Retention.Time)
						colnames(temp.aggregate.peptides )<- c("peptide.identifier","intensity","retention.time.in.min")
						temp.aggregate.peptides <- temp.aggregate.peptides[is.element( temp.i.pepid,uni.dupli),]
						temp.aggregate.peptides2 <- c()
						if(dim(temp.aggregate.peptides)[1] >0){
						temp.aggregate.peptides2	<- aggregate(temp.aggregate.peptides[,2],list(temp.aggregate.peptides[,1]),function(x){
				
							x <- as.numeric(as.character(x))
			
							if(action.dupli == "mean"& !all(is.na(x))){
								x <- mean(x,na.rm = TRUE)
							}
							if(action.dupli == "max"& !all(is.na(x))){
								x <- max(x,na.rm = TRUE)
							}
							if(action.dupli == "min"& !all(is.na(x))){
								x <- min(x,na.rm = TRUE)
							}
							if(action.dupli == "sum"& !all(is.na(x))){
								x <- min(x,na.rm = TRUE)
							}
							if(action.dupli == "exclude"& !all(is.na(x))){
								x <- mean(x,na.rm = TRUE)
							}	
				
								return(as.numeric(unique(x)))
						}

						)
						
						

						if(action.dupli !=  "exclude"){
							if(is.vector(temp.aggregate.peptides2)){
								temp.aggregate.peptides2 <- as.matrix(temp.aggregate.peptides2)
							}
							add.vec			   	<-	unlist(temp.aggregate.peptides2[,2])
							names(add.vec)		<-	temp.aggregate.peptides2[,1]
							temp.i.intensities <-	c(temp.i.intensities,add.vec) 													
						}
					#if(any(is.na(names(temp.i.intensities)))){stop()}
																
						dir.create("duplicates")	
						setwd("./duplicates")		
										
						if(length(temp.i.pepid[dupli.vec == TRUE]) != 0 ){
										write.csv(temp.aggregate.peptides,file = paste("duplicates-",unique(tolower(temp.s$rawfilename))[i],".csv",sep =""))
						}
						setwd("../")	
						
						}
						
						


	########## GUI 
	ratio.prog <- prog.max/length(unique(temp.my$rawfilename))

	iterator <- iterator +1
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*iterator, label=paste(text.out.add,"Creating matrix: ","Rf:",s,"/",length(unique(sam.id)), ".column: ",i,"/",length(unique(tolower(temp.s$rawfilename)))))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*iterator, label=paste(text.out.add,"Creating matrix: ","Rf:",s,"/",length(unique(sam.id)), ".column: ",i,"/",length(unique(tolower(temp.s$rawfilename)))))))
		}
	##############
									
						temp.a.pep <- cbind(names(temp.i.intensities),temp.i.intensities)
						colnames(temp.a.pep) <- paste(i,colnames(temp.a.pep))
						
temp.my.pepid <- merge(as.matrix(temp.my.pepid),temp.a.pep,by = 1,all = TRUE)
					
						print("done")
				
			}else{
						print("counting n")
						#empai for raw
						print("empai for every sample")
						if(empai.from.msms & !is.null(temp.i$MS.MS.Count)){
							prot.tryp.length 	<- aggregate(temp.i$MS.MS.Count,list(temp.i$code),length)
							prot.tryp.length1 	<- aggregate(temp.i$code,list(temp.i$code),length)
							
						}else{
								prot.tryp.length 	<- aggregate(temp.i$code,list(temp.i$code),length)
							if(!is.null(temp.i$MS.MS.Count)){
								prot.tryp.length1 	<- aggregate(temp.i$MS.MS.Count,list(temp.i$code),length)
							}
						}	
		
						#print(dim(prot.tryp.length))
						temp.my.pepid <- merge(as.matrix(temp.my.pepid), prot.tryp.length,by = 1,all = TRUE)
				
			}			
	}
	}else{
		#empai for exp
		print("empai of merged replicates")
		
		
		if(empai.from.msms& !is.null(temp.s$MS.MS.Count)){	
			prot.tryp.length 	<- aggregate(temp.s$MS.MS.Count,list(temp.s$code),length)
			prot.tryp.length1 	<- aggregate(temp.s$code,list(temp.s$code),length)

		}else{
			prot.tryp.length 	<- aggregate(temp.s$code,list(temp.s$code),length)
			if(!is.null(temp.s$MS.MS.Count)){
				prot.tryp.length1 	<- aggregate(temp.s$MS.MS.Count,list(temp.s$code),length)
			}
		}
		#print(dim(prot.tryp.length))
		temp.my.pepid <- merge(as.matrix(temp.my.pepid), prot.tryp.length,by = 1,all = TRUE)
				
	}

		
			# -- create .colnames / cut raw file path --
			.cols	<- unique(tolower(temp.s$rawfilename))
			colnames(raw.peptides)	<- c(1:dim(raw.peptides)[2])
			colnames(temp.my.pepid) <- c(1:dim(temp.my.pepid)[2])
			raw.peptides 			<- merge(raw.peptides, temp.my.pepid,by = 1, all = TRUE)
			rows					<- temp.my.pepid[,1]
			temp.my.pepid			<- as.matrix(temp.my.pepid[,2:dim(temp.my.pepid)[2]])
			rownames(temp.my.pepid) <- rows
			temp.my.mean			<- temp.my.pepid
			temp.my.mean 			<- cbind(rownames(temp.my.mean),temp.my.mean)
			# -- merge of raw-data --
			colnames(pep.all.mean)	<- c(1:dim(pep.all.mean)[2])
			colnames(temp.my.mean)	<- c(1:dim(temp.my.mean)[2])
			pep.all.mean			<- merge(pep.all.mean,temp.my.mean,by = 1,all = TRUE)
	
	}
		
		
		# -- For constructing experiment --
		
		
		graphics.off()# LC.pdf

		save.image(file="../matrix-binary.Rdata")

		rm(x,temp.my.mean, temp.my..col.uni, temp.my.pepid, temp.my.raw, temp.re,test, total2,temp,rows,s,result, quantified.identifier,lc.size,lc,.col.test,.col.ramp,.cols,lc.max)# cleaning namespace


		result		<- unique(tolower(as.character(temp.my$sam_id)))

		exp.set		<- c()
		for(re in 1:length(result)) {
			temp.re <- tolower(temp.my$rawfilename)[tolower(temp.my$sam_id )== result[re]]
			temp.re <- unique(temp.re)
			temp.re <- cbind(temp.re,rep(result[re],length(temp.re)))
			exp.set <- rbind(exp.set,temp.re)
		}
		.cols		<- exp.set[,1]
		exp.set 	<- cbind(.cols,exp.set[,2])
		.col.all 	<- c()
		
		for(tres in 1:length(unique(temp.my$sam_id))) {
			temp	<- temp.my[temp.my$sam_id == unique(temp.my$sam_id)[tres],]
			.cols	<- as.character(unique(temp$rawfilename))
			.col.all		<- c(.col.all,.cols)
		}
		.col.all			<- unique(tolower(.col.all))
	

	}


	save(pep.all.mean,exp.set,.col.all,file="../matrix-binary.Rdata")
#
}
			

	###
	# 2. merge peptide data to protein data
	###
	cat("\n")
		  cat("****************************************************************************************\n")
	cat(paste("************** 2. start merging peptide information to protein level********************\n"))
		  cat("****************************************************************************************\n")
	cat("\n")
	
	#stop()
	cbn.prot.data 			<- c() # for BSA matrix
	#### Naming:
	rows					<- pep.all.mean[,1]
	pep.all.mean			<- as.matrix(pep.all.mean[,2:dim(pep.all.mean)[2]])
	pep.all.mean			<- apply(pep.all.mean,2,as.numeric)
	
	
	

	#print(dim(pep.all.mean))
	#print(unique(sam.id))
	#print(length(rows))
	
	
	####
	## red. peptide new location 
	####
	
	#print(calc.empai == "TRUE"&Raw == FALSE)
	if(calc.empai == "TRUE"&Raw == FALSE){
		rownames(pep.all.mean)	<- rows
		colnames(pep.all.mean)	<- unique(sam.id)
	
	}else{
		rownames(pep.all.mean)	<- rows
		colnames(pep.all.mean)	<- .col.all
	}
	pep.all.mean[pep.all.mean == 0] <- zero.treat
	
	if(calc.empai){
		write.csv(pep.all.mean,"n-peptides.csv")
	}
	rpn.start.peptide.matrix <- pep.all.mean
	
	
	
	if(calc.empai == "TRUE"){
		print("calculating empai")
	####
	# empai calculations
	####
		rownames(pep.all.mean) <- gsub(" ","",rownames(pep.all.mean))
		
		write.csv(pep.all.mean,"n-empai.csv")

		
		.empai <- c()
		empai.aov <- c()
		ratio.prog <- prog.max/dim(pep.all.mean)[1]
		
		for(o in 1:dim(pep.all.mean)[1]){	

			
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb,o*ratio.prog, label=paste("Calculating emPAI",round(as.numeric(o)/dim(pep.all.mean)[1]*100, 0),"% done"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb,o*ratio.prog, label=paste("Calculating emPAI",round(as.numeric(o)/dim(pep.all.mean)[1]*100, 0),"% done"))))
		}
	##############
		
			temp.o 			<- grep(rownames(pep.all.mean)[o],.length.matrix[,2])
			temp.o.empai 	<- 10^(as.numeric(pep.all.mean[o,])/as.numeric(.length.matrix[o,(length.pep-1)]))-1
			.empai			<- rbind(.empai,temp.o.empai)
			temp.aov.data 	<- cbind(rownames(pep.all.mean)[o],colnames(pep.all.mean), temp.o.empai)
			empai.aov 		<- rbind(empai.aov, temp.aov.data)	
		}	
		
		colnames(empai.aov) = c("accession","experiment","intensity")

		rownames(.empai) <- rownames(pep.all.mean)
		#stop()
		sam.mean <- .empai
		if(Raw == FALSE){
			colnames(sam.mean) <- unique(sam.id)
		}else{
			colnames(sam.mean) <- .col.all
		}
		
		write.csv(sam.mean,"raw-empai.csv")
		
	####
	# empai normalization
	####
		
	n.sample <- apply(.empai,2,function(x){length(x[!is.na(x)])})
		
		
		if(empai.norm == "sum" & is.null(gui.input$cbn.prot)){
		sam.mean <-	hz.norm(sam.mean,2,norm = sum.of.total)$x
		}		
		if(empai.norm == "n"& is.null(gui.input$cbn.prot)){		
		sam.mean <- t(apply(sam.mean,1,function(x){x/n.sample}))
		}
		
		if(!is.null(gui.input$cbn.prot)){
			cbn.target.grep <- grep(gui.input$cbn.prot,rownames(sam.mean))
			cbn.target 		<- sam.mean[cbn.target.grep,]
			if(any(cbn.target == 0,is.na(cbn.target))){
				cbn.target[is.na(cbn.target)|cbn.target == 0] <- mean(cbn.target[!is.na(cbn.target)|cbn.target == 0])
			 cancel <- tkmessageBox(type = "okcancel",message = paste("No entry of reference protein found for sample(s):\n",paste(colnames(sam.mean)[is.na(cbn.target)|cbn.target == 0],collapse = "\n"),"\nUse average reference protein value?",collapse =  ""))
			 if(as.character(cancel)== "cancel"){stop("User stopped calculation")}
			
			}
			
			
			
			sam.mean <- t(apply(sam.mean,1,function(x){x/cbn.target}))
			
			
		}
		
		
		
		shape.loop <- 1
		while(shape.loop == 1){
			sam.mean.shape <- hz.shape(sam.mean,shape = shape,group.shape)$shape
			if(length(sam.mean.shape)< (2*dim(sam.mean)[2]-1) & shape != 1){
				shape <- shape -0.1
				gui.input$shape <- shape -shape*100
			}
			
			if(is.vector(sam.mean.shape)){
				shape <- shape -0.1
				gui.input$shape <- shape -shape*100
			}
			
			if(length(sam.mean.shape)< (2*dim(sam.mean)[2]-1) & shape == 1){
				
				 ui$messageBox(title="",message="No data points left after shahping!",icon="error",type="ok") ;stop()
			}
			if(length(sam.mean.shape)> (2*dim(sam.mean)[2]-1)){
				sam.mean<- sam.mean.shape	
				shape.loop <- 0
			}
			
		}		

	
		##############################################################################
		########## AOV.EXPORT.CREATION ###############################################
		##############################################################################
		
		if(empai.sd == TRUE & Raw == TRUE){exp.set.type <- "exp"}
		if(empai.sd == FALSE & Raw == TRUE){exp.set.type <- "raw"}
		if(empai.sd == FALSE & Raw == FALSE){exp.set.type <- "exp"}
		
		if(exp.set.type == "raw"){
			template.exp.set <- unique(exp.set[,1])
		}else{
			template.exp.set <- unique(exp.set[,2])	
		}
		
		
		for(exp.set.i in 1:length(template.exp.set)){
			if(empai.sd == TRUE & Raw == TRUE){
				aov.export <- c()
				for(change.name in template.exp.set){
					temp.change.name <- exp.set[exp.set[,2] == change.name ,1]
					temp.int <- c()
						for(temp.change.name.raw in temp.change.name){
							temp.int <- c(temp.int,sam.mean[,colnames(sam.mean)== temp.change.name.raw])	
						}
					change.name.output 	<- cbind(names(temp.int),change.name,temp.int)
					aov.export 			<- rbind(aov.export, change.name.output)
					
				}
				
			}else{
				aov.export <- c()
				
				for(column.size in 1:dim(sam.mean)[2]){
		
					temp.aov.export 			<- cbind(rownames(sam.mean),colnames(sam.mean)[column.size],sam.mean[,column.size])
					aov.export <- rbind(aov.export, temp.aov.export)
				
				}
				
				
			}
			
			
		}
		
		aov.export <- aov.export[!is.na(aov.export[,3]),]
		#print(dim(aov.export))
		
		####################################################################################
		####################################################################################
		####
		# create sds for empai
		####	
		sam.sd	<- matrix(NA,nrow = dim(sam.mean)[1],ncol = dim(.empai)[2])
		
		if(empai.sd == TRUE & Raw == TRUE){
					
			sam.sd.ia <- c()

			sam.mean.ia <- c()
			ratio.prog <- prog.max/length(unique(exp.set[,2]))
			for(ia in 1:length(unique(exp.set[,2]))){
		
				
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb,o*ratio.prog, label=paste("Calculating emPAI sd",round(as.numeric(ia)/ia*100, 0),"% done"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb,o*ratio.prog, label=paste("Calculating emPAI sd",round(as.numeric(ia)/ia*100, 0),"% done"))))
		}
	##############
		
				
				
				temp.ia		<- 	exp.set[exp.set[,2] == unique(exp.set[,2])[ia],]
				#print(temp.ia)
				temp.ib.grep <- c()
				if(is.vector(temp.ia)){temp.ia <- t(as.matrix(temp.ia))}
				for(ib in 1:length(temp.ia[,1])){
					temp.ib <- grep(temp.ia[,1][ib],colnames(sam.mean))
					temp.ib.grep[ib] <- temp.ib
				}
					
				temp.ia.m 	<- sam.mean[,temp.ib.grep]
			
			
			
			
				if(merge.method == "mean"){
					temp.ia.mean <- apply(as.matrix(temp.ia.m),1,function(x){mean(x,na.rm =TRUE)})
				}
				if(merge.method == "median"){
					temp.ia.mean <- apply(as.matrix(temp.ia.m),1,function(x){median(x,na.rm =TRUE)})
				}

				temp.ia.sd <- apply(as.matrix(temp.ia.m),1,function(x){sd(x,na.rm =TRUE)})
				temp.ia.sd <- temp.ia.sd/temp.ia.mean
	
				sam.sd.ia 	<- cbind(sam.sd.ia, temp.ia.sd)
				sam.mean.ia 	<- cbind(sam.mean.ia, temp.ia.mean)
	
			}	
		
			colnames(sam.sd.ia) <- unique(exp.set[,2])
			colnames(sam.mean.ia) <- unique(exp.set[,2])
	
			sam.mean <- sam.mean.ia
			sam.sd	 <- sam.sd.ia
	
			if(exists("gui.input")){gui.input$raw <- "FALSE"}
		}
		
	###
	# fill needed objects
	###
	
		write.pep.all.mean.n	<- matrix(0,nrow = dim(sam.mean)[1],ncol = dim(.empai)[2])
		all.n..col				<- pep.all.mean
			
			if(group.filter & gui.input$exp.design != ""){
				group.filter.order <- hz.merge.control(tolower(make.names(exp.group.shape[,1],allow = F)),colnames(all.n..col))
				group.shape <- exp.group.shape[group.filter.order,3]
		
			}else{
				group.shape <- rep(1,dim(sam.mean)[2])
			}
		
		all.n..col				<- hz.shape(all.n..col,shape = shape, group.shape)$shape


		
		aov.export.1			<- aov.export
		aov.export.2			<- aov.export
		
			colnames(aov.export.1) = c("accession","experiment","intensity")
			colnames(aov.export.2) = c("accession","experiment","intensity")

		temp.e.mod.mean.m		<- c()

	if(add.data == TRUE){
			###### info - data		
		info.data <- 	as.data.frame(cbind(
				as.character(temp.my$code),
				as.character(temp.my$Description),
				temp.my$Score,
				temp.my$mcr,
				temp.my$charge,
				temp.my$Calibrated.mass.relative.error..ppm,
				as.character(temp.my$sequence)
				))
		.colnames.info.data <- as.character(c(
			"code",
			if(length(temp.my$Description)>0){ "Description" },
			if(length(temp.my$Score)>0){ "Score" },
			if(length(temp.my$mcr)>0){ "mcr" },
			if(length(temp.my$charge)>0){ "charge" },
			if(length(temp.my$Calibrated.mass.relative.error..ppm)>0){ "Calibrated.mass.relative.error..ppm" },
			if(length(temp.my$sequence)>0){ "sequence" }
		))
		
		all.peptides <- as.matrix(
				temp.my$code)
			
		
		
	info.data <- cbind(all.peptides,info.data)
	colnames(info.data) <- c("id",.colnames.info.data)
	dec.vec <- c()



		info.data.protein 	<- c()
		total 				<-  length(rownames(sam.mean))
			info 				<- function(x){
								if(as.numeric(x[1]) %%10==0 & exists("pb")) {
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb,as.numeric(x[1])*ratio.prog, label=paste("Extracting info: ",round(as.numeric(x[1])/total*100, 0),"% done"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb,as.numeric(x[1])*ratio.prog, label=paste("Extracting info: ",round(as.numeric(x[1])/total*100, 0),"% done"))))
		}
	##############
		
									
								}
								
								x 		<- as.character(x[2])
								temp	<- grep(x,as.character(info.data[,1]),fixed = TRUE)
								temp.my.z		<- info.data[temp,]
								
								temp.score		<- as.numeric(as.vector(temp.my.z$Score[]))
								temp.my.z 		<- temp.my.z[as.numeric(as.vector(temp.score)) == max(as.numeric(as.vector(temp.score)),na.rm = TRUE),]
								temp.my.z  		<- temp.my.z[!is.na(temp.my.z$id),]
								
								if(dim(temp.my.z)[1] > 1 ){
									mass.accuracy	<- as.vector(temp.my.z$Calibrated.mass.relative.error..ppm)
									mass.accuracy <- as.numeric(as.vector(mass.accuracy)) == min(as.numeric(as.vector(mass.accuracy)),na.rm = TRUE)
									mass.accuracy[is.na(mass.accuracy)] <- FALSE
									
									
									temp.my.z 		<- temp.my.z[mass.accuracy,]
									temp.my.z		<- as.vector(as.matrix(temp.my.z))
								}
								#if(is.vector(temp.my.z) == FALSE){print(x)}
								if(length(temp.my.z) != 8){alarm();print(x);temp.my.z <- rep("error in data collection!",8)}
								
								return(as.vector(as.matrix(temp.my.z)))
							}
			
		info.temp 			<- cbind(c(1:dim(sam.mean)[1]),as.matrix(rownames(sam.mean)))
		ratio.prog 			<- prog.max/dim(info.temp)[1]
		info.data.protein 	<- t(apply(info.temp ,1,info))
		
				try(colnames(info.data.protein ) <- colnames(info.data)		)


	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, prog.max, label=paste("Extracting info: ","100 % done"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, prog.max, label=paste("Extracting info: ","100 % done"))))
		}
	##############
		
		
		info.data.protein 			<- rep("no.data..collected",dim(sam.mean)[1])

		try(info.data.matrix 	<- cbind(info.data.protein,sam.mean))
	
	}
	
	
	
	

	}else{
		print("Using IonIntensities")
		#### Normalisation:
		if( length(cbn.prot) == 0 &  N15 == FALSE){
			if(use.raw == FALSE){
				pep.all.mean.n			<- hz.norm(pep.all.mean,margin = 2,norm = sum.of.total)$x
	
				
					if(n.correction){
						hz.n.correction<- function(x){	
							factor.vec <- apply(x,2,function(x){x <- x[!is.na(x)];x <-  x[!x==0];length(x)})
							factor.vec <- factor.vec/max(factor.vec,na.rm =TRUE)
							x<- t(apply(x,1,function(x){x*factor.vec}))
							return(list(x=x,factor=factor.vec))
						}
						pep.all.mean.n		<- hz.n.correction(pep.all.mean.n)$x
						pep.all.mean.n.factor		<- hz.n.correction(pep.all.mean.n)$factor
	
					
					}
				
				
			}else{		
				pep.all.mean.n	 		<- pep.all.mean
				target.mean <- 1
	
			}
					
		}else{
			#### non N15:
			if(N15 == FALSE ){
				norm.type			<- rep("abs",dim(pep.all.mean)[2])
				target				<- grep(cbn.prot,rows,fixed = TRUE)
				
				if(length(target) == 0){
					print("ALERT: Target not found!")
				}
				#### Reference Protein:	
				target			<- pep.all.mean[target,]
				write.csv(target,file = "BSA-peptides.csv")
				target.bsa.pep 	<- target # backup
				if(script.shape){
					#input <- target
					#source(paste(path1,"scripts/script.shape.R",sep = ""))
					#input <- script.shape.fun(input)
					#target <- input$shape
				}else{
											
					target.shape	<- hz.shape(target,shape = prot.norm.shape)$shape # matrix with all Peptides of Reference protein
					# Loop to reduce prot.norm.shape, if target == 0
					while(length(target.shape) == 0){
						cat("\n prot.norm.shape to high!\n")
						prot.norm.shape  <- prot.norm.shape -0.05
						target.shape <- hz.shape(target,shape = prot.norm.shape)$shape
						print(paste("Reduced prot.norm.shape to", prot.norm.shape))
						input.add1 <- paste("\n Changed prot.norm.shape to",prot.norm.shape,".\n")
					}
					if(exists("inpud.add1")){gui.input <- c(gui.input,list(Warning= input.add1))
					}
					target <- target.shape
						if(is.null(target)){target <- rep(1,length(.col.all))
						}
				}
				write.csv(target,file = "CBN-peptides-shape.csv")
	
				test 			<- hz.norm(target,1,norm = norm.method)$x
				write.csv(target,file = "CBN-peptides-shape-norm.csv")
	
				if(merge.method == "mean"){
					target.mean		<- apply(test,2,function(x){mean(x,na.rm = TRUE)})
				}
				if(merge.method == "median"){
					target.mean		<- apply(test,2,function(x){median(x,na.rm = TRUE)})
				}
				write.csv(target.mean,file = "CBN-protein.csv")
				
				
				target.sd		<- apply(test,2,function(x){sd(x,na.rm = TRUE)})
				write.csv(target.sd,file = "CBN-sd.csv")
				target.sd.rel	<- apply(test,2,function(x){sd(x,na.rm = TRUE)/mean(x,na.rm = TRUE)})
				target.n		<- apply(test,2,function(x){length(x)})
			
				cbn.prot.data 				<- rbind(target.mean,target.sd,target.n, target.sd.rel,norm.type)
				colnames(cbn.prot.data) 	<- colnames(test)
				rownames(cbn.prot.data) 	<- c("mean","sd","n","sd.relative","norm.type")
				
				range.y <- range(c(target.mean,c(target.mean+target.sd)),na.rm = TRUE)
				## LC	
				# plot distribution of bsa peptides
	if(gui.input$graphic.type == "pdf"){
				pdf("CBN-distribution.pdf",family = "Palatino",pointsize = 14)
					par(mar = c(7,3.5,3,0))
				}else{
				postscript("CBN-distribution.pdf",family = "Palatino",pointsize = 14)
					par(mar = c(7,3.5,3,0))

					
				}				
					library(gplots)
					bp<- barplot(target.mean,las = 2, border = NA,col = "transparent",ylim = range.y,mgp = c(2.3,1,0),ylab = "normalised intensity",main = paste(cbn.prot))
					plotCI(bp,target.mean,ui = target.mean+target.sd,li = target.mean-target.sd,add = TRUE )
				graphics.off()
				# test if protein has been measured in all samples
				test <- grep("TRUE",is.na(target.mean),fixed = TRUE)
				if(length(test)!=0){
					cat("alert: Target protein is not existend in all rawfiles!")
					target.mean[test] <- sum(pep.all.mean[,test])
					cbn.prot.data[5,test] <- "rel"
				}		
				# Normalise based on Ref-Protein
				if(row.norm == FALSE|row.target.norm == FALSE){
					pep.all.mean.n				<- t(apply(as.matrix(pep.all.mean),1,function(x){as.numeric(x)/as.numeric(target.mean)}))
					rownames(pep.all.mean.n) 	<- rownames(pep.all.mean)
					colnames(pep.all.mean.n)	<- colnames(pep.all.mean)
				}else{
					pep.all.mean.n <- pep.all.mean
				}
				
				colnames(pep.all.mean)	<- colnames(pep.all.mean)
			}else{
				
				# N15 Normalization
				
				if(n15.log2){
	
				n15.correct.expect <- log2(n15.correct.expect)	
				n15.log <- log2(pep.all.mean)
				n15.log[n15.log == Inf] 	<- 9999 
				n15.log[n15.log == -Inf] 	<- -9999
				
			if(n15.correct.method == "none"){
				pep.all.mean <- n15.log
			}	
			}
			
			if(n15.correct.method != "none"& n15.log2){
			
				print("start 15N")
				
				pdf("15N-correction.pdf")
				
				n15.correct.value <- c()
				n15.correct.matrix <- c()
				

				for(.n15 in 1: dim(n15.log)[2]){
					temp.n15 <- n15.log[,.n15]
			
					if(n15.correct.method == "median"){ n15.mean <- median(temp.n15,na.rm = TRUE)}
					if(n15.correct.method == "mean"){ n15.mean <- mean(temp.n15,na.rm = TRUE)}
					temp.correct.value <- n15.mean-	n15.correct.expect
					temp.n15.correct <- temp.n15-	temp.correct.value 
					n15.correct.value <- c(n15.correct.value,temp.correct.value )
					n15.correct.matrix <- cbind(n15.correct.matrix,temp.n15.correct)
						
						
					#hist(temp.n15,xlim = range(temp.n15, temp.n15.correct,na.rm = TRUE),col = "#00009f90", freq = FALSE)
					plot.input <- c(1,2,3,4)
					#print(temp.n15)
					try(plot.input	<- density(temp.n15,na.rm = TRUE))
					
					plot(plot.input,xlim = range(temp.n15, temp.n15.correct,na.rm = TRUE),main=paste(colnames(n15.log)[.n15],"\nCorrection of log2 data"),type = 		"l",col = "blue",lwd =3)
		
					#hist(temp.n15.correct,add = TRUE,col = "#99000090",freq = FALSE)
					points(plot.input,main="Density estimate of data",type = "l",col = "red",lwd =3)
					
		
					legend("topleft",c("uncorrected","corrected"),fill= c(4,2),title = "legend")
					legend("topright",as.character(round(abs(n15.mean-						n15.correct.expect),digit = 4)),title = "Delta")
		
					abline(v = c(n15.mean,n15.correct.expect),col = c(4,2))
					
					#abline(v=n15.correct.expect,col = "red")
				}
				
				graphics.off()
		
				colnames(n15.correct.matrix ) <- colnames(n15.log)
				rownames(n15.correct.matrix ) <- rownames(n15.log)
				pep.all.mean <- n15.correct.matrix
			}
			
			norm.type		<- rep("N15",dim(pep.all.mean)[2])
			pep.all.mean.n	<- pep.all.mean
		}
	}
	
#	pdf("distribution.pdf")
	
#				try(hist(as.numeric(pep.all.mean.n),main = "untransformed data"))
	if(gui.input$n15.log2 & !gui.input$N15){		
		print("performing log2 transformation")	
				#assign("pep.all.mean.",pep.all.mean.n,envir=.GlobalEnv)
	
				pep.all.mean.n <- log2(pep.all.mean.n)
				pep.all.mean.n[pep.all.mean.n == Inf] 	<- 9999 
				pep.all.mean.n[pep.all.mean.n == -Inf] 	<- -9999
				print(n15.correct.method)
				if(n15.correct.method != "none"){
				print("start 15N")
				n15.correct.expect <- log2(n15.correct.expect)
				pdf("log2-correction-label-free.pdf")
				
				n15.correct.value <- c()
				n15.correct.matrix <- c()
				

				for(.n15 in 1: dim(pep.all.mean.n)[2]){
					
					temp.n15 <- pep.all.mean.n[,.n15]
			#assign("temp.n15",temp.n15,env = .GlobalEnv)
			
			
					if(n15.correct.method == "median"){ n15.mean <- median(temp.n15,na.rm = TRUE)}
					if(n15.correct.method == "mean"){ n15.mean <- mean(temp.n15,na.rm = TRUE)}
					if(n15.mean> 0){temp.correct.value <- n15.mean}else{
					temp.correct.value <- n15.mean*-1	}
					temp.n15.correct   <- as.numeric(temp.n15)+	as.numeric(temp.correct.value)
					n15.correct.value <- c(n15.correct.value,temp.correct.value )
					print("binding vectors")
					n15.correct.matrix <- cbind(n15.correct.matrix,temp.n15.correct)
						
						
					#hist(temp.n15,xlim = range(temp.n15, temp.n15.correct,na.rm = TRUE),col = "#00009f90", freq = FALSE)
					plot.input <- c(1,2,3,4)
					#print(temp.n15)
					try(plot.input	<- density(temp.n15,na.rm = TRUE))
					try(plot.input2	<- density(temp.n15.correct,na.rm = TRUE))

					plot(plot.input,xlim = range(temp.n15, temp.n15.correct,na.rm = TRUE),main=paste(colnames(pep.all.mean.n)[.n15],"\nCorrection of log2 data"),type = 		"l",col = "blue",lwd =3)
		
					#hist(temp.n15.correct,add = TRUE,col = "#99000090",freq = FALSE)
					points(plot.input,main="Density estimate of data",type = "l",col = "blue",lwd =3)
					
					points(plot.input2,main="",type = "l",col = "red",lwd =3)

					
		
					legend("topleft",c("uncorrected","corrected"),fill= c(4,2),title = "legend")
					legend("topright",as.character(round(abs(n15.mean-						n15.correct.expect),digit = 4)),title = "Delta")
		
					abline(v = c(n15.mean,n15.correct.expect),col = c(4,2))
					
					#abline(v=n15.correct.expect,col = "red")
				}
				
				graphics.off()
		
				colnames(n15.correct.matrix ) <- colnames(pep.all.mean.n)
				rownames(n15.correct.matrix ) <- rownames(pep.all.mean.n)
				pep.all.mean.n <- n15.correct.matrix
			}
										print("control point")

				
				
				#pdf("log2.pdf")
				#try(hist(as.numeric(pep.all.mean.n),"log2 transformed data"))
				#dev.off()
	}
						print("control point")

	
					
	####
	# group normalisation
	####
	#group.v <- NULL
	if(group.norm == TRUE & exists("exp.set.backup")){
		temp.order.cols <- colnames(pep.all.mean.n)
		temp.order <- c()
		for(i.ord in 1:length(temp.order.cols)){
			grep.i.ord	 <- as.character(temp.order.cols[i.ord])== tolower(as.character(exp.set.backup[,1]))
			if(length(grep.i.ord) == 0){
				temp.order[i.ord] <- "Error, not mapped"
				
			}else{
				
				temp.order[i.ord] <- exp.set.backup[grep.i.ord,3]
				
			}
			
			
			
		}
		
		exp.set.backup
			write.csv(temp.order,"yeah.csv")

	group.v <- 	temp.order
		
	}else{group.v <- NULL}
	



	# rownorm of all peptides
	if(row.norm == TRUE & calc.empai == FALSE){

			pep.all.mean 		<- hz.norm(pep.all.mean.n,1,norm = norm.method,group = group.v)$x
				
		if(length(cbn.prot) > 0 & row.target.norm){
			pep.all.mean			<- t(apply(as.matrix(pep.all.mean),1,function(x){as.numeric(x)/as.numeric(target.mean)}))
		}
	}else{
		pep.all.mean			 <- pep.all.mean.n
	}
	colnames(pep.all.mean.n) <- colnames(pep.all.mean)
	
	write.pep.all.mean.n	 <- pep.all.mean

	#### db exclusion
	if(exists("database")|maxq.exclu == TRUE){
		exclu=TRUE
	} else {
		exclu = FALSE;
		if(re.pep.ex == TRUE){
			ui$messageBox(title="",message="No Database loaded, continuing without redundance exclusion!",icon="error",type="ok")
		}
	}

	####
	#redundant peptide exclusion:
	####

	if(re.pep.ex == TRUE & exclu == TRUE & calc.empai == FALSE){
		#get seq from db:
		sequences 		<- strsplit(as.character(rows),"#")
		sequences.temp 	<- c()
		acc.temp 		<- c()
		acc.seq.temp	<- c()
		for(sequ in 1:length(rows)){
			temp 		<- sequences[[sequ]][2]
			temp.acc 	<- sequences[[sequ]][1]
			sequences.temp 	<- c(sequences.temp,temp)
			acc.temp 		<- c(acc.temp,temp.acc)
			acc.seq.temp	<- c(acc.seq.temp,paste(sequences[[sequ]][1],sequences[[sequ]][2],sep = "#"))
		}
		if(length(sequences.temp) != length(rows)){
			alarm()
		}
		#stop()
		if(maxq.exclu == FALSE){
			if(!exists("pb")){	
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				ratio.prog 			<- prog.max/length(unique(sequences.temp))
			}

			if(1==1){
				test <- 	apply(as.matrix(cbind(unique(sequences.temp),1:length(unique(sequences.temp)))),1,function(x){
					pb.check <- class(try(ui$setProgressBar(pb, as.numeric(x[2])*ratio.prog, label=paste("Pep.Red.Exclu.:", round(as.numeric(x[2])/length(unique(sequences.temp))*100),  "% done"))))
					test 	<- grep(x[1],database[,2],fixed = T)
					grep.r 	<-  as.character(database[test ,1])
					if(length(unique(grep.r)) > 1){
						return.vec <- paste(c(as.character(x[1]),unique(grep.r)),collapse = "|")
						return(return.vec)
					}else{return("")
					}
				})
			
			}

			print(time2-time1)
			test <- test[!test==""]
			if(length(test) != 0){
				write.csv(test,"redundant-peptides.csv")	
				test.split <- strsplit(test,"|",fixed = T)
				seq.temp <- c()
				for(i in 1: length(test)){
					
					temp.i <- test.split[[i]]
					seq.i <- temp.i[1]
					temp.i <- paste(temp.i[-1],collapse = "|")
					seq.temp <- rbind(seq.temp,c(seq.i,temp.i))	
				}
				
				non.redun.ident.all  <- seq.temp[,1]
				
							}else{
				non.redun.ident.all	 <- ""
			}
				#stop()
		}else{
			cat("exclude redundant peptides with maxquant info")
				try.error <- class(try(exclu.vec 	<- grep(import.list$Proteins.sep,temp.my$Proteins,fixed = TRUE)))			
		if(try.error == "try-error"){
		exclu.vec <- 0
		}
			exclu.vec 	<- setdiff(c(1:dim(temp.my)[1]), exclu.vec)
			non.redun.ident.all  <- unique(temp.my$sequence[exclu.vec])
			acc.seq.temp <- sequences.temp
		}

		if(!all(is.na(non.redun.ident.all))|!all(non.redun.ident.all[!is.na(non.redun.ident.all)] == "")){
		#exclude.vec <- grep(paste(non.redun.ident.all,collapse = "|"),sequences.temp)
		exclude.vec <- hz.merge.control(sequences.temp ,non.redun.ident.all)
		exclude.vec <- exclude.vec[!is.na(exclude.vec)] 
		
		pep.all.mean.n <- pep.all.mean.n[-exclude.vec,]
		}

	}#Redunant Pep End
				
	######
	# start of Calculating Proteinlevels
	######		
	sam.sd 			<- c()
	sam.mean 		<- c()
	sam.count 		<- c()
	sam.mean.row 	<- c()
	sam.sd.row 		<- c()
	all.n..col 		<- c()
	.col.mean.row	<- c()
	experiment 		<- unique(exp.set[,2])
	.col.names		<- colnames(pep.all.mean)
	.col.e 			<- NULL
	aov.data 		<- c()
	aov.data.2 		<- c()
	temp.e.mod.mean.m <- c()
	
	if(Raw == TRUE) {
		experiment <- c(1:dim(pep.all.mean)[2])
	}
	total <- length(experiment)
	
	if(any(outlier != "NA" , norm.tog.pep == TRUE) & 1==0){
	for(e in 1:length(experiment)){
		raw.l <- exp.set[exp.set[,2] == experiment[e],]
		if(is.vector(raw.l)) {
			raw.l <- t(as.matrix(raw.l))
		}
		cat(paste("Calculating Proteinlevels of .column",e,"from",length(experiment),"\n"))
		if(Raw == FALSE) {
			temp.t.grep <- c()
				for(t in 1:dim(raw.l)[1]) {
					
					temp.t <- c(1:length(.col.names))[raw.l[t,1]== .col.names ]#grep(raw.l[t,1],.col.names,fixed = TRUE)
					#print(temp.t)
					temp.t.grep <- c(temp.t.grep,temp.t)
				}
				temp.e <- as.matrix(pep.all.mean[,temp.t.grep] )
				
				.col.e <- NULL
				if(length(colnames(temp.e)) == 0){
					.col.e <- .col.names[e]
					}

		
			
			
			} else {
				temp.e <- as.matrix(pep.all.mean[,e])
				colnames(temp.e) <- .col.names[e]
				.col.e			<- .col.names[e]
			}
		
			##========================================================================================================
			if(norm.tog.pep == TRUE & dim(temp.e)[2] > 1 & ratios == FALSE & length(cbn.prot) == 0) {
				cat(" -- normalisation on base of common peptides --")
				temp.na <- temp.e
				temp.na[is.na(temp.na) == FALSE] <- 0
				temp.na[is.na(temp.na)] <- 1
				temp.na <- apply(temp.na,1,function(x){sum(as.numeric(x))})
				temp.ne <- temp.e[temp.na == 0,]
				
				if(is.vector(temp.ne)) {
					temp.na <- t(as.matrix(temp.ne))
					rownames(temp.na) <- rownames(temp.e)[temp.na == 0]
				} else {
					temp.na <- temp.ne
				}
				
				if(length(temp.na) != 0) {
					temp.mean.tog <- apply(temp.na,1,function(x) {
								sum(as.numeric(x),na.rm = TRUE)
					})
					
					temp.na.temp <- temp.na[temp.mean.tog != 0,]
					if(is.vector(temp.na.temp)) {
						temp.na.temp <- t(as.matrix(temp.na.temp))
					}
					
					tog.sum.all  <- apply(temp.na.temp,2,function(x){
								sum(as.numeric(x),na.rm = TRUE)
					})
					
					
					
					
					if(merge.method == "mean"){
						tog.mean.all <- mean(tog.sum.all,na.rm = TRUE)
					}
					if(merge.method == "median"){
						tog.mean.all <- median(tog.sum.all,na.rm = TRUE)
					}
					
					
					# 3. creating factor:
					factor.vec <- c()
					for(z in 1:dim(temp.e)[2]) {
						factor.vec[z] <- tog.mean.all/tog.sum.all[z]
					}
					
					# -- treat matrix with factor --
					norm.data <- c()
					for(u in 1:dim(temp.e)[2]){
						temp.u 	<- as.numeric(temp.e[,u])*factor.vec[u]
						norm.data <- cbind(norm.data,temp.u)
					}
					
					test.sd		<- apply(norm.data,1,function(x){
																				x <- mean(as.numeric(x),na.rm = TRUE);
								x.sd<- sd(as.numeric(x),na.rm = TRUE)/x;
							return(x.sd/x)
					})
					temp.e.sd	<- apply(pep.all.mean.n ,1,function(x){
																				x <- mean(as.numeric(x),na.rm = TRUE);
								x.sd<- sd(as.numeric(x),na.rm = TRUE)/x;
								return(x.sd/x)
					})
					new.sd <- sum(as.numeric(test.sd),na.rm = TRUE) #/sum(norm.data,na.rm = TRUE)
					old.sd <- sum(as.numeric(temp.e.sd),na.rm = TRUE) #/sum(as.numeric(temp.e),na.rm = TRUE)
					new.sd.r <- sum(as.numeric(new.sd),na.rm = TRUE)
					old.sd.r <- sum(as.numeric(old.sd),na.rm = TRUE)
					#text.out <- paste(#"Normalisation of peptide matrix:\n \n Sd changed in comparison to raw files from 100 % (",round(old.sd),") to", round(new.sd/old.sd* 100),"% (",round(new.sd),").\n",
					#		"\n \nSd changed in comparison to normalisation on sum of peptides from 100 % (",round(old.sd.r),") to", round(new.sd.r/old.sd.r* 100),"% (",round(new.sd.r),").\n")
					cat(paste(rep("*",nchar(text.out)/3),.collapse = "",sep = ""),"\n")
					#cat(text.out,"\n")
					cat(paste(rep("*",nchar(text.out)/3),.collapse = "",sep = ""),"\n")
					colnames(norm.data) <- colnames(temp.e)
					rownames(norm.data) <- rownames(temp.e)
					temp.e <- norm.data
				} else {
					alarm()
					text.out <- paste("ALERT: There are no common peptides available, \nswitching to relative normalisation (sum of peptide intensities)\n")
					cat(paste(rep("*",nchar(text.out)/1.65),.collapse = "",sep = ""),"\n")
				}
			}
			##======================================================================================================
			######### code
			all.temp.o 	<- c()
			all.mean.o 	<- c()
			all.sd.o 	<- c()
			all.n 			<- c()

			row.names.temp.i <- c()
			total3 			 <- length(unique(temp.my$code))
				
			temp.mean.row <- c()



			for(i in 1:length(unique(temp.my$code))) {
				if(i == 1) {
				
					########## GUI 
					pb.check <- class(try(ui$setProgressBar(pb, 0, label=paste( "Calculation: ","0% done"))))
					while(pb.check == "try-error"){
					print("Warning: User closed window!")
					pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
					pb.check <- class(try(ui$setProgressBar(pb, 0, label=paste( "Calculation: ","0% done"))))
					}
					##############
		
				}

				ratio.prog		<- prog.max/length(unique(temp.my$code))
				temp.grep		<- grep(unique(temp.my$code)[i],rownames(temp.e),fixed = TRUE)
				temp.i			<- temp.e[temp.grep,]

				if(is.vector(temp.i)) {
					temp.i 		<- as.matrix(temp.i)
					if(dim(temp.i)[1] > length(temp.grep)) {
						temp.i 	<- t(temp.i)
					}
					rownames(temp.i) <- rownames(temp.e)[temp.grep]
				}
			
				# -- outlier --
				if(outlier == "row"& row.norm == FALSE) {
					for(o in 1: dim(temp.i)[1]) {
						temp.o <- temp.i[o,]
						box.temp.o	<- boxplot.stats(temp.o)$out
						for(r in 1:length(box.temp.o)) {
							temp.r <- box.temp.o[r]
							temp.o[temp.o == temp.r] <- NA
						}
						if(length(box.temp.o) != 0) {
							box.ex <-1
						} else {
							box.ex <- NA
						}
						all.temp.o <- rbind(all.temp.o,temp.o)
					}
					temp.i <- all.temp.o
				}
				
				if(outlier == "row"& row.norm == TRUE) {
											
						box.temp.o	<- boxplot.stats(temp.i)$out
						for(r in 1:length(box.temp.o)) {
							temp.r <- box.temp.o[r]
							temp.i[temp.i == temp.r] <- NA
						}
						
						if(length(box.temp.o) != 0) {
							box.ex <-1
						} else {
							box.ex <- NA
						}
								
				}
				
				if(outlier == "all.below" &  row.norm == FALSE) {
					temp.o		<- temp.i
					box.temp.o	<- boxplot.stats(as.numeric(temp.o))$out
					box.temp.o	<- box.temp.o[box.temp.o < median(as.numeric(temp.i),na.rm = TRUE)]
					for(r in 1:length(box.temp.o)) {
						temp.r <- box.temp.o[r]
						temp.o[temp.o == temp.r] <- NA
					}
					if(length(box.temp.o) != 0){
						box.ex <- box.ex + 1
					} else {
						box.ex <- NA
					}
					temp.i <- temp.o
				}
				
				if(outlier == "top.3"){
				temp.i.test 	<- temp.i[order(temp.i,decreasing = TRUE)]
				temp.i.test  	<- temp.i.test[!is.na(temp.i.test)]
				temp.i.test 	<- temp.i.test[1:3]

				temp.i[setdiff(1:length(temp.i),grep(paste(temp.i.test,collapse = "|"),temp.i))] <- NA

				}
				
				
				# -- merge --
				row.names.temp.i <- c(row.names.temp.i,rownames(temp.i))
				temp.aov  <- temp.i[!is.na(temp.i)]
				temp.aov2 <- rep(as.character(unique(temp.my$code)[i]),length(temp.aov))
				temp.aov3 <- rep(experiment[e],length(temp.aov))
				temp.aov4 <- cbind(temp.aov2,temp.aov3,temp.aov)
				aov.data  <- rbind(aov.data,temp.aov4)
				
				if(merge.method == "median") {
					temp.mean <- median(as.numeric(temp.i),na.rm = TRUE)
				}
				if(merge.method == "mean") {
					temp.mean <- mean(as.numeric(temp.i),na.rm = TRUE) 
				}
					
				
				if(row.norm	== TRUE) {
					
					temp.mean.row 	<- rbind(temp.mean.row,temp.i)
					if(Raw == TRUE | length(.col.e) != 0){colnames(temp.mean.row) <- .col.e}
					
					
					
				} else {
					temp.sd		<- sd(as.numeric(temp.i),na.rm=TRUE)
					temp.sd 	<- temp.sd/mean(as.numeric(temp.i),na.rm=TRUE)
					temp.n		<- length(as.numeric(temp.i)[!is.na(as.numeric(temp.i))])
					all.n		<- c(all.n,temp.n)
					all.mean.o 	<- c(all.mean.o,temp.mean)
					all.sd.o 	<- c(all.sd.o,temp.sd)
					
					
				}

				if(i%%10==0){
					
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Calculation: ", ".col",e,"/",length(experiment),"|",round(i/total3*100, 0),  "% done"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Calculation: ", ".col",e,"/",length(experiment),"|",round(i/total3*100, 0),  "% done"))))
		}
	##############
		
					
				}

			}
			
			
			
#print(dim(all.mean.o))
		if(row.norm == FALSE){	
			
			
			#print(dim(sam.mean))
			#print(dim(sam.sd))
			
			sam.mean 	<- cbind(sam.mean,all.mean.o)
			sam.sd		<- cbind(sam.sd,all.sd.o) 
			all.n..col 	<- cbind(all.n..col,all.n)
		}else{
	
	
			temp..col.mean.row			<- colnames(temp.mean.row)
			names(temp..col.mean.row)	<- rep(e,length(temp..col.mean.row))
			
			.col.mean.row 		<- c(.col.mean.row, temp..col.mean.row)
			temp.mean.row		<- cbind(as.character(rownames(temp.mean.row)),temp.mean.row)
			temp.mean.row		<- unique(temp.mean.row)
			#print(dim(temp.mean.row))
			if(e==1){
				
				
				sam.mean 	<- merge(as.matrix(unique(rownames(pep.all.mean.n))),temp.mean.row,all = TRUE)}else{
					
					
		########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Merge data, Rows:",dim(sam.mean)[1],"..."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, i*ratio.prog, label=paste("Merge data, Rows:",dim(sam.mean)[1],"..."))))
		}
	##############
						
				
				sam.mean	<- merge(as.matrix(sam.mean),as.matrix(temp.mean.row),by = 1,all = TRUE )
			
			
			}
		}
				
	}
	}else{
		print("skip calculation")
	sam.mean <- cbind(rownames(pep.all.mean), pep.all.mean)	
	}	# skip calculation 
	
	#if(all(c(outlier == "NA"|norm.tog.pep == FALSE,row.norm == FALSE))){
	#	rownames(sam.mean) <- unique(temp.my$code)		
	#}
	
	#stop()

	print("Reached Control point!")

		if(all(outlier == "NA" & norm.tog.pep == FALSE) | 1==1){
			if(length(dupli.pep <-grep(TRUE,duplicated(sam.mean[,1]))) > 0){
				print("warning, found duplicated peptides.")
				sam.dupli <- sam.mean[dupli.pep,]
				sam.mean <- sam.mean[-dupli.pep,]
				
				#sam.dupli[,1] <- paste(sam.dupli[,1],"# unexpected.duplicated.peptide",c(1:length(sam.dupli[,1])))
				
			} 
			
			
			rownames(sam.mean) 	<- sam.mean[,1]
			sam.mean 			<- sam.mean[,-1]
		} else {
			sam.mean 			<- as.matrix(sam.mean)
			if(length(rownames(sam.mean) ==0)){
#				rownames(sam.mean) 	<- unique(temp.my$code)
			}
		}
		
		backup					<- sam.mean
		rows 					<- rownames(sam.mean)
		sam.mean 				<- as.matrix(apply(sam.mean,2,as.numeric))
		rownames(sam.mean) 		<- rows

sam.mean.phospho <- matrix()
			if(phospho){
				sam.mean.backup <- sam.mean
				.phospho 			<- grep(gui.input$phospho.string ,rownames(sam.mean),fixed = TRUE)
				if(length(.phospho)> 0){
				sam.mean.phospho	<- sam.mean[.phospho,]
				sam.mean			<- sam.mean[-.phospho,]
				
				if(exists("sam.sd")){
				sam.mean.phospho.sd	<- sam.sd[.phospho,]
				sam.mean.sd			<- sam.sd[-.phospho,]
				}
				all.n.phospho.col <- all.n..col[.phospho,]
				all.n..col <- all.n..col[-.phospho,]
				
				
				}
				
			}

		if(script.shape){
			#input 				<- sam.mean
			#source(paste(path1,"scripts/script.shape.R",sep = ""))
			#input 				<- script.shape.fun(input)
			#sam.mean.shape 		<- input
			#sam.mean		 	<- sam.mean.shape$shape
		}else{
					
			if(group.filter & gui.input$exp.design != ""){
				group.filter.order <- hz.merge.control(tolower(make.names(exp.group.shape[,1],allow = F)),colnames(sam.mean))
				group.shape <- exp.group.shape[group.filter.order,3]
		
			}else{
				group.shape <- rep(1,dim(sam.mean)[2])
			}
			sam.mean.shape 		<- hz.shape(as.matrix(sam.mean),shape = shape, group.shape)
			shape.vec			<- hz.merge.control(rownames(sam.mean),intersect(rownames(sam.mean),rownames(sam.mean.shape$shape)))
			shape.vec 			<- shape.vec[!is.na(shape.vec)]
						
			sam.mean		 	<- sam.mean.shape$shape
			all.n..col 			<- all.n..col[sam.mean.shape$n.vec,]

# Sam.sd kontrollieren
# all.n..col kontrollieren

			if(phospho & dim(sam.mean.phospho)[1] > 0){
				sam.mean 	<- rbind(sam.mean.phospho,sam.mean)
				
			try.error<- class(try(sam.mean.sd <- rbind(sam.mean.phospho.sd, sam.sd)))			
			if(try.error == "try-error"){
				try(sam.mean.sd <- rbind(matrix(0,nrow = dim(sam.mean.phospho)[1],ncol = dim(sam.mean.phospho)[2])))
			}
			
				all.n..col  <- rbind(all.n.phospho.col,all.n..col)
			}
			
			#if(any(outlier != "NA" ,norm.tog.pep == FALSE) & row.norm == FALSE){
			#	sam.sd <- sam.sd[sam.mean.shape$n.vec,]
			#}
			
			

		}
		
		if(all(outlier == "NA" & norm.tog.pep == FALSE) | 1==1){
			sam.sd		 		<- sam.sd[sam.mean.shape$n.vec,]
		} 
		
		if(is.vector(sam.mean)){
			sam.mean <- as.matrix(sam.mean)
		}

		choosen.peptides 	<- rownames(sam.mean)
		if(length(sam.mean) == 0 | dim(sam.mean)[1] == 0) {
			ui$messageBox(title="",message="No Peptides left after Shaping!\nPlease re-run cRacker with reduced shape value\nor include duplicate peptides!",icon="error",type="ok")
		}
		
		if(dim(sam.mean)[2]< 2) {
			ui$messageBox(title="An error has occured!",message="Row-Normalisation of IonIntensities of 1 .column is not meaningful!",icon="error",type="ok")
		}
		
		temp.row <- strsplit(choosen.peptides,"#",fixed = TRUE)	
		choosen.proteins <- c()
		for(z in 1 : length(choosen.peptides)) {
			temp 				<- temp.row[[z]][1]
			choosen.proteins 	<- c(choosen.proteins,temp)
		}

		if(Raw == FALSE){
			names..col.init <- c()
			for(re.i in 1:length(colnames(sam.mean))){
				temp.re.i<- colnames(sam.mean)[re.i]
				test.re.i 	<- exp.set[temp.re.i ==exp.set[,1],2]
				names..col.init  <- c(names..col.init,test.re.i) 
				
			}
			
		}else{
		names..col.init <- colnames(sam.mean)	
		}
#stop()
		rownames(sam.mean) <- gsub("  "," ",rownames(sam.mean),fixed = TRUE)

		if(all(outlier == "NA" & norm.tog.pep == FALSE) | 1==1) {	
			names..col 	<- names..col.init
			data.mean 	<- c()
			data.sd		<- c()
			data.sd.rel <- c()
			data.n		<- c()
			data.top3.rpn <- c()
			aov.data.2 <- c()
			acc.temp 	<- strsplit(rownames(sam.mean),"#", fixed = TRUE)
			acc 		<- c()

			for( i in 1:dim(sam.mean)[1]){
				acc <- c(acc,acc.temp[[i]][1])
			}
			uni.acc <- 	unique(acc)
			sys1 <- Sys.time()
			
			ratio.prog <- prog.max/(length(unique(names..col))*length(uni.acc))
			assign("cracker.counter.temp",1,envir = .GlobalEnv)
			assign("cracker.ratio.prog.temp",ratio.prog,envir = .GlobalEnv)

			
			for( i in 1:length(unique(names..col))){
				temp 		<- as.matrix(sam.mean[,names..col == unique(names..col)[i]])
				
				test		<- as.matrix(aggregate(as.vector(temp),list(rep(acc,dim(temp)[2])),function(x){x <- hz.agg.fun(x,outlier,row.norm,merge.method,ui,pb,prog.max,c(i,length(unique(names..col))))}))
				test.order  <- hz.merge.control(test[,1],uni.acc)


				temp.split <- strsplit(test[,8],"#",fixed = T)
				aov.data <- c()
				for(s in 1:dim(test)[1]){
						temp.s <- as.numeric(temp.split[[s]])
						temp.s <- temp.s[!is.na(temp.s)]
						#print(s)
						if(length(temp.s) > 0){
						temp.s <- cbind(test[s,1],unique(names..col)[i],as.matrix(temp.s))
						aov.data <- rbind(aov.data,temp.s)
						}
				}
#print(i)

				data.mean 	<- cbind(data.mean,test[test.order,3])
				data.sd 	<- cbind(data.sd,test[test.order,4])
				data.sd.rel <- cbind(data.sd.rel,test[test.order,5])
				all.n..col	<- cbind(all.n..col,test[test.order,6])
				data.top3.rpn 	<- cbind(data.top3.rpn,test[test.order,7])
				aov.data.2 <- rbind(aov.data.2,aov.data)
				
			ratio.prog 	<- prog.max/length(unique(names..col))
								
			########## GUI 
				pb.check <- class(try(ui$setProgressBar(pb,i*ratio.prog, label=paste(i,"\\",length(unique(names..col)),"Peptide-Mean",round(as.numeric(i)/length(unique(acc))*100, 0),"% done"))))
	
				while(pb.check == "try-error"){
					print("Warning: User closed window!")
					pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
					pb.check <- class(try(ui$setProgressBar(pb,as.numeric(i)*ratio.prog, label=paste(i,"\\",length(unique(names..col)),"Peptide-Mean",round(as.numeric(x[2])/length(unique(acc))*100, 0),"% done"))))
				}
			##############
				
				
				
			}
		colnames(aov.data.2) <- c("accession","experiment","intensity")
		all.n..col <- as.matrix(all.n..col)
		sam.sd 		<- data.sd.rel
		sam.mean 	<- data.mean
	}
	
	sys2 <- Sys.time()
	print("hurray")
	# -- give names --

	if(Raw == TRUE) {
		colnames(sam.mean) 	<- unique(.col.all)
		rownames(sam.mean) 	<- unique(choosen.proteins)
		sam.sd 				<- as.matrix(sam.sd)
		colnames(sam.sd) 	<- paste("SD",unique(.col.all))
		rownames(sam.sd) 	<- unique(choosen.proteins)
			
	}
		
	if(Raw == FALSE){
		colnames(sam.mean) 	<- paste(unique(temp.my$sam_id),unique(temp.my$sam_name))
		rownames(sam.mean) 	<- unique(choosen.proteins)
		sam.sd 				<- as.matrix(sam.sd)
		colnames(sam.sd) 	<- paste("SD" ,unique(temp.my$sam_id),unique(temp.my$sam_name))
		rownames(sam.sd) 	<- unique(choosen.proteins)

	
	}	
	all.n..col 				<- as.matrix(all.n..col)
	colnames(all.n..col) 	<- colnames(sam.mean)
	rownames(all.n..col) 	<- rownames(sam.mean)
	if(score > 0){
		temp.my <- temp.my[as.numeric(temp.my$Score) > score,]
	}
	
	###### Phospho - ratios
	


	###### info - data		
	info.data <- 	as.data.frame(cbind(
			as.character(temp.my$code),
			as.character(temp.my$Description),
			temp.my$Score,
			temp.my$mcr,
			temp.my$charge,
			temp.my$Calibrated.mass.relative.error..ppm,
			as.character(temp.my$sequence)
			))
	.colnames.info.data <- as.character(c(
			"code",
			if(length(temp.my$Description)>0){ "Description" },
			if(length(temp.my$Score)>0){ "Score" },
			if(length(temp.my$mcr)>0){ "mcr" },
			if(length(temp.my$charge)>0){ "charge" },
			if(length(temp.my$Calibrated.mass.relative.error..ppm)>0){ "Calibrated.mass.relative.error..ppm" },
			if(length(temp.my$sequence)>0){ "sequence" }
		))
		
	all.peptides <- as.matrix((paste(
			temp.my$code,
			temp.my$sequence,
			round(as.numeric(temp.my$mcr),digits = 0),
			temp.my$charge,sep = "#"
		)))
		
	if(all(outlier == "NA" & norm.tog.pep == FALSE) | 1==1) {
		all.peptides			<-  as.matrix((paste(
									temp.my$code,
									temp.my$sequence,
									round(as.numeric(temp.my$mcr),digits = 0),
									temp.my$charge,sep = "#"
			)))
		all.peptides.merge 		<- cbind(all.peptides,c(1:dim(all.peptides)[1]))
		all.peptides.merge		<- merge(all.peptides.merge,choosen.peptides,by = 1 )
	}
		
	info.data <- cbind(all.peptides,info.data)
	colnames(info.data) <- c("id",.colnames.info.data)
	dec.vec <- c()
	if(add.data == TRUE) {
		
	info.data$Score	<- as.numeric(as.character(info.data$Score))
	info.data$Score[is.na(info.data$Score)]	<- 0	
		info.data.protein 	<- c()
		total 				<-  length(rownames(sam.mean))
			info 				<- function(x){
								if(as.numeric(x[1]) %%10==0 & exists("pb")) {
									
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb,as.numeric(x[1])*ratio.prog, label=paste("Extracting info: ",round(as.numeric(x[1])/total*100, 0),"% done"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb,as.numeric(x[1])*ratio.prog, label=paste("Extracting info: ",round(as.numeric(x[1])/total*100, 0),"% done"))))
		}
	##############	
									
								}
								
								x 		<- as.character(x[2])
								temp	<- grep(x,as.character(info.data[,1]),fixed = TRUE)
								temp.my.z		<- info.data[temp,]
								
								temp.score		<- as.numeric(as.vector(temp.my.z$Score[]))
								temp.my.z 		<- temp.my.z[as.numeric(as.vector(temp.score)) == max(as.numeric(as.vector(temp.score)),na.rm = TRUE),]
								temp.my.z  		<- temp.my.z[!is.na(temp.my.z$id),]
								
								if(dim(temp.my.z)[1] > 1 ){
									mass.accuracy	<- as.vector(temp.my.z$Calibrated.mass.relative.error..ppm)
									mass.accuracy <- as.numeric(as.vector(mass.accuracy)) == min(as.numeric(as.vector(mass.accuracy)),na.rm = TRUE)
									mass.accuracy[is.na(mass.accuracy)] <- FALSE
									
									
									temp.my.z 		<- temp.my.z[mass.accuracy,]
									temp.my.z		<- as.vector(as.matrix(temp.my.z))
								}
								#if(is.vector(temp.my.z) == FALSE){print(x)}
								if(length(temp.my.z) != 8){alarm();print(x);temp.my.z <- rep("error in data collection!",8)}
								
								return(as.vector(as.matrix(temp.my.z)))
							}
			
		info.temp 			<- cbind(c(1:dim(sam.mean)[1]),as.matrix(rownames(sam.mean)))
		ratio.prog 			<- prog.max/dim(info.temp)[1]
		info.data.protein 	<- t(apply(info.temp ,1,info))
			
			
		try(colnames(info.data.protein ) <- colnames(info.data)		)
		

		
			########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, prog.max, label=paste("Extracting info: ","100 % done"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, prog.max, label=paste("Extracting info: ","100 % done"))))
		}
	##############
		
		
		
		
				

					
		if(length(info.data.protein) == 0){					info.data.protein 		<- rep("no.data..collected",dim(sam.mean)[1])
		}		
					
		if(dim(info.data.protein)[1] != dim(sam.mean)[1] ){
			info.data.protein 			<- rep("mistake in data .collection",dim(sam.mean)[1])
		}	
	}else{ 
		info.data.protein 			<- rep("no.data..collected",dim(sam.mean)[1])
	}
	info.data.matrix <- cbind(info.data.protein,sam.mean)
	# -- Finished Read additional data --
	#################################################################sam.mean
	# -- SD of proteinmatrix --
	
	if(outlier!= "NA") {
		if(is.na(box.ex)) {
			print(paste("WARNING: no outliers detected!"))
		} else {
			print(paste("Outliers detected for",box.ex,"proteins!"))
		}
	}
	cat("\n")
	cat(	"****************************************************************************************\n")
	cat(paste("************** FINISHED matrix.creation ************************************************\n"))
	cat("****************************************************************************************\n")
	cat("\n")

	#stop()
if(all(outlier == "NA" & norm.tog.pep == FALSE) | 1==1){
		if(length(aov.data.2) >0){
			colnames(aov.data.2) = c("accession","experiment","intensity")
		}
		if(length(aov.data) >0){
			colnames(aov.data) = c("accession","experiment","intensity")
		}
	aov.export.1 = aov.data.2
	aov.export.2 = aov.data	
	}else{
	colnames(aov.data) = c("accession","experiment","intensity")

	aov.export.1 = aov.data
	aov.export.2 = "no data"
	}


	}
	
		used.peptides <- c()
		if(group.filter & gui.input$exp.design != ""){
				group.filter.order <- hz.merge.control(tolower(make.names(exp.group.shape[,1],allow = F)),colnames(write.pep.all.mean.n))
				group.shape <- exp.group.shape[group.filter.order,3]
		
			}else{
				group.shape <- rep(1,dim(sam.mean)[2])
			}
		
			
		
		try(used.peptides <- hz.shape(write.pep.all.mean.n,shape = shape,group.shape)$shape)

	if(!exists("info.data.matrix")){info.data.matrix <- "not collected"}
	if(!exists(data.top3.rpn)){		data.top3.rpn <- NULL
	}
	return(list(	x=sam.mean, 
					x.sd =sam.sd,
					peptidelist=write.pep.all.mean.n,
					proteinlist.info = info.data.matrix,
					used.peptides = used.peptides,
					method = paste("type: ",type,"| Raw:",as.character(Raw),"| merge.method",merge.method,"| outlier",outlier),
					prot.norm = cbn.prot.data,
					prot.n = all.n..col,
					aov.export.1 = aov.export.1, 
					aov.export.2 = aov.export.2, 
					mod.peptides.experiment = temp.e.mod.mean.m,
					exp.design 	= exp.set,
					rpn 		= rpn.start.peptide.matrix 
					))
}
}
