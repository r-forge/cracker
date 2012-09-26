hz.script <-
function(path1= NA , path2.set = list("NA","maxquant","default") , import.list=NULL,.data = NA){
require("tcltk2")
tk2font.set("TkDefaultFont",settings= "-family Tahoma -size 10 -weight normal")   



path2 			<-	normalizePath(path2.set$path)
	ratio.prog <- 10000
path1 			<- normalizePath(path1)

path2.test 		<- class(try(setwd(path2)))
if(path2.test == "try-error"){
 	path2.input.file	<- basename(path2)
	path2				<- dirname(path2)
	path2.set$path 		 <-	path2
	path2.set$input.file <- path2.input.file	
	
	
}else{
	path2.input.file <- ""
}	



if(is.na(path1)){
	try(path1 <- paste(path.package("cRacker"),"data",sep = "/"))
}
	if(!exists("ratio.prog")){ratio.prog <- 1000}

	
data(cracker.ui.tk)
try(require(gplots))
path2 			<-	normalizePath(path2.set$path)
print("starting Script")

print(getwd())
pb.check <- "numeric"

path.check 		 	<- list.files(path2)
#print(test)
path.check 			<- grep("matrix-binary.Rdata",path.check)
if(length(path.check) == 0){build.matrix <- "1"}else{build.matrix <- "0"}


ui <- cracker.ui.tk;	# Use TclTk for ui
ui$init();				# Init (loads library)
prog.max <- 10000
pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)


import.list <- import.list[import.list$file.type ==path2.set$engine,]

if(!exists(".data")|!is.data.frame(.data)){
print("Loading .data")

assign("import.list",list(import.list = import.list, path.data = path2.set$path, path2.input.file = path2.input.file,prog.max=prog.max,ui=ui,pb=pb),envir = .GlobalEnv)

if(path2.set$data!= "default"){
	assign("path2.set.test",path2.set,envir = .GlobalEnv)
	load(path2.set$data)
	
}else{
	
	try(	.data 		<- hz.import(import.list = import.list, path.data = path2.set$path, path2.input.file = path2.input.file,prog.max=prog.max,ui=ui,pb=pb))
try(print(dim(.data)))

}

try(print(dim(.data)))
#assign(".data",.data,envir = .GlobalEnv)

#assign("import.list", import.list, envir=globalenv())
#assign("path2.set", path2.set, envir=globalenv())
#assign("path2.input.file", path2.input.file, envir=globalenv())
#assign(".data.a", .data, envir=globalenv())
#if(all())
try(.data$code <- make.names(.data$code) )
#print(path2.set$path)
}

if(exists("pb")){
close(pb)
rm(pb)
}
if(exists(".data")){
	if(all(is.na(.data))){
			ui$messageBox(title="Abort",message="Import of data failed!\nPlease check format of your file and rerun cRacker!\nCurrent cRacker session  can be only used for importing new libraries.",icon="error",type="ok") 
	}
}else{
	.data <- NA
}

try(rm(ui))
#assign(".data.test", .data, envir=globalenv())


gui.input 		<- hz.read.parameters(image.path = NULL, build.matrix = build.matrix,path2 = path2, path2.set = path2.set,.data=.data,path1=path1)





	if(all(is.na(.data))){
			ui$messageBox(title="Abort",message="No data loaded, please restart cRacker!",icon="error",type="ok") 
}
	
loop.control <- 1
if(length(gui.input)==1){
	if(gui.input == "reload"){loop.control <- 2}
}

if(length(gui.input)==1){
while(gui.input== "switch"|gui.input == "reload"){
	
	if(loop.control == 1 & length(gui.input)==1){
		loop.control 	<- 2
		print("hu")
		gui.input 		<- hz.read.empai.parameters(image.path = NULL, build.matrix = build.matrix,path2 = path2, path2.set = path2.set,.data=.data,path1 = path1)
	if(gui.input == "reload"){loop.control <- 1}

	}
	
if(loop.control == 2 & length(gui.input) == 1){
	if(gui.input != "stopped"){
	print("re")
		gui.input 		<- hz.read.parameters(image.path = NULL, build.matrix = build.matrix,path2 = path2, path2.set = path2.set,.data=.data,path1 = path1)
			loop.control 	<- 1
			if(gui.input == "reload"){loop.control <- 2}

		}
	}
}




if( gui.input == "stopped"){
	ui <- cracker.ui.tk
	ui$messageBox(icon="warning",message="Abort by user!")
	}
}else{
	if(substr(gui.input$cracker,(nchar(gui.input$cracker)),nchar(gui.input$cracker)) == "/"){
	
	gui.input$cracker <- substr(gui.input$cracker,1,nchar(gui.input$cracker)-1)
}
}	
if(gui.input$N15){
	if(!exists(".cRacker.check.N15.loaded")|1==1){
	#assign(".cRacker.check.N15.loaded", TRUE, envir=.GlobalEnv)
		 
print("Loading .data")
	.data 		<- hz.import(import.list = import.list, path.data = path2.set$path, path2.input.file = path2.input.file,prog.max=prog.max,ui=ui,pb=pb,N15=TRUE)
		print(dim(.data))
print(path2.set$path)
}

}
#.data <- .data[!is.na(.data$code),]


path2 			<-	normalizePath(path2.set$path)
	ratio.prog <- 10000
path1 <- normalizePath(path1)

path2.test 		<- class(try(setwd(path2)))
if(path2.test == "try-error"){
 	path2.input.file	<- basename(path2)
	path2.set$path		<- dirname(path2)
path2 <- dirname(path2)	
}else{
	path2.input.file <- ""
}

####
## exclude samples
#####
if(gui.input$exp.design!=""){
try(exp.design.temp <- read.table(gui.input$exp.design, stringsAsFactors = FALSE,header = T))


if(any(unique(exp.design.temp$Include)>1)){
	print(text.warning<- "error in experimental design, exclusion string is not binary!")
	tkmessageBox(title="Message",message=text.warning,icon="warning",type="ok")

	
	exp.design.temp$Include[exp.design.temp$Include!= 0] <- 1

	write.table(exp.design.temp,gui.input$exp.design,sep = "\t")

}
}

if(exists("exp.design.temp")){
	exclude.string <- exp.design.temp$Name[exp.design.temp$Include == 0 ]
	if(length(exclude.string)> 0){
		for(i in 1:length(exclude.string)){
			print(paste("Excluding samples:",exclude.string[i]))
			print(dim(.data))

			.data <- .data[!make.names(.data$rawfilename,allow = F) == exclude.string[i],]
		}
	}
	
}


gui.input$shape				= (100-gui.input$shape)/100
gui.input$shape.prot.norm	= gui.input$shape.prot.norm/100#


if(gui.input$time.grouped & gui.input$exp.design== ""){
	gui.input$time.grouped <- FALSE
}


if(gui.input$norm.method == "z-score"){
	gui.input$norm.method <- "z-pos"	
}


if(gui.input$color.plots == "rainbow"){
	gui.input$color.plots 	= TRUE
	colorblind.set			= FALSE
}
if(gui.input$color.plots == "colorblind"){
	gui.input$color.plots 	= TRUE
	colorblind.set			= TRUE
}
if(gui.input$color.plots == "greytone"){
	gui.input$color.plots 	= FALSE
	colorblind.set			= FALSE
}


calc.empai 	= gui.input$calc.empai
empai.sd 	= gui.input$empai.sd

calc.empai.list = gui.input$calc.empai.list


gui.backup 		<- gui.input

if(gui.input$calc.empai == "TRUE"){
	

	
	print("loading empai reference data, this can take up to a minute")
#.length.matrix <- read.csv(paste(path1,"empai/",,sep = ""))
print(gui.input$empai.reference)
	load(paste(gui.input$cracker,"/cRackerEmPAI-", gui.input$empai.reference,sep = ""))
	.length.matrix <- .length.matrix
	object <- .length.matrix

}else{
	.length.matrix <- NULL
}	

error.try <- class(.error<- try(hz.script.y.lab.return <- hz.script.y.lab(gui.input = gui.input)))

if(error.try == "try-error"){
		tkmessageBox(title="Message",message=paste("Error in ylab string initiationn!\n",.error),icon="warning",type="ok")
	

}	
if(gui.input$exclu == TRUE){

	.wd <- getwd()
	setwd(gui.input$cracker)

	if(nchar(gui.input$db) !=0){	
		load(paste(gui.input$db,sep = ""))
	}else{
print("Database could not be loaded!")

	}
setwd(.wd)	
}
#progress
ui <- cracker.ui.tk;	# Use TclTk for ui
ui$init();				# Init (loads library)
prog.max <- 10000
pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)


#progress


if(gui.input$plot.only != "" | gui.input$plot.only == FALSE ){
	
		##############	GUI
	.label <- "Loading results from previous analysis"
	pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))

	while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))
	}
	##############
	
	try.out <- try(load(gui.input$plot.only))
	
	path.data 	<- try(dirname(gui.input$plot.only))
	foldername 	<- basename(path.data)
	path.dec <- function(path){
	wd <- getwd()
	setwd(path)
	setwd("../")
	path <- getwd()
	setwd(wd)
	return(path)
	}
	gui.input$path.data <- path.dec(path.data)
	
	
	if(try.out == "try-error"){
		print("could not find binary.Rdata\nPlease process data first!")
		tkmessageBox(title="Message",message="Rdata could not be loaded!\nApplication stopped.",icon="warning",type="ok")
		stop()
	}
	
	if(!exists(".data2")){
		print("The loaded Rdata file contains no accessible objects")
		tkmessageBox(title="Message",message="Rdata could not be loaded!\nApplication stopped.",icon="warning",type="ok")
		stop()
	}
	
	gui.input$raw 		<- .data2$gui.input$raw 	
	gui.input$empai.sd 	<- .data2$gui.input$empai.sd 	
}




#####################
# Folder Creation
setwd(normalizePath(gui.input$path.data))

	##############	GUI
	.label <- "creating analysis folder"
	pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))

	while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))
	}
	##############

print("start creating folder")
wd 	<- getwd()
.add <- ""
if(gui.input$N15){.add <- paste("N15",.add,sep = "-")}
if(length(gui.input$prot.norm) >0){.add <- paste(.add,"Ref-Prot",sep = "-")}
if(gui.input$calc.empai& empai.sd){.add <- paste(.add,"exp",sep = "")}
if(gui.input$calc.empai& !empai.sd){.add <- paste(.add,"raw",sep = "")}
if(gui.input$raw == FALSE & gui.input$calc.empai == FALSE){.add <- paste(.add,"exp",sep = "")}else{.add <- paste(.add,"raw",sep = "")}
if(gui.input$calc.empai){.add <- paste(.add,"empai",sep = "-")}else{.add <- paste(.add,"IonIntensity",sep = "-")}


if(gui.input$plot.only ==  ""){
foldername 	<- paste(gui.input$expname,.add,Sys.Date(), sep="-")
}

print(foldername)
dir.create(.setpath <- paste(path2,foldername,sep = "/")) 
print(.setpath)
setwd(.setpath)

wd.write 	<- getwd()
#######
print("Writing parameters!")

parameters.write <- function(){
zz <- file("parameters.txt", open="wt")
sink(zz)
sink(zz, type="message")
print(gui.input)
sink(type = "message")
sink()
settings <- gui.input$settings
save(settings,gui.input,file = "parameters.Rdata")
}
parameters.write()
## back to the console
console <- file("console.log", open="wt")
sink(console)
sink(console,append = TRUE,type = "message")

henrik = FALSE


color.blind <- list(	

		yellow = rgb(240,228,66,max = 255),
		orange = rgb(230,159,0,max = 255),
		vermillion = rgb(213,94,0,max = 255),
		bluishgreen = rgb(0,158,115,max = 255),
		reddishpurple = rgb(204,121,167,max = 255), 		skyblue = rgb(86,180,233,max = 255),
		blue = rgb(0,133,178,max = 255)

)


####	
# reading import definitions
####

	##############	GUI
	.label <- "importing import-definitions"
	pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))

	while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))
	}
	##############





if(gui.input$plot.only == ""){
	
	

		##############	GUI
	.label <- "normalizing names"
	pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))

	while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))
	}
	##############
	
	.make.names <- c("sam_id","rawfilename")
	for(i in 1:length(.make.names)){
		temp.i <- grep(.make.names[i],colnames(.data))
		.data[,i] <- make.names(tolower(.data[,i]),allow = FALSE)
		
	##############	GUI
	.label <- "normalizing names"
	pb.check	<- class(try(ui$setProgressBar(pb, prog.max/length(.make.names)*i, label=.label)))

	while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))
	}
	##############
		
	
	}

	##############	GUI
	.label <- "Saving import-binary"
	pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))

	while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- class(try(ui$setProgressBar(pb, 0, label=.label)))
	}
	##############

	save(.data,file = "import-binary.Rdata")
	pb.check	<- class(try(ui$setProgressBar(pb, prog.max, label=.label)))

	
	#if(henrik){
	#	template.grep <- grep("col-n15-bsa-4",as.character(.data$rawfilename),fixed = TRUE)
	#	.data <- .data[-template.grep,]
	#}

	# exclude contaminants and reverse peptides
	if(gui.input$ex.conrev == TRUE){
			try(		conrev.vec <- grep(as.character(import.list$Exclusionstring),.data$code,ignore.case = TRUE))
			if(!exists("conrev.vec")){ conrev.vec <- c()}
		if(is.null(gui.input$cbn.prot) == FALSE){
			
			cbn.target 	<- grep(gui.input$cbn.prot,.data$code)
			if(length(cbn.target) == 0){cbn.target <-  0}
			conrev.vec	<- setdiff(conrev.vec,cbn.target)
		}
		
		if(length(conrev.vec) != 0 ){
			.convec <- .data[conrev.vec,]
			.data 	<- .data[-conrev.vec,]
		}
		
		if(dim(.data)[1] == 0){
			tkmessageBox(title="Message",message=paste("No peptides left after exclusion of contaminants.\nAborting process!"),icon="warning",type="ok")
			stop()
			
		}
	
	}
	#stop()
	####
	# Advanced:
	####
	exclude.samples <- FALSE
	if(exclude.samples){
		try(
		exclude.raw <- read.table(paste(path2,"/excludelist.tab",sep = ""))[,1]
		)	

		if(exists("exclude.raw")){
			exclude.data <- c()
			#.data <- backup
			for( i in 2:length(exclude.raw)){
				temp.i <- exclude.raw[i]
				data.i <- grep(tolower(as.character(temp.i)),tolower(as.character(.data$rawfilename)),fixed = TRUE)	
				exclude.data <- rbind(exclude.data,.data[data.i,])
				if(length(data.i) >0){
				.data 	<- .data[-data.i,]
				}else{	print(paste("Target not found",temp.i))
				}
					print(temp.i)
					print(dim(.data))
			}
		}
		
	}
	#stop()
	####
	# end of Advanced
	####
	
	# running hz.matrix.creator
	gui.input$phospho.string <- import.list$Modifications.identifier
	
	#gui.input$phospho <- TRUE
	gui.input$build.matrix <- TRUE

	if(gui.input$phospho){
		gui.input$build.matrix <- TRUE
		phospho.grep <- grep(gui.input$phospho.string ,.data$Modifications)
		unphospho.pep			<- unique(paste(.data$code[phospho.grep],.data$sequence[phospho.grep],sep = "..")	)
		
		if(all(is.na(.data$Modified.Sequence))){
		.data$code[phospho.grep] <- paste(.data$code[phospho.grep],.data$sequence[phospho.grep],.data$Modifications[phospho.grep],sep = "..")
		}else{
		.data$code[phospho.grep] <- paste(.data$code[phospho.grep],.data$Modified.Sequence,gui.input$phospho.string,sep = "..")
	
		}
	}
	
	
	 sink(type = "message")
	 try(sink())


rm(.design)
#stop()
#stop()
hz.write.unique.prots.seq <- function(.data){
	temp.seq  <- unique(cbind(.data$sequence,.data$rawfilename))

.sequences <- aggregate(temp.seq[,1],list(temp.seq[,2]),length)

temp.prots  <- cbind(.data$code,.data$rawfilename)
.proteins <- aggregate(unique(temp.prots)[,1],list(unique(temp.prots)[,2]),length)


if(all(.proteins[,1]==.sequences[,1])){
	


if(max(nchar(.proteins[,1])) > 5){
temp.mai <- 1+(max(nchar(.proteins[,1]))-1)*0.075

}else{
	temp.mai <- 1
}

temp 		<- rbind(.sequences[,2],.proteins[,2])
colnames(temp) <- .proteins[,1]
temp.range 	<- range(temp)
temp.range[2] <- max(temp.range)+0.15*max(temp.range)
temp.range[1] <- 0


#.proteins[,1] <- paste(c(.proteins[,1],.proteins[,1]),collapse = "")
#stop()

if(dim(temp)[2] > 30){
	width.counts <- 7+(dim(temp)[2]-30)*0.18
	
}else{width.counts <- 7}

if(gui.input$graphic.type == "pdf"){
	pdf("counts.pdf",width = width.counts)
}else{
	postscript("counts.eps",width = width.counts,height =7, paper = "special",onefile = FALSE,horizontal = FALSE)
}

par(mai = c(temp.mai,1,1,1))
barplot2(temp,beside = T,col = c(4,2),main = "unique proteins/sequences per sample" ,las = 2,plot.grid = T,ylim = temp.range,ylab = "n")

legend("topright",c("sequences","proteins"),fill = c(4,2),bg = "#FFFFFF99")

dev.off()

}
}

try(hz.write.unique.prots.seq(.data))
save.image("test.Rdata")
print(dim(.data))
	.error <- class(try(
	.data2 	<- hz.matrix.creator(	.data,
									Raw 			= gui.input$raw,
								type			= gui.input$quant,
								score 			= gui.input$score,
								merge.method 	= gui.input$peptidemerge,	
								outlier 		= gui.input$outlier,
								norm.tog.pep 	= gui.input$peptidenorm,
								shape 			= gui.input$shape,
								row.norm 		= gui.input$row.norm,
								add.data 		= FALSE,#add.data,
								ui 				= ui,
								re.pep.ex 		= gui.input$exclu,
								cbn.prot 		= gui.input$cbn.prot,
								maxq.exclu 		= gui.input$maxq.exclu,
								phospho 		= as.logical(gui.input$phospho),
								action.dupli 	= gui.input$dupli.val,
								prot.norm.shape = gui.input$shape.prot.norm,
								use.raw 		= gui.input$use.raw,
								N15				= gui.input$N15,
								row.target.norm = gui.input$row.target.norm,
								zero.treat		= gui.input$zero.treat ,
								path.design		= gui.input$exp.design,
								n15.correct.method = gui.input$n15.correct.method,
								n15.correct.expect = gui.input$n15.correct.expect,
								n15.log2		= gui.input$n15.log2,
								n.correction	= gui.input$n.correction,
								
								build.matrix 	= gui.input$build.matrix,
								calc.empai		= gui.input$calc.empai,
								empai.sd		= gui.input$empai.sd,
								empai.from.msms	= gui.input$empai.from.msms,
								empai.norm		= gui.input$empai.norm,
								length.pep 		= as.numeric(gui.input$empai.pep.length),
								group.norm		= FALSE, #gui.input$group.norm,
								group.filter	= gui.input$group.filter,
								norm.method		= gui.input$norm.method,
								sum.of.total	= gui.input$sum.of.total,
								graphic.type	= gui.input$graphic.type,
								gui.input		= gui.input,
								prog.max=prog.max,pb=pb,
								.length.matrix 	= .length.matrix
								
								
								)
								
	))
	#assign(".data2",.data2,envir=.GlobalEnv)


if(exists(".data2")){
	print("calculated data successfully")
}	
	
#stop()
	if(.error == "try-error"){
		
	
		ui$messageBox(title="Abort",message="Error in matrix.creator.function.\nCalculation failed!",icon="error",type="ok") ;stop()
	
	}
	if(gui.input$calc.empai){
	#all.empai.backup 	<- .data2$x
	#.data2$x		<- 	hz.shape(all.empai.backup,shape)$shape
		
	}
	
	#stop()
	
	######
	# group scaling
	######
	


if(gui.input$exp.design != ""){
	
	.design  <- read.table(gui.input$exp.design,header = TRUE,sep = "\t")
	.design  <- .design[.design$Include == 1,]
	.design[,2] <- tolower(make.names(.design[,2],allow = F))
	.design[,1] <- tolower(make.names(.design[,1],allow = F))
}	


	


error.control <- class(try(temp.order <- hz.temp.order.fun(gui.input = gui.input,.data2,.design)))
if(error.control == "try-error"){
		ui$messageBox(title="Abort",message="Error in parsing group vector. Continuing without group scaling.",icon="error",type="ok")

}


	if(gui.input$group.norm){
		#.data2$x <-  hz.norm(.data2$x,1,norm = "mean",group = temp.order)$x
	}


	
	try(.info.data <- hz.info.search(.data2,.data,prog.max = prog.max,ui,pb))
	
	if(!exists(".info.data")){
		.info.data <- c()
	}
	
	.data2$proteinlist.info <- .info.data
	########
	########
	.data2$gui.input <- gui.input
	

		#########
		# order tables
		#########
	
	if(exists(".design")){


		if(!as.logical( gui.input$raw)){
							order.vec <- c()
							for(f in 1:length(colnames(.data2$x))){
		
								grep.f <- grep(gsub(" ","",colnames(.data2$x)[f]),.design$Experiment)						
								temp.f <- min(.design$Order[grep.f])
								order.vec <- c(order.vec,temp.f)
							}
							order.vec <- order(order.vec)
		

		}else{
					
					order.vec <- c()
					for(f in 1:length(colnames(.data2$x))){

						grep.f <- grep(gsub(" ","",colnames(.data2$x)[f]),.design$Alternative.name)						
						temp.f <- min(.design$Order[grep.f])
						order.vec <- c(order.vec,temp.f)
					}
					order.vec <- order(order.vec)

		}
				
	}else{
	
	order.vec <- order(colnames(.data2$x))
	}

	.data2$x 		<- .data2$x[,order.vec]
	.data2$x.sd 	<- .data2$x.sd[,order.vec]
	.data2$prot.n 	<- .data2$prot.n[,order.vec]
		print("hui")

	if(exists(".design")){
		
		save(.data2,.data,temp.order,.design,file = "results-binary.Rdata")								
	}else{
		save(.data2,.data,temp.order,file = "results-binary.Rdata")								

	}
	print("hui")
	cat(rep("*",30),"\n")
	
	cat("Starting with Calculation.\n ")
	cat(rep("*",30),"\n")
	
	setwd(wd.write)
	write.csv(.data,file = "raw-peptidelist.csv")
	write.csv(.data2$x,"proteinlist.csv")
	write.csv(.data2$x.sd,"proteinlist-sd.csv")
	write.csv(.data2$peptidelist,"peptidelist.csv")
	write.csv(.data2$proteinlist.info,"INFO-proteinlist.csv")
	#write.csv(.data2$x.raw.merge.peptide.sd,"sd-peptide-raw-merge.csv")
	write.csv(.data2$prot.n,"n-protein.csv")
	
	
	
	if(gui.input$phospho){
		#write.csv(.data2$phospho.peptides.info,"info-phospho-peptides.csv")
		#write.csv(.data2$phospho.ratios,"ratios-phospho-peptides.csv")
		#write.csv(.data2$phospho.peptides,"intensities-phospho-peptides.csv")
	}


}else{
	print("Loaded binary.Rdata")
	
	setwd(dirname(gui.input$plot.only))

}

gui.input <- gui.backup

close(pb)	
ui <- cracker.ui.tk;	# Use TclTk for ui
ui$init();				# Init (loads library)
ratopx <- 10000
pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)



plot.loop <- 2

if(gui.input$calc.empai){
	plot.loop <- 1
}


#try(	
	
#.data2 <- backup
# color code

if(gui.input$plot.only != "" & exists(".design")){
	
	if(gui.input$exp.design != ""){
	
	.design  <- read.table(gui.input$exp.design,header = TRUE,sep = "\t")
	.design  <- .design[.design$Include == 1,]

	.design[,2] <- tolower(make.names(.design[,2],allow = F))
	.design[,1] <- tolower(make.names(.design[,1],allow = F))
}	

exp.design <- .design	
}else{
exp.design <- .data2$exp.des	
}
set.seed(2)

error.try <- class(.error	<- try(results.script.exp.design<- hz.script.exp.design(exp.design = exp.design,gui.input = gui.input, colorblind.set = colorblind.set, color.blind = color.blind,.data2)))
try(.col 						<- results.script.exp.design$.col)
try(hz.exp.des.parse.data2  	<- results.script.exp.design$hz.exp.des.parse.data2)

print("checked experimental design")

if(error.try == "try-error"){
		yesno.answer <- tkmessageBox(title="Message",message=paste("Error in experimental design!",.error,"Do you like to proceed?"),icon="warning",type="yesno")
	if(as.character(yesno.answer) == "no"){stop("stopped Analysis:",.error)}

}


if(dim(hz.exp.des.parse.data2)[2] >3 ){
	
order.control <- hz.merge.control(hz.exp.des.parse.data2[,2],colnames(.data2$x))
if(any(!is.na(order.control))){
	hz.exp.des.parse.data2 <- hz.exp.des.parse.data2[order.control,]
}
 
	if(all(hz.exp.des.parse.data2[,4] != "")){
	order.templates <- c("x","x.sd","prot.n","phospho.ratios","phospho.peptides")
	for(f in 1:length(order.templates)){
		.data2[[order.templates[f]]] <- .data2[[order.templates[f]]][,order(as.numeric(hz.exp.des.parse.data2[,4]))]
		print(order.templates[f])
				print(colnames(.data2[[order.templates[f]]]))
	}#
	}

}

plot.type 	<- 1
if(!gui.input$color.plots & gui.input$barpl){
	#hz.exp.des.parse.data2[,1] <- "white"
	#.col	<- "white"
}else{
	if(length(unique(hz.exp.des.parse.data2[,1])) == 1){
	#		hz.exp.des.parse.data2[,1] <- "lightgrey"
	#.col	<- "lightgrey"
		
		}

	
}
assign("hz.exp.des.parse.data2",hz.exp.des.parse.data2,envir = .GlobalEnv)
error.try <- class(.error<- try(hz.script.plot.main.return <-  hz.script.plot.main(.data2,.data,gui.input, hz.exp.des.parse.data2,.col,.design,y.lab.input = hz.script.y.lab.return,prog.max,ratio.prog,pb,ui, plot.loop,path.data= gui.input$path.data, foldername=foldername, colorblind.set= colorblind.set, color.blind = color.blind)))


if(error.try == "try-error"){
		tkmessageBox(title="Message",message=paste("Error in plotting!",.error),icon="warning",type="ok")
	

}
graphics.off()	

if(gui.input$phospho){
	gui.input$phospho.string <- "Phospho"
	error.try <- class(.error<- try(hz.script.phospho(.data2,.data,gui.input, hz.exp.des.parse.data2,.col,.design,y.lab.input = hz.script.y.lab.return,prog.max,ratio.prog,pb,ui, plot.loop,path.data= gui.input$path.data, foldername,  colorblind.set, color.blind, hz.script.plot.main.return$hz.cracker.anova.return,plot.type,import.list)))
	
if(error.try == "try-error"){

		tkmessageBox(title="Message",message=paste("Error in phospho-peptide analysis!",.error),icon="warning",type="ok")
	}
}	


graphics.off()	
console <- file("console.log", open="wt")

error.try <- class(.error<- try(hz.parameter.report(gui.input,.data2,.data, hz.script.plot.main.return$.report,.design, foldername)))



if(error.try == "try-error"){
		tkmessageBox(title="Message",message=paste("Error in writing report!",.error),icon="warning",type="ok")
	

}	
	
	if(exists(".design")){
try(		save(.data2,.data,temp.order, hz.exp.des.parse.data2,.design,hz.script.plot.main.return,file = paste(path2, foldername,"results-binary.Rdata",sep = "/"))								
)	
}else{
try(		save(.data2,.data,temp.order, hz.exp.des.parse.data2,file = paste(path2, foldername,hz.script.plot.main.return,"results-binary.Rdata",sep = "/")								
))
	}

	
	
close(pb)

if(gui.input$raw){
	.raw <- "raw"	
}else{
	.raw <- "exp"
}

##sink(type="message")
##sink() 

print("ended script")



print("done")	
return(list(.data2 =.data2,path2 = path2,.data=.data,gui.input = gui.input,statistics = hz.script.plot.main.return$hz.cracker.anova.return))

}
