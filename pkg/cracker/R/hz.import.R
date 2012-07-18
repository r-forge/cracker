hz.import <-
function(import.list, path.data,input.file,N15 = FALSE,ED.ui = FALSE, path2.input.file,prog.max,ui,pb){
test.envi <- environment()
templates <- list.files(path.data,pattern = ".txt")	
if(path2.input.file != ""){
	templates <- path2.input.file
}	
	
	if(import.list$sep == "\\t"){
		import.list$sep = "\t"
	}	
	if(!exists("ui")){
	ui <- cracker.ui.tk;	# Use TclTk for ui
	ui$init();				# Init (loads library)

	require(tcltk)
	prog.max <- 10000
	}
	if(!exists("pb")){
	pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
	}
	
	try(	ui$setProgressBar(pb,0, label=paste("Reading out data,", "file: 0 of ",length(templates),"files | 0 % done"))
	)	
	
		
	
	if(!exists("import.list")){
try(		ui$messageBox(title="Abort",message="import.list is missing.",icon="error",type="ok") 
)
;stop()
	}		
	
	
	
	if(length(templates) > 1){
			message = "More than one file is accessable for input.It is recommended to have only files in your working directory of the same type, containing your peptide lists.\n\nDo you like to proceed the process?"

				value <- tclvalue(tkmessageBox(message =message,
    			icon = "warning", type = "yesno", default = "yes",title = "Warning!"))
    			

    		
    		
    		if(value == "no"){	
    			
try(    			ui$messageBox(title="Abort",message="User aborted .",icon="error",type="ok"))
    			stop()
    			
    			
    			}
}   

if(length(templates) ==0){
		ui$messageBox(title="Abort",message="Raw Data could not be loaded.\nCheck if the raw data is in the right format and within the working directory.",icon="error",type="ok") 

				stop()
    			  			
}

	all.data <- c()
	data.file <- c()
	data.import.all <- c()
	
	
	for(i in 1:length(templates)){
		
		if(length(list.files(path.data)) == 0){
			read.in <- path.data
				
		}else{
			read.in <- paste(path.data,templates[i],sep = "/")
		}
		
		print(read.in)
		
		try(data.file <-	read.csv(read.in,sep = import.list	$sep,dec =(import.list$dec),skip = import.list$skip ,stringsAsFactors = FALSE))
		print(dim(data.file))
		
		
		try(all.data <- rbind(all.data,data.file))
		all.data <- data.file
		if(exists("ui")){
			ratio.prog <- prog.max/length(templates)
try(			ui$setProgressBar(pb, i*ratio.prog, label=paste( "Reading out data:",i,"of",length(templates),"files |",round(i/length(templates))*100,  "% done"))
)		}
			
	
	print(paste("data sucessfully loaded"))
	## data preparation
	grep.all <- c()
	col.start <- 7
	
	print(colnames(all.data))
	
	for(a in col.start:dim(import.list)[2]){
		print(a)
		grep.i	<- import.list[,a] == colnames(all.data)
		grep.i 	<- unique(c(1:dim(all.data)[2])[grep.i])
		if(length(grep.i) == 0){
			grep.i <- NA	
			
		#ui$messageBox(title="Error!",message=paste("Could not find",import.list[import.list$file.type == "maxquant",a],"in",colnames(data.file)[a],".\nPlease check, if your import settings are correct.")	,icon="error",type="ok") ;
		print(colnames(all.data))
		print(import.list[,a])
		
		}
		if(length(grep.i)!=1){
			print(a)
		}
		grep.all 	<- c(grep.all,grep.i)
	}
	print("finished ")
	grep.all.na <- grep.all
	grep.all[is.na(grep.all)] <- 1
	
	data.import <- all.data[grep.all]
	data.import[,is.na(grep.all.na)] <- NA
	
	colnames(data.import) <- colnames(import.list)[col.start:dim(import.list)[2]]
	
	if(N15){
		if(all(is.na(data.import$intensity.2))){
			require(tcltk)
				message = "No second intensity defined for this experiment!\nUse label free (fraction of total) normalisation instead?"

				value <- tclvalue(tkmessageBox(message =message,
    			icon = "question", type = "yesnocancel", default = "yes",title = "Warning!"))
    			
    			
    		if(value == "yes"){	gui.input$quant.method <- "lf"
    							}
    		
    		
    		if(value == "no"){	gui.input$quant.method <- "lf"
    							gui.input$use.raw
    			}
    		
    		
    		if(value == "cancel"){	if(exists("ui")){
				ui$messageBox(title="Abort",message="User stopped application.",icon="error",type="ok") ;stop()
			};stop()}		

			
			
				
		}else{
			intensity.choose 		<- data.import$intensity.2/data.import$intensity.1
			data.import$intensity.1 <- intensity.choose
		
		}
	}
	
	try(data.import.all  <- rbind(data.import.all ,data.import))
	}
	
	
	data.import.all$code 		<- tolower(data.import.all$code)
		
	data.import.all$sam_id 		<- tolower(basename(gsub("\\\\","/", data.import.all$sam_id)))
	data.import.all$rawfilename <- tolower(basename(gsub("\\\\","/", data.import.all$rawfilename)))

	data.import.all$sam_id 		<- gsub("\\\\","/", data.import.all$sam_id)

	data.import.all$rawfilename  	<- sub(".mgf.htm","", data.import.all$rawfilename ,fixed = TRUE)
	data.import.all$rawfilename  	<- sub(".msm.htm","", data.import.all$rawfilename ,fixed = TRUE)

	data.import.all$sam_id 	<- sub(".raw","", data.import.all$sam_id,fixed = TRUE)
	data.import.all$sam_id 	<- sub(".msm.htm","", data.import.all$sam_id,fixed = TRUE)


	if(ED.ui){
	try(close(pb))
    rm(ui,pb)
		
	}
	return(data.import.all)

}
