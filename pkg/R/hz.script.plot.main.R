hz.script.plot.main <-
function(.data2,.data,gui.input, hz.exp.des.parse.data2,.col,.design,y.lab.input,prog.max, ratio.prog,pb,ui, plot.loop,path.data, foldername, colorblind.set,color.blind){
	.report <- list()
	#save(.data2,.data,gui.input, hz.exp.des.parse.data2,.col,.design,y.lab.input,prog.max, ratio.prog,pb,ui, plot.loop,path.data, foldername, colorblind.set,color.blind,file = "script.plot.main.Rdata")
	if(!exists("ratio.prog")){ratio.prog <- 1000}
for(plot.type in 1:plot.loop){

	if(plot.type == 1){
		dir.create(.setpath<- paste(paste(path.data,foldername,sep = "/"),"/plot-all","-",gui.input$expname,sep = ""))
		setwd(paste(.setpath,sep = "/"))	
	}
	
	if(plot.type == 2){
		setwd("../")
		if(plot.type == 2){
		dir.create(.setpath<- paste(paste(path.data,foldername,sep = "/"),"/plots-significant-changes","-",gui.input$expname,sep = ""))
		setwd(paste(.setpath,sep = "/"))	
	}
		if(exists("hz.cracker.anova.return")){
			try(p.aov <- hz.cracker.anova.return$p.aov)
		}
		if(!exists("p.aov")){
			p.aov <- matrix("empty",nrow = 1,ncol = 3)
		}
		
		aov.exclude <- as.numeric(p.aov[,3]) < gui.input$p.value
		aov.exclude[is.na(aov.exclude)] <- FALSE
		 
		p.all			<- p.aov[aov.exclude,]


		if(is.vector(p.all)){
			p.all 		<- t(as.matrix(p.all))	
		}
		
		collapsed.aov 	<- paste(p.all[,1],collapse = "|")
		grep.aov 		<- grep(collapsed.aov,rownames(.data2$x))
		if(collapsed.aov == ""){ grep.aov <- NULL}

		
		collapsed.data2 <- paste(rownames(.data2$x)[grep.aov],collapse = "|")
		grep.data2 		<- grep(collapsed.data2,p.all[,1])
		p.all			<- p.all[grep.data2,]
		
		if(gui.input$phospho == FALSE){
			.data2.backup 		<- .data2
		}
		.data2$x 			<- .data2$x[grep.aov,]
		.data2$x.sd 		<- .data2$x.sd[grep.aov,]
		.data2$prot.n 		<- .data2$prot.n[grep.aov,]
		try(.data2$proteinlist.info	 <- .data2$proteinlist.info[grep.aov,])
		if(is.vector(.data2$x )){
			.data2$x  					<- t(as.matrix(.data2$x))
			.data2$x.sd  				<- t(as.matrix(.data2$x.sd))
			.data2$prot.n 				<- t(as.matrix(.data2$prot.n))
			.data2$proteinlist.info		<- t(as.matrix(.data2$proteinlist.info))
			
			rownames(.data2$x) <- rownames(.data2.backup$x)[grep.aov]
			rownames(.data2$prot.n) <- rownames(.data2.backup$prot.n)[grep.aov]
			rownames(.data2$x.sd) <- rownames(.data2.backup$x.sd)[grep.aov]
			rownames(.data2$proteinlist.info	) <- rownames(.data2.backup$proteinlist.info)[grep.aov]

			
		}
	
	
	
		if(is.vector(.data2$x)){
			.data2$x <- t(as.matrix(.data2$x))
		}


		if(dim(.data2$x)[1]== 1){
			write("Only one protein with significant changes, plot error!","Error.txt")
			write.csv(t(as.matrix(.data2$proteinlist.info)),"protein.csv")
			write.csv(t(as.matrix(.data2$x.sd)),"protein-sd.csv")
		}
	}	

	if(!exists("pb")){
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
	}

if(dim(.data2$x)[1] > 1){
	
	if(plot.type ==1){
	error.try <- class(.error<- try(hz.cracker.anova.return <-  hz.script.anova(.data2,gui.input, plot.type,as.numeric(prog.max),pb,ui)))
	}
	
		if(error.try == "try-error"){
		print(.error)
		tkmessageBox(title="Message",message=paste("Error in creating barplots!",.error),icon="warning",type="ok")
		stop()
		}
	
	
#	if(!exists("hz.cracker.anova.return")){stop()}

	if(!exists("pb")){
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
	}
	error.try <- class(.error <- try(hz.script.row.plot.space <- hz.script.row.plot(.data2,gui.input,y.lab.input, hz.cracker.anova.return$.aov.new,hz.exp.des.parse.data2,colorblind.set,.col,prog.max,ratio.prog,pb,ui)))

print("Finished hz.script.row.plot")
	if(error.try == "try-error"){
		print(.error)
		tkmessageBox(title="Message",message=paste("Error in creating barplots!",.error),icon="warning",type="ok")
	}

	error.try <- class(.error <- try(hz.script.bp(.data2,gui.input,.col,prog.max,ui,pb)))

	if(error.try == "try-error"){
		print(.error)
		tkmessageBox(title="Message",message=paste("Error in creating sd boxplots!",.error),icon="warning",type="ok")
	}
print("Finished hz.script.bp")

##############	GUI
	pb.check	<- class(try(ui$setProgressBar(pb, 0, label=paste( "Starting additional Statistics. Idle..."))))
	pb.check	<- class(try(ui$setProgressBar(pb, ratio.prog*1, label=paste( "data preparation..."))))
loop.control <- 1
while(pb.check == "try-error"&loop.control < 10){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- try(ui$setProgressBar(pb, ratio.prog*1, label=paste( "data preparation...")))		
		loop.control <- loop.control +1
}

##############	
print("heatmap")

	error.try <- class(.error<- try(hz.script.heatmap.return <- hz.script.heatmap(.data2,gui.input,prog.max,pb,ui,ratio.prog)))
	
print("Finished hz.script.heatmap")	
	
	if(error.try == "try-error"){
				print(.error)

		tkmessageBox(title="Message",message=paste("Error in creating heatmap matrix!",.error),icon="warning",type="ok")
	}
	

	
	#####	
print("pca")
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*3, label=paste( "PCA calculation. Idle..."))))
	
	
	
while(pb.check == "try-error"){

		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*3, label=paste( "PCA calculation. Idle..."))))

}


			
#####
	error.try <- class( .error <- try(hz.script.pca.return <- hz.script.pca(.data2,gui.input, hz.exp.des.parse.data2,prog.max,pb,ui,ratio.prog)))
	.report$report.pca.bpca <- hz.script.pca.return$report.pca.bpca
	
print("Finished hz.script.pca")	

	if(error.try == "try-error"){
				print(.error)

		tkmessageBox(title="Message",message=paste("Error in prinicipal component analysis!",.error),icon="warning",type="ok")
	}

	
	
	error.try <- class(.error<- try(hz.script.heatmap2.return <- hz.script.heatmap2(.data2,gui.input,hz.cracker.anova.return$p.aov, hz.exp.des.parse.data2,.col,colorblind.set,prog.max,pb,ui, plot.type= plot.type,color.blind= color.blind, ratio.prog = ratio.prog)))
	
print("Finished hz.script.heatmap2")		
	
	if(exists("hz.script.heatmap2.return")){.report$report.heatmap.norm <- hz.script.heatmap2.return$heatmap.norm}
	if(error.try == "try-error"){
		print(.error)
		tkmessageBox(title="Message",message=paste("Error in creating heatmap!",.error),icon="warning",type="ok")
	}


		#assign("test",hz.script.heatmap2.return,envir=.GlobalEnv)

	error.try <- class( .error <- try(hz.script.kmeans.return  <- hz.script.kmeans(.data2,gui.input,.design,y.lab.input,colorblind.set,color.blind,hz.script.heatmap2.return$hz.script.hiercl.return$plot.clustering,.col,prog.max,pb,ui)))
	print("Finished hz.script.kmeans")	
	
	
	if(error.try == "try-error"){
				print(.error)

		tkmessageBox(title="Message",message=paste("Error in kmeans analysis!",.error),icon="warning",type="ok")
	}
	

	if(gui.input$do.cor != "none"){	
	

	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*5, label=paste( "Correlation List..."))))

		while(pb.check == "try-error"){

				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*5, label=paste( "Correlation List..."))))
		}
	##############
			
	
	
	
	error.try <- class( .error <- try(hz.script.graph(.data2,gui.input,prog.max,pb,ui)))
	
	print("Finished hz.script.graph")	
	
	if(error.try == "try-error"){
				print(.error)

		tkmessageBox(title="Message",message=paste("Error in correlation analysis!",.error),icon="warning",type="ok")
	}
	}


	
	if(gui.input$go.library != "none"){
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*6, label=paste( "Started Mapping..."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*6, label=paste( "Started Mapping..."))))
		}
	##############
		
	if(exists("hz.script.go.terms.return")){backup.go.input.agg <- hz.script.go.terms.return$backup.go.input.agg}else{backup.go.input.agg <- list()}

	error.try <- class(.error <- try(hz.script.go.terms.return <- hz.script.go.terms(.data2, hz.script.kmeans.return$kmeans.cluster.output, .data2$proteinlist.info,gui.input, hz.script.kmeans.return$kmeans.col, hz.script.kmeans.return$kmeans.at, hz.script.kmeans.return$kmeans.list,prog.max,pb,ui,plot.type,.col,colorblind.set,color.blind, backup.go.input.agg)))
	
	print("Finished hz.script.go.terms")	
						

	if(error.try == "try-error"){
		tkmessageBox(title="Message",message=paste("Error in mapping!",.error),icon="warning",type="ok")
	}
	if(exists("hz.script.go.terms.return")){
	try(		extended.info <- hz.script.go.terms.return$extended.info)
	}


	}
	
	if(plot.type == 1&gui.input$volcano){
		
		error.try <- class(.error<- try(hz.script.volcano(.data2,gui.input, extended.info, colorblind.set, color.blind ,hz.cracker.anova.return,prog.max,pb,ui)))
		
	print("Finished hz.script.volcano")	

		
	if(error.try == "try-error"){
		print(.error)
	stop("data")

		tkmessageBox(title="Message",message=paste("Error in volcano!",.error),icon="warning",type="ok")
	}
	}
	
	
	ui$setProgressBar(pb, ratio.prog*8, label=paste( "Done..."))

	
}
# plot.type


backup <- .data2	
}
return(list(.report=.report, hz.cracker.anova.return= hz.cracker.anova.return))
}
