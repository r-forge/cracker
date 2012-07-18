hz.design <-
function(	path=NA,
	id 		= date(),
	exp_id 	= "non.db",
	sam_id	= "NA",
	exp_name = "NA",
	species = "NA",
	stepbystep = FALSE,
	N15 = FALSE, # for N15 Quantitation
	tr = NULL,
	henrik = FALSE	,
	ui = ui,
	prog.max = 10000,
	.data
){
		if(!exists("prog.max")){prog.max <- 10000}
	#import.list <- import.list[import.list$file.type == path2.set$engine,]


	temp.all <- .data#hz.import(import.list = import.list, path = path,ED.ui = T)
	


	temp.all <- cbind(temp.all$rawfilename, temp.all$sam_id)

print(dim(temp.all))
temp.all <-unique((temp.all))
#### 
# cut paths 

	cut.path <-	function(x){
	
			.col		<- strsplit(as.character(x),"\\",fixed = TRUE)
			.cols	<- c()

			ratio.prog <- prog.max/length(.col)
			for(i in 1: length(.col)) {
				temp.f		<- .col[[i]]
				temp.f		<- temp.f[length(temp.f)]
				.cols[i] 	<- temp.f
			}

			return(.cols)
		}
		
temp.all[,1] <- make.names(cut.path(as.vector(temp.all[,1])),allow = F)
temp.all[,2] <- make.names(cut.path(as.vector(temp.all[,2])),allow = F)



temp.all <- cbind(temp.all,rep("1",dim(temp.all)[1]),rep("1",dim(temp.all)[1]),temp.all[,1])

temp.all <- temp.all[order(temp.all[,1]),]
temp.all <- cbind(temp.all,c(1:dim(temp.all)[1]),rep(1,dim(temp.all)[1]),rep(1,dim(temp.all)[1]))

colnames(temp.all) <- c("Name","Experiment","Group","Group.Filter","Alternative.name","Order","Time","Include")

   # try(close(pb.temp))


	text.out <- "Finished Read Out! Used function >hz.msquant.csv<"		
	cat(paste(rep("*",nchar(text.out)),collapse = "",sep = ""),"\n")
	cat(	text.out		,"\n")
	cat(paste(rep("*",nchar(text.out)),collapse = "",sep = ""),"\n")

	ratio.prog <- prog.max/length(files)
	return(temp.all)
}
