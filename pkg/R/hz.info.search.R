hz.info.search <-
function(.data2,.data, prog.max = prog.max,ui,pb){
	if(!exists("ratio.prog")){ratio.prog <- 1000}

x <- rownames(.data2$peptidelist)	
all.peptides <- as.matrix((paste(
			.data$code,
			.data$sequence,
			round(as.numeric(.data$mcr),digits = 0),
			.data$charge,sep = " # "
		)))
uni.prot <-unique(rownames(.data2$x))

if(is.null(rownames(.data2$peptidelist))|all(is.na(rownames(.data2$peptidelist)))){
	x <- paste(rownames(.data2$x)," ",sep = "")

	all.peptides <- as.matrix((paste(
			.data$code)))
	uni.prot <- paste(uni.prot," ",sep = "")

}


	
raw.data <- 	cbind(all.peptides,.data$code,.data$Description,.data$Score,
				.data$Calibrated.mass.relative.error..ppm,
				.data$sequence,.data$rawfilename)
raw.data <- unique(raw.data)

raw.data 	<- apply(raw.data,2,function(x){as.character(x)})
merge.vec 	<- hz.merge.control(raw.data[,1],x)
	
.merged.data <- as.data.frame(raw.data[merge.vec,])
colnames(.merged.data) <- c("peptide.species","code","description","score","calibrated.mass.relative.error..ppm","sequence","best in rawfile")



.data.info <- c()

ratio.prog <- prog.max/length(uni.prot)
for(i in 1:length(uni.prot)){
	
	########## GUI 
				
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*i, label=paste( "Collecting Protein Information"))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*i, label=paste( "Collect additional Information"))))
		}
	##############
	
	
	temp.i 	<- .merged.data[as.character(.merged.data[,2]) == uni.prot[i],]
	
	
	temp.i.score 	<- as.numeric(as.character( temp.i$score))
	temp.i.score[is.na(temp.i.score)] <- 0
	
	temp.i.best  <- temp.i[temp.i.score == max(as.numeric(temp.i.score),na.rm = TRUE),]
	
	
	
	if(dim(temp.i.best)[1] > 1){
		print(paste("check mass error for",uni.prot[i]))
		temp.i.ppm 	<- as.numeric(as.character( temp.i.best$calibrated.mass.relative.error..ppm))
		temp.i.ppm[is.na(temp.i.ppm)] <- 0
		temp.i.best <- temp.i.best[temp.i.best$calibrated.mass.relative.error..ppm  == min(as.numeric(temp.i.ppm),na.rm = TRUE),]
		if(dim(temp.i.best)[1] > 1){
			temp.i.best <- temp.i.best[1,]
			print(paste("Warning for protein",uni.prot[i],"! Could not extract single best peptide\nThe first peptide with the highest Score is chosen instead!"))
			
		}
	}
	
	#temp.i.best <- cbind(temp.i.best)
	.data.info <- rbind(.data.info,temp.i.best)
}

if(dim(.data.info)[1] > length(uni.prot)){
	print("Warning! More entries in .data.info than expected!")	
}

.data.info <- .data.info[order(.data.info[,1]),]
.data.info <- apply(.data.info,2, as.character)
.data.info[duplicated(.data.info[,2]),2] <- as.character(paste(.data.info[duplicated(.data.info[,2]),2],paste("search.error.",1:length(.data.info[duplicated(.data.info[,2]),2]),sep = ""),sep = "."))

try(rownames(.data.info) <- .data.info[,2])

return(.data.info)

}
