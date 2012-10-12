hz.info.search <-
function(.data2,.data, prog.max = prog.max,ui,pb){
prog.max <- 10000
	
	
	if(!exists("ratio.prog")){ratio.prog <- 1000}

sys1 <- Sys.time()
x <- rownames(.data2$peptidelist)	
all.peptides <- as.matrix((paste(
			.data$code,
			.data$sequence,
			round(as.numeric(.data$mcr),digits = 0),
			.data$charge,sep = "#"
		)))
	
	pb.check <- class(try(ui$setProgressBar(pb,prog.max, label=paste( "Preparing collecting protein information 1"))))


uni.prot <-unique(rownames(.data2$x))


if(is.null(rownames(.data2$peptidelist))|all(is.na(rownames(.data2$peptidelist)))){
	x <- paste(rownames(.data2$x)," ",sep = "")

	all.peptides <- as.matrix((paste(
			.data$code)))
	uni.prot <- paste(uni.prot," ",sep = "")

}
	pb.check <- class(try(ui$setProgressBar(pb,prog.max, label=paste( "Preparing collecting protein information 2"))))

raw.data <- 	paste(all.peptides,.data$code,.data$Description,.data$Score,
				.data$Calibrated.mass.relative.error..ppm,
				.data$sequence,.data$rawfilename,sep = "+#+")

raw.data <- unique(raw.data)
	pb.check <- class(try(ui$setProgressBar(pb,prog.max, label=paste( "Preparing collecting protein information 3"))))

#raw.data 	<- apply(raw.data,2,function(x){as.character(x)})
merge.vec 	<- hz.merge.control(all.peptides,x)
	pb.check <- class(try(ui$setProgressBar(pb,prog.max, label=paste( "Preparing collecting protein information 4"))))
	

.merged.data <-raw.data[merge.vec]
ratio.prog <- prog.max/length(unique(.data$code[merge.vec]))
assign("cracker.counter.temp",1,envir = .GlobalEnv)
assign("cracker.ratio.prog.temp",ratio.prog,envir = .GlobalEnv)

	pb.check <- class(try(ui$setProgressBar(pb,prog.max, label=paste( "Preparing collecting protein information 5"))))

test <- aggregate(.merged.data,by=list(.data$code[merge.vec]),function(x,ratio.prog=ratio.prog){
	
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb,cracker.ratio.prog.temp*cracker.counter.temp, label=paste( "Collecting Protein Information",cracker.counter.temp))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, cracker.ratio.prog.temp*i, label=paste( "Collect additional Information"))))
		}
	##############
	assign("cracker.counter.temp",cracker.counter.temp+1,envir = .GlobalEnv)
	
	
	tes <- apply(as.matrix(x),1,function(x){
		x <- unlist(strsplit(as.character(x),"+#+",fixed = T))
		})
	tes <- t(tes)
	tes <- as.data.frame(tes)
colnames(tes) <- c("peptide.species","code","description","score","calibrated.mass.relative.error..ppm","sequence","best in rawfile")
	temp.i <- tes
	
	temp.i.score 	<- as.numeric(as.character(tes$score))
	temp.i.score[is.na(temp.i.score)] <- 0
	
	temp.i.best  <- temp.i[temp.i.score == max(as.numeric(temp.i.score),na.rm = TRUE),]
	
	
	
	if(dim(temp.i.best)[1] > 1){
		print(paste("check mass error for",unique(tes[,2])))
		temp.i.ppm 	<- as.numeric(as.character( temp.i.best$calibrated.mass.relative.error..ppm))
		temp.i.ppm[is.na(temp.i.ppm)] <- 0
		temp.i.best <- temp.i.best[temp.i.best$calibrated.mass.relative.error..ppm  == min(as.numeric(temp.i.ppm),na.rm = TRUE),]
		if(dim(temp.i.best)[1] > 1){
			temp.i.best <- temp.i.best[1,]
			print(paste("Warning for protein",uni.prot[i],"! Could not extract single best peptide\nThe first peptide with the highest Score is chosen instead!"))
			
		}
	}
		
#print(temp.i.best)
#temp.i.best <- as.character(temp.i.best)
#print(temp.i.best)
#stop()
print(dim(temp.i.best))
temp.i.best <- as.character(unlist(temp.i.best))
if(length(temp.i.best) !=7){
	l.temp <- 7-length(temp.i.best)
	if(l.temp > 0){
	temp.i.best <- c(temp.i.best,rep(NA,l.temp))
	}else{
		temp.i.best <- temp.i.best[1:7]
	}
	
}

return(temp.i.best)
})

test <- as.matrix(test)
test <- test[,-1]

colnames(test) <- c("peptide.species","code","description","score","calibrated.mass.relative.error..ppm","sequence","best in rawfile")
.data.info <- test
if(dim(.data.info)[1] > length(uni.prot)){
	print("Warning! More entries in .data.info than expected!")	
}

.data.info <- .data.info[order(.data.info[,1]),]
.data.info <- apply(.data.info,2, as.character)
try(.data.info[duplicated(.data.info[,2]),2] <- as.character(paste(.data.info[duplicated(.data.info[,2]),2],paste("search.error.",1:length(.data.info[duplicated(.data.info[,2]),2]),sep = ""),sep = "."))
)
try(rownames(.data.info) <- .data.info[,2])

sys2 <- Sys.time()

return(.data.info)

}
