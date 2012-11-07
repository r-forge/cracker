hz.info.search <-
function(.data2,.data, prog.max = prog.max,ui,pb){
if(!exists("ratio.prog")){ratio.prog <- 1000}

sys1 <- Sys.time()
x <- rownames(.data2$peptidelist)	
all.peptides <- as.matrix((paste(
			.data$code,
			.data$sequence,
			round(as.numeric(.data$mcr),digits = 0),
			.data$charge,sep = "#"
		)))




temp <- hz.merge.control(all.peptides,x)
if(is.null(rownames(.data2$peptidelist))|all(is.na(rownames(.data2$peptidelist)))){
	x <- paste(rownames(.data2$x)," ",sep = "")

	all.peptides <- as.matrix((paste(
			.data$code)))

}

pb.check <- class(try(ui$setProgressBar(pb,prog.max, label=paste( "Preparing collecting protein information 2"))))


all.peptides <- all.peptides[temp ]
proteins <- .data$code[temp]
.data <- .data[temp,]
raw.data <- 	paste(all.peptides,.data$code,.data$Description,.data$Score,
				.data$Calibrated.mass.relative.error..ppm,
				.data$sequence,.data$rawfilename,sep = "+#+")

na.test <- is.na(proteins)
raw.data <- raw.data[!na.test]
proteins <- proteins[!na.test]


temp.vec <- c(1:length(raw.data))
temp.list <- list()
i=1
test.all <- c()

factor.value <- as.numeric(as.factor(proteins))
factor.value.uni <- unique(factor.value)
a = 1

temp.vec.count <- 0

ratio.prog <- prog.max/length(unique(.data$code))
assign("cracker.counter.temp",1,envir = .GlobalEnv)
assign("cracker.ratio.prog.temp",ratio.prog,envir = .GlobalEnv)

	pb.check <- class(try(ui$setProgressBar(pb,prog.max, label=paste( "Preparing collecting protein information 5"))))

while(temp.vec.count < length(proteins)){
	
temp.vec.start <- temp.vec.count
	pb.check <- class(try(ui$setProgressBar(pb,cracker.ratio.prog.temp*cracker.counter.temp, label=paste( "Counting next batch."))))

while(temp.vec.count-temp.vec.start< 50000){
temp.vec.count <- temp.vec.count +length(factor.value[factor.value == factor.value.uni[a]])
a = a+1
}
print(paste("check peptides",temp.vec.start,"to",temp.vec.count))

temp.vec <- temp.vec[-c(1:temp.vec.count)]
temp.list <- c((temp.vec.start+1):temp.vec.count)
				
test <- aggregate(raw.data[temp.list],by=list(proteins[temp.list]),function(x,ratio.prog=ratio.prog){
	
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb,cracker.ratio.prog.temp*cracker.counter.temp, label=paste( "Collecting Protein Information",cracker.counter.temp))))

		if(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, cracker.ratio.prog.temp*i, label=paste( "Collect additional Information"))))
		}
	##############
	
	assign("cracker.counter.temp",cracker.counter.temp+1,envir = .GlobalEnv)
	if(length(grep("at1g01300.1",x)) > 0){
		assign("x",x,envir = .GlobalEnv)
		#stop()
	}
	
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
			print(paste("Warning for protein",unique(proteins)[i],"! Could not extract single best peptide\nThe first peptide with the highest Score is chosen instead!"))
			
		}
	}
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
test.all <- rbind(test.all,as.matrix(test))
}

test <- as.matrix(test.all)
test <- test[,-1]

colnames(test) <- c("peptide.species","code","description","score","calibrated.mass.relative.error..ppm","sequence","best in rawfile")
.data.info <- test
if(dim(.data.info)[1] > length(unique(proteins))){
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
