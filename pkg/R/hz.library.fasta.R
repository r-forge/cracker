hz.library.fasta <-
function(path1){
peptide.length <- 6
####

.exists<- grep("seqinr",library()$results)
if(length(.exists) == 0){
install.packages("seqinr")

}

library(seqinr)
Filters 			<- matrix(c("Fasta file", ".fasta", "All files", "*"),
                  4, 2, byrow = TRUE)

data.input.name 	<- hz.path.empai(path1)#tk_choose.files("Users/henno/downloads/",filters = Filters)
digest 				<- data.input.name[[2]]
data.input.name		<- data.input.name[[1]]


if(length(data.input.name) == 0){stop()}


ui <- cracker.ui.tk;	# Use TclTk for ui
ui$init();				# Init (loads library)
prog.max 	<- 10000
pb 			<- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)

	ui$setProgressBar(pb, 0, label=paste( "Loading Fasta!"))
library(seqinr)

database.1 			<- read.fasta(data.input.name, seqtype = "AA", as.string = TRUE)
#stop()


	ui$setProgressBar(pb, prog.max/2, label=paste( "Saving database!"))

hz.seq 		<-function(x){
temp 		<- x
attr 		<- attributes(temp)
temp 		<- unlist(temp)
temp 		<- paste(temp,collapse="")
return(c(attr,"Seq" =temp))}

database <- lapply(database.1,hz.seq)
name.database<- basename(data.input.name)
name.database <- gsub(".fasta","", name.database)


all.n <- names(database)
seq.vec	<- c()
ratio.prog <- prog.max/length(all.n)
for(i in 1:length(all.n)){
	temp.i <- database[[i]]$Seq
	seq.vec <- c(seq.vec,temp.i)
	
	
	##############	GUI
	pb.check	<- class(try(ui$setProgressBar(pb, ratio.prog*i, label=label.label<- paste("Preparing Sequences" ,i,"/",length(all.n)))))

while(pb.check == "try-error"){
		print("Warning: User closed window!")
		pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
		pb.check	<- try(ui$setProgressBar(pb, ratio.prog*i, label=label.label))		
}
##############	
	
}
database <- cbind(all.n,seq.vec)

print("starting hz.save.control")
name.database <- paste("cRackerSequence-", name.database,sep = "")
hz.save.control(name.database, path =paste(path1,"/",sep = ""),database,"database")
print("hz.save.control")

database <- database.1

	try(ui$setProgressBar(pb, prog.max, label=paste( "Saving database!")))


if(digest == 1){

###
# database rebuild
###
.length.matrix <- c()
for( i in 1:length(names(database))){
	#if(i == 4){ stop()}
try(ui$setProgressBar(pb,prog.max/length(names(database))* i, label=paste( "Digesting Sequences!")))


	temp.i 		<- unlist(database[i][[1]])
	temp.i.seq 	<- tolower(temp.i[1])
	temp.i.split <- strsplit(as.vector(temp.i.seq),"k|r")
	temp.i.nchar <- nchar(unlist(temp.i.split))
	temp.i.nchar[1:(length(temp.i.nchar)-1)] <- temp.i.nchar[1:(length(temp.i.nchar)-1)]+1
		
	.length.vec <- c()
	for(p in c(peptide.length:15)){
		temp.p <- length(temp.i.nchar[temp.i.nchar >=p])
		.length.vec <- c(.length.vec,temp.p)
	}
	names(.length.vec) <- paste("Pep.Length.>=",c(peptide.length:15))
	
	temp.i.final 	<- c(attributes(temp.i)[1:3],.length.vec,nchar(temp.i.seq))
	.length.matrix 	<- rbind(.length.matrix,temp.i.final)
}

close(pb)

name.database <- basename(data.input.name)
name.database <- gsub(".fasta","", name.database[length(name.database)])
name.database <- gsub(".fas","", name.database[length(name.database)])
name.database <- paste("cRackerEmPAI-", name.database,sep = "")

if(exists(".length.matrix")){
hz.save.control(name.database, path =paste(path1,"/",sep = ""),.length.matrix,".length.matrix")
}
}
}
