hz.go.input <-
function(path1){
Filters 		<- matrix(c("Text", ".txt", "All files", "*"),
                  4, 2, byrow = TRUE)

data.input.name 	<- tk_choose.files(filters = Filters)
print(data.input.name)
if(length(data.input.name) == 0){stop()}



ui <- cracker.ui.tk;	# Use TclTk for ui
ui$init();				# Init (loads library)
prog.max <- 10000
pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
	ui$setProgressBar(pb, 0, label=paste( "Reading mapping library"))


if(length(data.input.name)>1){data.input.name <- data.input.name[2]}
try(
data.input 		<- read.delim(data.input.name,header = FALSE,stringsAsFactors = FALSE,sep = "\t")
)
if(!exists("data.input")|!dim(data.input)[1] >1| !dim(data.input)[2]> 1){
	try(
	
data.input 		<- read.delim(data.input.name,header = FALSE,stringsAsFactors = FALSE,sep = ";")
)}

if(!exists("data.input")|!dim(data.input)[1] >1| !dim(data.input)[2]> 1){
	try(
	
data.input 		<- read.delim(data.input.name,header = FALSE,stringsAsFactors = FALSE,sep = ",")
)}



if(!exists("data.input")|!dim(data.input)[1] >1| !dim(data.input)[2]> 1){
print(paste("could not import data",data.input.name))
}else{
	
if(exists("data.input")){	
print(paste("Successful import of",data.input.name))
print(data.input[1,])
}else{
print("Data could not be readed.")
}
}


ui$setProgressBar(pb,prog.max/2, label=paste( "Saving mapping library!"))
print(data.input.name)

data.input.name <- basename(data.input.name)
data.input.name <- gsub(".txt","", data.input.name)
data.input.name <- gsub(".csv","", data.input.name)
data.input.name <- gsub(".tab","", data.input.name)

print(data.input.name)
.go.list 	<- hz.identifier.import(data.input,data.input.name)
.go 		<- .go.list$array
data.input.name <- .go.list$name
messagetext <- "ok"
if(!is.na(.go.list$name)){
	grep.name <- grep(as.character(data.input.name),list.files(paste(path1,"mapping-library/",sep = "")))
	if(length(grep.name) > 0){
		print("Warning, filename already exists.")
		messagetext <- tk_messageBox("Warning",message = "Filename already exists.\nOverwrite existent file? ",type ="okcancel")
	}

	if(exists("data.input") & messagetext == "ok"){
	try.error <- class(.error <- try(	save(.go,file =paste(path1,"/cRackerMapping-",basename(data.input.name),"",sep = ""), precheck = TRUE)
))
	if(try.error == "try-error"){
		messagetext <- tk_messageBox("Warning",message = paste("Could not save data:\n",.error),type ="ok")

	}
	}
	ui$setProgressBar(pb,prog.max, label=paste( "Done"))
		
}else{
	ui$setProgressBar(pb,prog.max, label=paste( "Abort"))	
}
close(pb)
}
