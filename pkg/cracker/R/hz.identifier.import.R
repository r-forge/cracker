hz.identifier.import <-
function(.array, data.input.name = "mapping"){
.vec.array 	<- as.character(as.vector(.array[1,]))

.vec.array 	<- paste(c(1:length(.vec.array)),.vec.array,sep =": ")

tt3 		<- tktoplevel()


#### NAME
entryInit 	<- data.input.name
  	entryWidth 	<- 20
  	question 	<- "Mapping name"
  
  	tt3.var.name 	<- tclVar(paste(entryInit))
  	textEntryWidget 	<- tkentry(tt3,width=paste(entryWidth),textvariable= tt3.var.name)
  	
tkgrid(tklabel(tt3,text=question),textEntryWidget,pady = 5,padx = 5)
####



tt3.acc		<- tclVar()  
tclvalue(tt3.acc) 	<- .vec.array[1]

tt3.identifier		<- tclVar()  
tclvalue(tt3.identifier) 	<- .vec.array[2]

tt3.type		<- tclVar()  
tclvalue(tt3.type) 	<- .vec.array[3]

comboBox.acc 				<- ttkcombobox(tt3,values= .vec.array,textvariable = tt3.acc,width = 17,state = "readonly")

comboBox.identifier.class 			<- ttkcombobox(tt3,values= .vec.array,textvariable = tt3.identifier,width = 17,state = "readonly")

comboBox.function 			<- ttkcombobox(tt3,values= c(.vec.array,"not available"),textvariable = tt3.type,width = 17,state = "readonly")

tkgrid(tklabel(tt3,text = "Choose columns:"))
tkgrid(tklabel(tt3,text = "Accession"),tklabel(tt3,text = "Mapping"),tklabel(tt3,text = "Type"))
tkgrid(comboBox.acc, comboBox.identifier.class, comboBox.function)

done <- tclVar(0)

Cancel.but 	<- tk2button(tt3,text="Stop",command=function() {tclvalue(done)<-1;tkdestroy(tt3);return()})
OK.but 		<- tk2button(tt3,text="Done",command=function() {tclvalue(done)<-2;tkdestroy(tt3)})
	
	tkbind(tt3, "<Return>",function(x){tclvalue(done)<-2 ; tkdestroy(tt3)})
	tkbind(tt3, "<Escape>",function(x){tclvalue(done)<-1 ; tkdestroy(tt3)})

	
   		
   		

tkgrid(OK.but,Cancel.but,columnspan = 1,pady=5)	


tkwait.window(tt3)

if(as.integer(tclvalue(done)) == 2){

tt3.acc		<- tclvalue(tt3.acc)
tt3.identifier 		<- tclvalue(tt3.identifier)
tt3.type 	<- tclvalue(tt3.type)

.vec <- strsplit(c(tt3.acc,tt3.identifier,tt3.type),": ")
.vec <- c(.vec[[1]][1],.vec[[2]][1],.vec[[3]][1])
.vec[.vec == "not available"] <- NA

print(.vec)
print(as.numeric(.vec[!is.na(.vec)]))
.array <- .array[,as.numeric(.vec[!is.na(.vec)])]
if(dim(.array)[2] == 2){
.array <- cbind(.array,"NA")	
}

.gsub <- "'"
.array <- t(apply(.array,1,function(x){gsub(.gsub,"",x) }))


return(list(array=.array,name = tclvalue(tt3.var.name)))}else{return(list(name = NA))}
}
