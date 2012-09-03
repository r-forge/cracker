hz.write.design <-
function(x=NULL,input.file,.data,prog.max = 10000){
	if(!exists("prog.max")){prog.max <- 10000}

N15 <- FALSE




.design <-as.data.frame(hz.design(x, ui = ui,.data = .data)) 

wd 		<- getwd()
.try <- try(setwd(x))
.try <- grep("cannot change working directory",fixed = TRUE,.try)
print(.try)
if(length(.try) != 0){
	setwd(x)#; tkm <- tkmessageBox(title = "",
   # message = paste("Please set true data path first!"), icon = "info", type = "ok");return(.try)
    x <- dirname(x)
}else{setwd(wd)}

print("Write Exp")

ED.name <- "experimental-design-cRacker.tab"
it <- 1
while(length(grep(ED.name,list.files(x)))){
	ED.name <- paste("experimental-design-cRacker-",it,".tab",sep = "")
	it <- it+1
}

write.table(.design,normalizePath(paste(x,ED.name,sep = "/")),sep = "\t",
row.names
 = FALSE)
#tclvalue(path4) <- paste(x,"experimental-design-cRacker.txt",sep ="")
#tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))

tkmessageBox(title = "",
    message = paste("Wrote Experimental Design File into",x), icon = "info", type = "ok")
    
    
 temp.fun <- function(){wd <- getwd()
setwd(x)
hz.show.path(x)
setwd(wd)
} 
    
try(

temp.fun()

)    
    
setwd(wd)
#return(as.data.frame(hz.design(x)))

}
