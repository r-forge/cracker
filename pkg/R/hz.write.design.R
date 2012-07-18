hz.write.design <-
function(x=NULL,.data,prog.max = 10000){
	if(!exists("prog.max")){prog.max <- 10000}

N15 <- FALSE
wd 		<- getwd()
.try <- try(setwd(x))
.try <- grep("cannot change working directory",fixed = TRUE,.try)
print(.try)
if(length(.try) != 0){setwd(wd); tkm <- tkmessageBox(title = "",
    message = paste("Please set true data path first!"), icon = "info", type = "ok");return(.try)}else{setwd(wd)}

print("Write Exp")



.design <-as.data.frame(hz.design(x,ui = ui,.data = .data)) 

ED.name <- "experimental-design-cRacker.tab"
it <- 1
while(length(grep(ED.name,list.files(path2)))){
	ED.name <- paste("experimental-design-cRacker-",it,".tab",sep = "")
	it <- it+1
}
write.table(.design,paste(path2,ED.name,sep = "/"),sep = "\t",
row.names
 = FALSE)
#tclvalue(path4) <- paste(x,"experimental-design-cRacker.txt",sep ="")
#tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))

tkmessageBox(title = "",
    message = paste("Wrote Experimental Design File into",x), icon = "info", type = "ok")
    
    
 temp.fun <- function(){wd <- getwd()
setwd(path2)
hz.show.path(x)
setwd(wd)
} 
    
try(

temp.fun()

)    
    
setwd(wd)
#return(as.data.frame(hz.design(x)))

}
