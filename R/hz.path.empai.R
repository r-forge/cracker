hz.path.empai <-
function(path1){	
	
	
	done <- 1
	tt <- tktoplevel()

 	tkwm.title(tt,"Choose your path")
 	path.temp<- path1

if(!is.character(path.temp)){path.temp <- tclvalue(path.temp)}
	
#	if(is.null(image.path) == FALSE){
	
	try(.cra.img 	<- tkimage.create("photo", file=paste(path.temp,"/cRacker-logo-grey.pnm",sep = "")))
	if(!exists(".cra.img")){.cra.img <- ""}


#	}


 	locationFrame 		<- tkframe(tt)

 	question4			<- "Choose"
    path4 				<- tclVar("/Users/")
    path.width  <- 16
    entryWidth <- 40
    pad.val = 5








    directoryFrame 		<- tkframe(locationFrame)
    
  
  
    buttonRcmdr <- function(..., borderwidth, fg, foreground, relief) ttkbutton(...)

    onBrowse 			<- function(){
    	Filters 			<- matrix(c("Fasta file", ".fasta", "All files", "*"),
                  4, 2, byrow = TRUE)

        tclvalue(path4) <- tk_choose.files("Users/",filters = Filters)[2]
       }
       


   	browseButton <- buttonRcmdr(directoryFrame, text=gettext(question4, domain="R-Rcmdr"), width= path.width, command=onBrowse, borderwidth=3)
    

    intro 		<- tkframe(tt)
    digest 		<- tkframe(tt)
    startstop	<- tkframe(tt)

     	
 
    cb3 					<- tkcheckbutton(digest)
    
	tt.val.digest.seq 		<- tclVar(0)
	tkconfigure(cb3,variable=tt.val.digest.seq)  
   	tkgrid(tklabel(digest,text="Digest sequence for emPAI reference:\n(Time consuming process but evident for emPAI!)" ),cb3)
    
 

browseButton <- buttonRcmdr(directoryFrame, text=gettext(question4, domain="R-Rcmdr"), width= path.width, command=onBrowse, borderwidth=3)
    locationField <- ttkentry(directoryFrame, width=entryWidth, textvariable=path4)
    test<- function(){tclvalue(path4) <- tclvalue(path4)  }
	test()
    locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
        command=function(...) tkxview(locationField, ...))
    tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))
tkgrid(browseButton,locationField, sticky="w",padx = pad.val)
tkgrid(directoryFrame, sticky="nw",padx = pad.val)

tkgrid(tklabel(intro,image= .cra.img) ,tklabel(intro,text = "Please choose your protein sequence library (fasta format).")	)


tkgrid(intro,columnspan = 6)
tkgrid(digest,columnspan = 6)

tkgrid(locationFrame,sticky="w"			,columnspan = 3,padx = pad.val,pady = pad.val)

	Cancel.but 	<- tk2button(tt,text="Stop",command=function() {tclvalue(done)<-1;tkdestroy(tt)})
	OK.but 		<- tk2button(tt,text="Go!",command=function() {tclvalue(done)<-2;tkdestroy(tt)})
	
	tkbind(tt, "<Return>",function(x){tclvalue(done)<-2 ; tkdestroy(tt)})
	tkbind(tt, "<Escape>",function(x){tclvalue(done)<-1 ; tkdestroy(tt)})

tkgrid(OK.but,Cancel.but,columnspan = 2,pady=5)	

tkwait.window(tt)

if(tclvalue(done) == 2){
	return(list(path = tclvalue(path4), digest=  tclvalue(tt.val.digest.seq)))
	
	
}else{
	value <- tkmessageBox(message ="User stopped application!",
    icon = "question",title = "Abort")

stop()	
}


}
