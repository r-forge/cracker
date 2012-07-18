hz.path.set <-
function(import.list = import.list,path1){	
	.bg <- "#efefef"
	
	done <- 1
	tt <- tktoplevel(bg = .bg)
	tkwm.resizable(tt, "FALSE","FALSE")

 	tkwm.title(tt,"cRacker")

	
#	if(is.null(image.path) == FALSE){
		error.try 		<- class(try(	.cra.img 	<- tkimage.create("photo", file=paste(path1,"/cRacker-logo-grey.pnm",sep = ""))))
		if(error.try == "try-error"){.cra.img <- ""}



	try(load(paste(path1,"/settings.Rdata",sep = "")))
	if(is.null(settings$import.type)){
		settings$import.type <- import.list$file.type[2]
	}

#	}


 	locationFrame 		<- tkframe(tt,bg = .bg)

 	question4			<- "Browse"
    path4 				<- tclVar("/Users/")
    path.width  <- 10
    entryWidth <- 50
    pad.val = 5


	tt.var.import.list 				<- import.list$file.type
	tt.val.import.list				<- tclVar()  
	tclvalue(tt.val.import.list) 	<- settings$import.type
	comboBox 						<- ttkcombobox(tt,values=tt.var.import.list,textvariable = tt.val.import.list,width = 17,state = "readonly")



    directoryFrame 		<- tkframe(locationFrame,bg = .bg)
  	buttonRcmdr <- function(..., borderwidth, fg, foreground, relief) ttkbutton(...)

    onBrowse1 			<- function(){
        tclvalue(path4) <- tclvalue(tkchooseDirectory(title="Choose your working directory", initialdir = "/Users/",parent = tt))
       }
    onBrowse2 			<- function(){
        tclvalue(path4) <- tclvalue(tkgetOpenFile(title = "Select your file!",parent = tt))
       }


   	browseButton1 <- buttonRcmdr(directoryFrame, text=gettext("Folder", domain="R-Rcmdr"), width= path.width, command=onBrowse1, borderwidth=3	)
   	browseButton2 <- buttonRcmdr(directoryFrame, text=gettext("File", domain="R-Rcmdr"), width= path.width, command=onBrowse2, borderwidth=3)
 

    intro 		<- tkframe(tt,bg = .bg)
    startstop	<- tkframe(tt,bg = .bg)

       
       
 

    locationField <- ttkentry(directoryFrame, width=entryWidth, textvariable=path4,	background = .bg)
    test<- function(){tclvalue(path4) <- tclvalue(path4)  }
	test()
    locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
        command=function(...) tkxview(locationField, ...))
    tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))
    
    
tkgrid(browseButton1,browseButton2,locationField, sticky="w",padx = pad.val,pady = 1)
tkgrid(directoryFrame,padx = pad.val)

tkgrid(tklabel(intro,image= .cra.img,bg = .bg) ,tklabel(intro,text = "Please choose the folder or file containing your peptidelist(s).\n\nThe folder option reads in all files with a defined file.ending (e.g. '.txt').\nTherefor it is recommended to put all files for the folder analysis in a separate folder!"	,bg = .bg))

tkgrid(intro,columnspan = 2)
tkgrid(locationFrame,comboBox,sticky="w"			,columnspan = 1,padx = pad.val,pady = pad.val)





######
#####
######

locationFrame 					<- tkframe(tt,bg = .bg)

 	question4					<- "Browse"
    path.width  <- 10
    entryWidth <- 80
    pad.val = 5


	tt.val.settings				<- tclVar()  
	tclvalue(tt.val.settings) 	<- "default"

	

    directoryFrame 		<- ttklabelframe(locationFrame,text = " Optional: Load parameters from parameters.Rdata ")

  	buttonRcmdr <- function(..., borderwidth, fg, foreground, relief) ttkbutton(...)

    onBrowse2 			<- function(){
        tclvalue(tt.val.settings) <- tclvalue(tkgetOpenFile(title = "Select your file!",parent = tt, initialdir = "anova-p-values.csv/Users"))
       }


   	browseButton2 <- buttonRcmdr(directoryFrame, text=gettext("Browse", domain="R-Rcmdr"), width= path.width, command=onBrowse2, borderwidth=3)
   	
   	locationField <- ttkentry(directoryFrame, width=entryWidth, textvariable= tt.val.settings,	background = .bg)
    test<- function(){tclvalue(tt.val.settings) <- tclvalue(tt.val.settings)  }
	test()
    locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
        command=function(...) tkxview(locationField, ...))
    tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))
    
   	
   	
   	tkgrid(browseButton2,locationField, sticky="w",padx = pad.val,pady = 4)
	tkgrid(directoryFrame, sticky="nw",padx = pad.val,columnspan = 1)
	tkgrid(locationFrame,padx = pad.val,pady = pad.val,columnspan = 2)

#####
#####
#####

cancel.ok.line <- tkframe(tt,bg = .bg)

	Cancel.but 	<- tk2button(cancel.ok.line,text="Stop",command=function() {tclvalue(done)<-1;tkdestroy(tt)})
	OK.but 		<- tk2button(cancel.ok.line,text="Go!",command=function() {tclvalue(done)<-2;tkdestroy(tt)})
	
	tkbind(tt, "<Return>",function(x){tclvalue(done)<-2 ; tkdestroy(tt)})
	tkbind(tt, "<Escape>",function(x){tclvalue(done)<-1 ; tkdestroy(tt)})

tkgrid(OK.but ,tklabel(cancel.ok.line,text = "       ",bg = .bg),Cancel.but,pady = pad.val)

tkgrid(cancel.ok.line,columnspan = 2)
tkwait.window(tt)
settings$import.type <- tclvalue(tt.val.import.list)
try(save(settings,file = paste(path1,"/settings.Rdata",sep = "")))


if(tclvalue(done) == 2 ){
	if(tclvalue(tt.val.settings)!="default"){
			error.return <- class(try(load(tclvalue(tt.val.settings))))

	}else{
	error.return <- "try-error"	
	}
if(error.return != "try-error"){
	#error.return <- class(try(save(	settings,file = paste(path1,"settings/settings.Rdata",sep =""))))
	if( error.return == "try-error"){
		tkmessageBox(text = "Error in loading parameters from parameters.Rdata\nLoading previous saved parameters.")
	}else{print("loaded saved parameters")}
}
	
	
	return(list(path = tclvalue(path4), engine=  tclvalue(tt.val.import.list),settings = tclvalue(tt.val.settings)))
	
	
}else{
	value <- tkmessageBox(message ="User stopped application!",
    icon = "question",title = "Abort")
stop()	
}


}
