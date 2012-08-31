hz.read.parameters <-
function(	image.path = NULL,
								build.matrix = "1",
								.bg = "lightgrey",
								path2,
								path2.set,
								.data,
								path1 
								
								) {
					##d9d9d9
					#.bg = "red"	
					
					require(tcltk)
					require(tcltk2)

					print("test")								
	done 			<- tclVar(0)	
	print("starting GUI")		

(	path1.backup 	<- path1)
(	path2.backup 	<- path2)
	print("test")
	print(getwd())
	
    
license.text <-  paste( "Copyright (C) <2011>  \n<Henrik Zauber; Max Planck Institute for Molecular Plant Physiology>\n
 This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>")
    
 try(	function.file <- read.csv(paste(path1,"/help-file.csv",sep = ""), stringsAsFactors = FALSE)	
)
if(!exists("function.file")){	
	try(	function.file <- read.csv(paste(path1,"/help-file.csv.gz",sep = ""), stringsAsFactors = FALSE)	
)
}

if(!exists("function.file")){	function.file <- ""
}else{
	print("Read function.file successfully.")
}

 tk2font.set("TkDefaultFont",settings= "-family Tahoma -size 10 -weight normal")   
 
    #.bg.col <- "#efefef"
    fontHeading <- tkfont.create(family = "Tahoma",size=10,weight="bold")

	hz.path.set  <-  function(title.q){
    fileDir<-tclvalue(tkchooseDirectory(title=title.q))
    return(ReadAffy(celfile.path = fileDir))}
     
     
########################################################################################
### GUI START #############################################################################
########################################################################################
	
	spacer 	<- "                                                            "
	tt3 	<- tktoplevel(bg =  .bg)
	tkwm.title(tt3,"cRacker")
	tkwm.resizable(tt3, "FALSE","FALSE")
	topMenu <- tk2menu(tt3,background = .bg)
	tkconfigure(tt3,menu=topMenu,background = .bg)
	fileMenu <- tk2menu(topMenu,tearoff=FALSE,background = .bg)
	go.libraryMenu <- tk2menu(topMenu,tearoff=FALSE,background = .bg)	


	tt2 <- tkframe(tt3,bg =  .bg)
	
	
	if(is.null(image.path) == FALSE& 1==0){
		.wd.image 		<- getwd()
		setwd(image.path)
		error.try 		<- class(try(	.cra.img 	<- tkimage.create("photo", file=paste(path1,"/cRacker-logo-grey.pnm",sep = ""))))
		if(error.try == "try-error"){.cra.img <- ""}
		.cra 			<- tk2label(tt2,image = .cra.img,background = .bg)
		setwd(.wd.image)
		image.file <- TRUE
	}else{image.file <- FALSE}
	
		     ########
			# loading:
			#######
			tk.loading <- tktoplevel(background = .bg)

			tkfocus(tk.loading)
			tkgrid(tklabel(tk.loading,text="Loading, please be patient",background = .bg))
			if(image.file){
				tkgrid(tklabel(tk.loading,image=.cra.img,background = .bg))
			}
	
#load settings-files
	if(path2.set$settings == "default" | !exists("path2.set")|exists(".data2")){
		print("loaded settings")
	try(load(paste(path1,"/settings.Rdata",sep = "")))
	print("test")
	}else{
		print("loaded selected settings")

	try(load(path2.set$settings))		
	}
	#tkwm.geometry(tt2,"445x410+20+0")
	if(!exists("settings")){
		settings <- hz.cracker.default.settings()
	}


	info <- paste(license.text)

#tkgrid(image.button(tt2,"cRacker",info,image = .cra.img,fgcol = "white",bgcol = "white",size = 15,bf = "bold"),columnspan = 2,pady = 5)
label.width = 35
#####
# top down level
#####
# libraries


# GO

list.files.mapping <- list.files(paste((path1),sep = ""))



list.files.mapping <- list.files.mapping[grep("^cRackerMapping",list.files.mapping)]

tkadd(	go.libraryMenu,
		"command",
		label="Add",
		command=function(){
			try(hz.go.input(path1))
		tkmessageBox(title="Message",message="Finished Import\nPlease restart cRacker!",icon="warning",type="ok")
		},background = .bg)

tkadd(go.libraryMenu,"command",label="Remove",
	command=function(){hz.remove.data(type = "^cRackerMapping",path1 =tclvalue(path1))},background = .bg)
tkadd(fileMenu,"cascade",label="Mappings",menu=go.libraryMenu,background = .bg)
#emPAI
empai.libraryMenu <- tk2menu(topMenu,tearoff=FALSE)
tkadd(	empai.libraryMenu,
		"command",
		label="Add", 
		command=function(){
			try(hz.library.fasta(path1))
			tkmessageBox(title="Message",message="Finished Import\nPlease restart cRacker!",icon="warning",type="ok")

		})
tkadd(empai.libraryMenu,"command",label="Remove", command=function(){hz.remove.data(type = "^cRackerSequence|^cRackerEmPAI",path1 =tclvalue(path1))},background = .bg)
tkadd(fileMenu,"cascade",label="Protein Sequence library",menu=empai.libraryMenu,background = .bg)


tkadd(topMenu,"cascade",label="Libraries",menu=fileMenu,background = .bg)





######################
# choose quantitation type
###################### 
  	entryInit 	<- "ms-analysis"
  	entryWidth 	<- 20
  	question 	<- "Experiment Name"
  
  	tt2.var.expname 	<- tclVar(paste(entryInit))
  	textEntryWidget 	<- tkentry(tt2,width=paste(entryWidth),textvariable= tt2.var.expname)
  	
  	if(image.file){image.input<- image.button(tt2,"cRacker",info,image = .cra.img,fgcol = "white",bgcol = "white",size = 3,bf = "bold")}else{
  		
  		image.input <- help.button(tt2,"cRacker",buttontext = "cRacker")
  	}
  	
tkgrid(tk2label(tt2,text=question,font = fontHeading,background = .bg,width = label.width),textEntryWidget, image.input,pady = 5,padx = 5)

	
    
    
      
    tt2.experiment.name <- tclVar("1.x")
	
	
	tkgrid(tk2separator(tt2),sticky = "we",columnspan = 3,pady = 5)
######
# Quantitation method
######    



quant.method.button <- "switch to emPAI"
empai.button 	<- tk2button(tt2,text=quant.method.button ,command=function() {tkdestroy(tt3);path1 <- path1.backup;path2 <- path2.backup;tclvalue(done) <- 3})


	tkgrid(tklabel(tt2,text= hz.function.file(function.file,"QM")[2] ,font = fontHeading,width = label.width ,background = .bg,justify = "center",compound = "center"),empai.button,help.button(tt2,hz.function.file(function.file,"QM")[2], hz.function.file(function.file,"QM")[3]))
	tkgrid(tklabel(tt2,text=">> IonIntensities <<",fg = "red",bg = .bg),sticky = "s",columnspan = 3)

########################################################################################################################
### TAB 
########################################################################################################################
tabs <- c(
	tab1 <- "Main",
	tab2 <- "Peptides -> Protein",
	tab3 <- "Extra",
	tab4 <- "Quantitation Mode",
	tab5 <- "Statistics",
	tab5.1 <- "Plotting",
	tab6 <- "Paths"
	)


nb 	<- tk2notebook(tt2, tabs = tabs) # starts notebook

# Make non selected tabs a little darker
#tcl("ttk::style", "configure", "TNotebook.Pane", forground = "blue",background="skyblue3")
#tcl("ttk::style", "map", "TNotebook.Tab", background=c("active", "skyblue"))
#tcl("ttk::style", "map", "TNotebook.Tab", background=c("active", "skyblue"))
#tcl("ttk::style", "configure", "TNotebook.Pane", forground=c("red"))
#tcl("ttk::style", "configure", "TNotebook.Pane", forground=c("red"))

tkgrid(nb,columnspan = 3,padx=5)

# first tab:

	tb1 <- tk2notetab(nb, tab1)
	tb2 <- tk2notetab(nb, tab2)
	tb3 <- tk2notetab(nb, tab3)
	tb4 <- tk2notetab(nb, tab4)
	tb5 <- tk2notetab(nb, tab5)
	tb5.1 <- tk2notetab(nb, tab5.1)
	tb6 <- tk2notetab(nb, tab6)

   pad.val 	<- 2
   pad.y	<- 2
tcl("ttk::style", "configure", "TNotebook", background=.bg)
tcl("ttk::style", "configure", "TNotebook.Tab",foreground = "black", background=.bg)
#tcl("ttk::style", "map", "TNotebook.Tab",foreground = c("focus", "red"))
tcl("ttk::style", "map", "TNotebook.Tab",foreground = c("active", "red"))

#tcl("ttk::style", "configure", "TNotebook", foreground = .bg,background=.bg)
#tcl("ttk::style", "configure", "TNotebook.Pane", background=.bg)

# Make non selected tabs a little darker
#tcl("ttk::style", "configure", "TNotebook.Pane", forground = "blue",background="skyblue3")
#tcl("ttk::style", "map", "TNotebook.Tab", background=c("active", "skyblue"))
#tcl("ttk::style", "map", "TNotebook.Tab", background=c("active", "skyblue"))
tcl("ttk::style", "configure", "TNotebook.Tab", foreground=c("black"))
tcl("ttk::style", "configure", "TNotebook.Pane", foreground=c("red"))

######
# Replicates
###### 


    cb3 					<- tk2checkbutton(tb1)
    
	tb1.val.replicates 		<- tclVar(settings$tb1.val.replicates)
	tkconfigure(cb3,variable=tb1.val.replicates)


tkgrid(tk2label(tb1,text= hz.function.file(function.file,"AR")[2],font = fontHeading,width = label.width ),cb3,help.button(tb1,hz.function.file(function.file,"AR")[2] ,hz.function.file(function.file,"AR")[3]),padx = pad.val, pady = pad.y,sticky = "we")

######
# common peptides:
###### 
    
	help.tb1.replicate.cp 		<- "Only active if \"Replicates\" option is used. An additional normalization factor is created based on the mean relative change of common peptides between replicate samples."
      	
  	cb4 						<- tk2checkbutton(tb1)
	tb1.val.replicates.cp 		<- tclVar(0)#settings$tb1.val.replicates.cp
	tkconfigure(cb4,variable=tb1.val.replicates.cp)
	
#tkgrid(tk2label(tb1,text=title.t <- "normalise on common peptides",width = label.width),cb4,help.button(tb1,title.t ,help.tb1.replicate.cp),padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(tk2label(tb1,text=spacer, width = label.width ),padx = pad.val, pady = pad.y,columnspan = 3)
	

######
# merging Peptides
######     	
####
# protein.intensity
####
	


#tkgrid(tk2label(tb1,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")

	tb1.var.protein.intensity 					<- c("mean","median","sum")
	tb1.val.protein.intensity				<- tclVar()  
	tclvalue(tb1.val.protein.intensity) 	<- settings$tb1.var.protein.intensity	
	comboBox 							<- ttkcombobox(tb1,values=tb1.var.protein.intensity,textvariable = tb1.val.protein.intensity,width = 17,state = "readonly")

tkgrid(tk2label(tb1,text = hz.function.file(function.file,"AM")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb1,hz.function.file(function.file,"AM")[2],  hz.function.file(function.file,"AM")[3]),padx = pad.val, pady = pad.y,sticky = "we")



####
# duplicate peptides
####

#tkgrid(tk2label(tb1,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")

	tb1.duplicates.var 					<- c("mean","sum","max","min","exclude")
	tb1.val.pep.duplicates 				<- tclVar()  
	tclvalue(tb1.val.pep.duplicates) 	<- settings$tb1.duplicates.var	
	comboBox 							<- ttkcombobox(tb1,values=tb1.duplicates.var,textvariable = tb1.val.pep.duplicates,width = 17,state = "readonly")

tkgrid(tk2label(tb1,text=hz.function.file(function.file,"PD")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb1,hz.function.file(function.file,"PD")[2] , hz.function.file(function.file,"PD")[3]),padx = pad.val, pady = pad.y,sticky = "we")

######
# contaminants
###### 



    cb4 				<- tk2checkbutton(tb1)
	tb1.val.conrev 		<- tclVar(settings$tb1.val.conrev)
	tkconfigure(cb4,variable=tb1.val.conrev)


tkgrid(tk2label(tb1,text= hz.function.file(function.file,"ExCR")[2],font = fontHeading,width = label.width ),cb4,help.button(tb1,hz.function.file(function.file,"ExCR")[2] ,hz.function.file(function.file,"ExCR")[3]),padx = pad.val, pady = pad.y,sticky = "we")




#######
# tab 2
#######
pad.val <- pad.val
######
# score
###### 
  	entryWidth <- 40
  	tb2.val.score <- tclVar(settings$tb2.val.score)
  	textEntryWidget <- tkentry(tb2,width=paste(entryWidth),textvariable=tb2.val.score,width = 20)
tkgrid(tk2label(tb2,text=hz.function.file(function.file,"PST")[2],font = fontHeading,width = label.width ),textEntryWidget,help.button(tb2,hz.function.file(function.file,"PST")[2],hz.function.file(function.file,"PST")[3]),padx = pad.val, pady = pad.y,sticky = "we")



######
# Tab2
#####

####
# normalization
####


	cb <- tk2checkbutton(tb2)
	tb2.var.rownorm <- tclVar(settings$tb2.var.rownorm)
	tkconfigure(cb,variable=tb2.var.rownorm)
tkgrid(tk2label(tb2,text=hz.function.file(function.file,"SaS")[2],font = fontHeading,width = label.width ),cb,help.button(tb2,hz.function.file(function.file,"SaS")[2] ,hz.function.file(function.file,"SaS")[3]),padx = pad.val, pady = pad.y,sticky = "we")


####
# norm.method
####
	tb2.var.norm.method					<- c("mean","median","z-score")
	tb2.val.norm.method				<- tclVar()  
	tclvalue(tb2.val.norm.method) 	<- settings$tb2.norm.method	
	comboBox 							<- ttkcombobox(tb2,values=tb2.var.norm.method,textvariable = tb2.val.norm.method,width = 17,state = "readonly")

tkgrid(tk2label(tb2,text=hz.function.file(function.file,"SM")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb2,hz.function.file(function.file,"SM")[2] , hz.function.file(function.file,"SM")[3]),padx = pad.val, pady = pad.y,sticky = "nwe")


####
# shape
####


tb2.var.shape 		<- tclVar(settings$tb2.var.shape)
SliderValueLabel 	<- tk2label(tb2,text=tclvalue(tb2.var.shape))


tkconfigure(SliderValueLabel,textvariable=  tb2.var.shape)
slider <- tkscale(tb2, from=0, to=100,
                    variable= tb2.var.shape,
                   orient="horizontal",
                   resolution = 1,
                   showvalue = T,
                   length = 136,
                   bg = .bg
                  	)
tkgrid(tk2label(tb2,text=hz.function.file(function.file,"PoaN")[2],font = fontHeading,width = label.width ),slider,help.button(tb2,hz.function.file(function.file,"PoaN")[2] ,hz.function.file(function.file,"PoaN")[3]),padx = pad.val, pady = pad.y,sticky = "we")


   	  	
#####
# redundant pep
#####



	
	
	cb 				<- tk2checkbutton(tb2)
	tb2.var.red.pep <- tclVar(settings$tb2.var.red.pep)
	tkconfigure(cb,variable=tb2.var.red.pep)

tkgrid(tk2label(tb2,text=hz.function.file(function.file,"EoRP")[2] ,font = fontHeading,width = label.width),cb,help.button(tb2,hz.function.file(function.file,"EoRP")[2] ,hz.function.file(function.file,"EoRP")[3]),sticky = "s",padx = pad.val, pady = pad.y,sticky = "we")


	tb2.var.db 				<- tclVar()  
	tclvalue(tb2.var.db) 	<- settings$tb2.var.db
	list.files.data <- list.files(paste(path1,sep = ""))
	list.files.data <- list.files.data[grep("cRackerSequence-",list.files.data,fixed = T)]
	list.files.data <- substring(list.files.data,17,nchar(list.files.data))
	
	db 						<- as.character(list.files.data)
	comboBox 				<- ttkcombobox(tb2,values=db,textvariable = tb2.var.db,width = 17,state = "readonly")

tkgrid(tk2label(tb2,text="search duplicates for exclusion"),comboBox,padx = pad.val, pady = pad.y,sticky = "we")

	cbmq 				<- tk2checkbutton(tb2)
	tb2.var.red.pep.mq 	<- tclVar(settings$tb2.var.red.pep.mq)
	tkconfigure(cbmq,variable=tb2.var.red.pep.mq)

tkgrid(tk2label(tb2,text="use list info for exclusion" ),cbmq,sticky = "we")


# tab3
######


####
# PHospho
####
	
	tb4.4 				<- tk2checkbutton(tb3)
	tb3.var.phospho 	<- tclVar(settings$tb3.var.phospho)
tkconfigure(tb4.4,variable=tb3.var.phospho)
tkgrid(tk2label(tb3,text=hz.function.file(function.file,"PP")[2],font = fontHeading,width = label.width ),tb4.4,help.button(tb3,hz.function.file(function.file,"PP")[2] ,hz.function.file(function.file,"PP")[3]),padx = pad.val, pady = pad.y,sticky = "we")


######
# outlier
######      	



#tkgrid(tk2label(tb3,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")

	tb3.outlier 					<- c("none","outside both whiskers","outside lower whisker","top 3")
	tb3.var.outlier 				<- tclVar()  
	tclvalue(tb3.var.outlier) 		<- settings$tb3.outlier
	comboBox 							<- ttkcombobox(tb3,values=tb3.outlier,textvariable = tb3.var.outlier,width = 17,state = "readonly")

tkgrid(tk2label(tb3,text=hz.function.file(function.file,"OE")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb3,hz.function.file(function.file,"OE")[2] , hz.function.file(function.file,"OE")[3]),padx = pad.val, pady = pad.y,sticky = "we")

  	entryWidth 			<- 5
  	tb3.val.zero.treat 		<- tclVar(paste(settings$tb3.val.zero.treat ))
onArgEdit <- function () {
            if (is.integer(tb3.val.zero.treat) == FALSE){
               tclvalue(tb3.val.zero.treat) <- tclvalue(tb3.val.zero.treat)
               tkconfigure(textEntryWidget, textvarible= tb3.val.zero.treat)
            }
            return(tclVar(is.integer(tb3.val.zero.treat))) # Must return TRUE to accept edition!
        }
        
textEntryWidget 	<- tkentry(tb3,width=paste(entryWidth),textvariable= tb3.val.zero.treat,width = 19,validate = "none",validatecommand = onArgEdit)
tkgrid(tk2label(tb3,text=hz.function.file(function.file,"ZH")[2],font = fontHeading,width = label.width),textEntryWidget,help.button(tb3,hz.function.file(function.file,"ZH")[2], hz.function.file(function.file,"ZH")[3]),padx = pad.val, pady = pad.y,sticky = "we" )

####		
# TAB4 
pad.val = pad.val
pad.y = 0.5
######
# raw values
######	

	tb4.rb1 <- tk2radiobutton(tb4)
	tb4.rb2 <- tk2radiobutton(tb4)
	tb4.rb3 <- tk2radiobutton(tb4)
	tb4.rb4 <- tk2radiobutton(tb4)

	
	tb4.quant.method <- tclVar(settings$tb4.quant.method)
    tkconfigure(tb4.rb1,variable=tb4.quant.method,value="lf")
	tkconfigure(tb4.rb2,variable=tb4.quant.method,value="cbn")
	tkconfigure(tb4.rb3,variable=tb4.quant.method,value="15n")

entryWidth 	<- 17

	help.tb4.raw.values <-
	

	tb4.val.raws  					<- c(	"fraction of total (fot)",
											"fot plus n correction (fotn)",
											#"fot 50",
											#"fotn 50",
											"no normalization")
	tb4.val.raw.values  			<- tclVar()  
	tclvalue(tb4.val.raw.values ) 	<- settings$tb4.val.raw.values 
	if(length(grep(settings$tb4.val.raw.values,tb4.val.raws,fixed = TRUE))==0){tclvalue(tb4.val.raw.values ) 	<- tb4.val.raws[1]}

	
	
	comboBox						<- ttkcombobox(tb4,values=tb4.val.raws ,textvariable = tb4.val.raw.values, width = entryWidth ,state = "readonly")

tkgrid(tk2label(tb4,text=hz.function.file(function.file,"LF")[2],font = fontHeading,width = label.width),tb4.rb1,help.button(tb4,hz.function.file(function.file,"LF")[2], hz.function.file(function.file,"LF")[3]),padx = pad.val, pady = pad.y,sticky = "we")
tkgrid(tk2label(tb4,text = ""),comboBox,tk2label(tb4,text = ""))  
#####
# CBN
#####
 
	tb4.p 			<- tk2checkbutton(tb4)
	tb4.var.cbn 	<- tclVar(settings$tb4.var.cbn)
	tkconfigure(tb4.p,variable= tb4.var.cbn)	
tkgrid(tk2label(tb4,text=hz.function.file(function.file,"RPN")[2],font = fontHeading,width = label.width),tb4.rb2,help.button(tb4,hz.function.file(function.file,"RPN")[2],hz.function.file(function.file,"RPN")[3]),padx = pad.val, pady = pad.y ,columnspan = 1,sticky = "we")

  	tb4.help1.1 	<- "Choose if factor is used on raw data, or on scaled peptide data."
	tb4.p2 		<- tk2checkbutton(tb4)
	tb4.var2.cbn 	<- tclVar(settings$tb4.var2.cbn)
	tkconfigure(tb4.p2,variable= tb4.var2.cbn)
	#tkgrid(tk2label(tb4,text="use on scaled peptides"),tb4.p2,help.button(tb4,"row norm",tb4.help1.1),padx = pad.val, pady = pad.y,sticky = "we")	
  	
  	question 	<- "Protein:"

  	tb4.help2.cbn 		<- "Define the Accession code of your reference Protein. Default is BSA."
  	tb4.val.cbn.prot 	<- tclVar(paste(settings$tb4.val.cbn.prot))
  	tb4.textEntryWidget <- tkentry(tb4,width=paste(entryWidth+3),textvariable=tb4.val.cbn.prot)
tkgrid(tk2label(tb4,text=hz.function.file(function.file,"RP")[2]),tb4.textEntryWidget,help.button(tb4,hz.function.file(function.file,"RP")[2],hz.function.file(function.file,"RP")[3]),padx = pad.val, pady = pad.y,sticky = "we")
   	  	

####
# norm.prot shape
####
	
	tb4.val.cbn.shape <- tclVar(settings$tb4.val.cbn.shape)
	SliderValueLabel.tb4 <- tk2label(tb4,text=as.character(tclvalue(tb4.val.cbn.shape)))
	tkconfigure(SliderValueLabel.tb4,textvariable=tb4.val.cbn.shape)
	slider.tb4 <- tkscale(tb4, from=0, to=100,
                    variable=tb4.val.cbn.shape,
                   orient="horizontal",
                   resolution = 1,
                   showvalue = T
                   , length =paste(entryWidth+120),
                   bg = .bg)
tkgrid(tk2label(tb4,text = hz.function.file(function.file,"RPNS")[2] ),slider.tb4,help.button(tb4,hz.function.file(function.file,"RPNS")[2] ,hz.function.file(function.file,"RPNS")[3]  ))
#tkgrid(tk2label(tb4,text=spacer),columnspan = 3)      
            
####### 
# tb4.val.n15
######


	tb4.tb4.val.n15.p 		<- tk2checkbutton(tb4)
	tb4.val.n15 	<- tclVar(settings$tb4.val.n15)
	tkconfigure(tb4.tb4.val.n15.p,variable= tb4.val.n15)
tkgrid(tk2label(tb4,text=hz.function.file(function.file,"LN")[2] ,font = fontHeading,width = label.width  ),tb4.rb3,help.button(tb4,hz.function.file(function.file,"LN")[2] ,hz.function.file(function.file,"LN")[3] ),padx = pad.val, pady = pad.y,sticky = "we")
#####
# log
#####

	tb4.tb4.val.n15.log.p 		<- tk2checkbutton(tb4)
	tb4.val.n15.log 	<- tclVar(settings$tb4.val.n15.log)
	tkconfigure(tb4.tb4.val.n15.log.p,variable= tb4.val.n15.log)
tkgrid(tk2label(tb4,text=hz.function.file(function.file,"LOG")[2] ),tb4.tb4.val.n15.log.p,help.button(tb4,hz.function.file(function.file,"LOG")[2],hz.function.file(function.file,"LOG")[3]),padx = pad.val, pady = pad.y,sticky = "we")

####
# n15.log.correction
####
	
	tb4.var.n15.log.correction 					<- c("none","mean","median")
	tb4.val.n15.log.correction				<- tclVar()  
	tclvalue(tb4.val.n15.log.correction) 	<- settings$tb4.var.n15.log.correction	
	comboBox 							<- ttkcombobox(tb4,values=tb4.var.n15.log.correction,textvariable = tb4.val.n15.log.correction,width = 17,state = "readonly")

tkgrid(tk2label(tb4,text=hz.function.file(function.file,"CorLOG")[2]),comboBox,help.button(tb4,hz.function.file(function.file,"CorLOG")[2] , hz.function.file(function.file,"CorLOG")[3]),padx = pad.val, pady = pad.y,sticky = "we")



######
# ratio
######

###### 


  	entryWidth 		<- 5
  	tb4.val.ratio 		<- tclVar(paste(settings$tb4.val.ratio ))	
onArgEdit <- function () {
            if (is.integer(tb4.val.ratio) == FALSE){
               tclvalue(tb4.val.ratio) <- settings$tb4.val.ratio 
               tkconfigure(textEntryWidget, textvariable= tb4.val.ratio)
            }
            return(tclVar(is.integer(tb4.val.ratio))) # Must return TRUE to accept edition!
        }
        
 textEntryWidget 	<- tkentry(tb4,width=paste(entryWidth),textvariable= tb4.val.ratio,width = 19,validate = "none",validatecommand = onArgEdit)
tkgrid(tk2label(tb4,text=hz.function.file(function.file,"ERLU")[2]),textEntryWidget,help.button(tb4,hz.function.file(function.file,"ERLU")[2], hz.function.file(function.file,"ERLU")[3]),padx = pad.val, pady = pad.y,sticky = "we" )




#tb5.1	<- tk2notebook(tb5, tabs = nb.tb5.names<- c("settings","plots")) # starts notebook


#tb5	<- tk2notetab(tb5.1, nb.tb5.names[1] )
#tb5.1.2	<- tk2notetab(tb5.1, nb.tb5.names[2] )

#tkgrid(tb5.1)


####
# cluster.method
####
tkframe.kmeans <- tk2frame(tb5 )

  	entryWidth 		<- 5
  	tb5.val.kmeans 		<- tclVar(paste(settings$tb5.val.kmeans ))
  	
  	
	onArgEdit <- function () {
            if (is.integer(tb5.val.kmeans) == FALSE){
               tclvalue(tb5.val.kmeans) <- settings$tb5.val.kmeans
               tkconfigure(textEntryWidget, textvariable= tb5.val.kmeans)
            }
            return(tclVar(is.integer(tb5.val.kmeans))) # Must return TRUE to accept edition!
        }
	  	textEntryWidget 	<- tkentry(tb5,width=paste(entryWidth),textvariable= 		tb5.val.kmeans,width = 5,validate = "none",validatecommand = onArgEdit)


	tb5.var.cluster.method 			<- c("k-Means","hclust")
	tb5.val.cluster.method			<- tclVar()  
	tclvalue(tb5.val.cluster.method) 	<- settings$tb5.val.cluster.method
	comboBox 							<- ttkcombobox(tkframe.kmeans,values=tb5.var.cluster.method,textvariable = tb5.val.cluster.method,width = 8,state = "readonly")


tkgrid(tk2label(tkframe.kmeans,text = "Clustering ",font = fontHeading,),comboBox,tk2label(tkframe.kmeans,text = " , cluster:"),columnspan = 3)
tkgrid(tkframe.kmeans, textEntryWidget,help.button(tb5,hz.function.file(function.file,"MTC")[2], hz.function.file(function.file,"MTC")[3]),padx = pad.val, pady = pad.y,sticky = "we")
#stop()

#tkgrid(tk2label(tb5,text=hz.function.file(function.file,"MTC")[2]),comboBox,help.button(tb5,hz.function.file(function.file,"MTC")[2], hz.function.file(function.file,"MTC")[3]),padx = pad.val, pady = pad.y,sticky = "we")

######
# anova.p
###### 

tkframe.anova.log2 <- tk2frame(tb5 )

 	anova.log.cb 			<- tk2checkbutton(tkframe.anova.log2)
    
	tb5.val.anova.log2 		<- tclVar(settings$tb5.val.anova.log2)
	tkconfigure(anova.log.cb,variable=tb5.val.anova.log2)


	entryWidth 		<- 5
	
  	tb5.val.anova.p 		<- tclVar(paste(settings$tb5.val.anova.p))
  	
  	
onArgEdit <- function () {
            if (is.integer(tb5.val.anova.p) == FALSE){
               tclvalue(tb5.val.anova.p) <- settings$tb5.val.anova.p
               tkconfigure(textEntryWidget, textvariable= tb5.val.anova.p)
            }
            return(tclVar(is.integer(tb5.val.anova.p))) # Must return TRUE to accept edition!
        }
        
  	textEntryWidget 	<- tkentry(tb5,width=paste(entryWidth),textvariable= tb5.val.anova.p,width = 19,validate = "none",validatecommand = onArgEdit)
  	
tkgrid(tk2label(tkframe.anova.log2,text = paste(hz.function.file(function.file,"PANOVA")[2],"(log2"),font = fontHeading),anova.log.cb,tk2label(tkframe.anova.log2,text = "):",font= fontHeading),columnspan = 3,sticky = "we" ) 
tkgrid(tkframe.anova.log2,textEntryWidget,help.button(tb5,hz.function.file(function.file,"PANOVA")[2],hz.function.file(function.file,"PANOVA")[3]),padx = pad.val, pady = pad.y,sticky = "we" )

####
# p.adjust
####



	tb5.var.p.adjust 			<- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "fdr", "none")
	tb5.val.p.adjust			<- tclVar()  
	tclvalue(tb5.val.p.adjust) 	<- settings$tb5.var.p.adjust	
	comboBox 							<- ttkcombobox(tb5,values=tb5.var.p.adjust,textvariable = tb5.val.p.adjust,width = 17,state = "readonly")

tkgrid(tk2label(tb5,text=hz.function.file(function.file,"MTC")[2]),comboBox,help.button(tb5,hz.function.file(function.file,"MTC")[2], hz.function.file(function.file,"MTC")[3]),padx = pad.val, pady = pad.y,sticky = "we")

######
# onetailed
###### 

 
    cb3 					<- tk2checkbutton(tb5)
    
	tb5.val.onetailed.ttest 		<- tclVar(settings$val.onetailed.ttest)
	tkconfigure(cb3,variable=tb5.val.onetailed.ttest)


tkgrid(tk2label(tb5,text="Include one tailed t-test",width = label.width ),cb3,help.button(tb5,"Include one tailed t-test" ,""),padx = pad.val, pady = pad.y,sticky = "we")


####
# do.go
####


	tb5.do.go.m				<- tk2checkbutton(tb5)
	tb5.var.do.go 			<- tclVar(settings$tb5.var.do.go)

tkconfigure(tb5.do.go.m,variable=tb5.var.do.go)

#tkgrid(tk2label(tb5,text=hz.function.file(function.file,"MAP")[2],font = fontHeading,width = label.width ),tb5.do.go.m,help.button(tb5,hz.function.file(function.file,"MAP")[2],hz.function.file(function.file,"MAP")[3]),padx = pad.val, pady = pad.y,sticky = "we")



####
# path Goterm 
####
path1 <- tclVar(path1)


	tb5.var.go.term.list 		<-c("none",as.character(gsub("^cRackerMapping-","",list.files.mapping)))
	
	tb5.val.go.term.list				<- tclVar()  
	tclvalue(tb5.val.go.term.list) 	<- settings$tb5.go.term.list	
	comboBox 							<- ttkcombobox(tb5,values=tb5.var.go.term.list,textvariable = tb5.val.go.term.list,width = 17,state = "readonly")

tkgrid(tk2label(tb5,text=hz.function.file(function.file,"MAPL")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb5,hz.function.file(function.file,"MAPL")[2], hz.function.file(function.file,"MAPL")[3]),padx = pad.val, pady = pad.y,sticky = "we")
####
# Correlation Matrix
####
tb5.do.cor.help <-
	"Choose the type of correlation for creating a cross correlation list."


#tkgrid(tk2label(tb5,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")

	tb5.do.cor 					<- c("pearson", "kendall", "spearman","none")
	tb5.val.do.cor				<- tclVar()  
	tclvalue(tb5.val.do.cor) 	<- settings$tb5.do.cor	
	comboBox 							<- ttkcombobox(tb5,values=tb5.do.cor,textvariable = tb5.val.do.cor,width = 17,state = "readonly")

tkgrid(tk2label(tb5,text=hz.function.file(function.file,"CA")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb5,hz.function.file(function.file,"CA")[2], hz.function.file(function.file,"CA")[3]),padx = pad.val, pady = pad.y,sticky = "we")

############
# network
############
do.network <- 0
	tb5.do.network.m				<- tk2checkbutton(tb5)
	tb5.var.do.network 			<- tclVar(0)#tclVar(settings$tb5.var.do.network)
tkconfigure(tb5.do.network.m,variable=tb5.var.do.network)

#tkgrid(tk2label(tb5,text=hz.function.file(function.file,"PN")[2]),tb5.do.network.m,help.button(tb5,hz.function.file(function.file,"PN")[2],hz.function.file(function.file,"PN")[3]),padx = pad.val, pady = pad.y,sticky = "we")




######
# volcano
###### 

 
    cb3 					<- tk2checkbutton(tb5)
    
	tb5.val.volcano 		<- tclVar(settings$tb5.val.volcano)
	tkconfigure(cb3,variable=tb5.val.volcano)


tkgrid(tk2label(tb5,text=hz.function.file(function.file,"VP")[2],font = fontHeading,width = label.width ),cb3,help.button(tb5,hz.function.file(function.file,"VP")[2] ,hz.function.file(function.file,"VP")[3]),padx = pad.val, pady = pad.y,sticky = "we")


##########
# log2 ratio threshold
#########


tb5.var.log2.ratio.thres 	<- tclVar(settings$tb5.var.log2.ratio.thres)
SliderValueLabel 			<- tk2label(tb5,text=tclvalue(tb5.var.log2.ratio.thres))


tkconfigure(SliderValueLabel,textvariable=  tb5.var.log2.ratio.thres)
slider <- tkscale(	tb5, from=0, to=5,
                   variable= tb5.var.log2.ratio.thres,
                   orient="horizontal",
                   resolution = 0.1,
                   showvalue = T,
                   length = 136,
                   bg = .bg
                  	)
tkgrid(tk2label(tb5,text=hz.function.file(function.file,"LRT")[2]),slider,help.button(tb5,hz.function.file(function.file,"LRT")[2] ,hz.function.file(function.file,"LRT")[3]),padx = pad.val, pady = pad.y,sticky = "we")


##########################################################################
########## tb5.1.2.fl.bp <- tk2labelframe(tb5.1.2) statistics
###########################################################################
label.width.temp <- 22
pad.y <- 2
pad.val <- 2
tb5.1.2.fl.bp <- tk2labelframe(tb5.1,text = "type of single protein plots")
tb5.1.2.fl.bp.frame <- tk2frame(tb5.1.2.fl.bp )
#############
	parent.placeholder <- tb5.1.2.fl.bp.frame

	tb5.1.2.fl.bp.frame.plot.type.var 					<- c("barplot","grouped barplot (based on time column)","lineplot","time series (all in one)","time series (single plot)")
	tb5.nb.2.fl.bp.frame.val.plot.type 			<- tclVar()  
	tclvalue(tb5.nb.2.fl.bp.frame.val.plot.type) 	<- settings$tb5.nb.2.fl.bp.frame.val.plot.type	
	
	comboBox 							<- ttkcombobox(parent.placeholder,values=tb5.1.2.fl.bp.frame.plot.type.var,textvariable = tb5.nb.2.fl.bp.frame.val.plot.type,width = .width.temp<-  30,state = "readonly")
print("test")

tkgrid(tk2label(parent.placeholder,text=hz.function.file(function.file,"SPPT")[2],font = fontHeading,width = label.width.temp ),comboBox,help.button(parent.placeholder,hz.function.file(function.file,"SPPT")[2] , hz.function.file(function.file,"SPPT")[3]),padx = pad.val, pady = pad.y,sticky = "wen")

##############
  
  	tb5.1.2.fl.bp.frame.var.xlab 	<- tclVar(settings$tb5.nb.2.fl.bp.frame.var.xlab)
  	textEntryWidget 	<- tk2entry(parent.placeholder,width=.width.temp,textvariable= tb5.1.2.fl.bp.frame.var.xlab)
  	
  	
tkgrid(tk2label(parent.placeholder,text=hz.function.file(function.file,"XL")[2],font = fontHeading,width = label.width.temp),textEntryWidget, help.button(parent.placeholder,hz.function.file(function.file,"XL")[2] , hz.function.file(function.file,"XL")[3]),padx = pad.val, pady = pad.y,sticky = "wen")

##############
#############

tkgrid(tb5.1.2.fl.bp.frame,sticky = "wen",padx = pad.val, pady = pad.y,sticky = "we")
tkgrid(tb5.1.2.fl.bp,columnspan = 3, sticky = "wen",padx = pad.val, pady = pad.y,sticky = "we")

########
########
##########################################################################################
############
# color plots
############

 #settings$tb5.color.plots <- "rainbow"
	tb5.help.color.plots 	<- "Select the color mode for plot output."


#tkgrid(tk2label(tb5,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")

	tb5.color.plots 					<- c("rainbow","colorblind","greytone")
	tb5.var.color.plots				<- tclVar()  
	tclvalue(tb5.var.color.plots) 	<- settings$tb5.var.color.plots	
	comboBox 							<- ttkcombobox(tb5.1,values=tb5.color.plots,textvariable = tb5.var.color.plots,width = 17,state = "readonly")

tkgrid(tk2label(tb5.1,text=hz.function.file(function.file,"CM")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb5.1,hz.function.file(function.file,"CM")[2], hz.function.file(function.file,"CM")[3]),padx = pad.val, pady = pad.y,sticky = "we")

############
# graphic.type 
############

 #settings$tb5.graphic.type <- "rainbow"
	tb5.help.graphic.type 	<- "Graphic type defines the graphic type, that is used for graphic output."
print("test")


#tkgrid(tk2label(tb5,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")

	tb5.graphic.type 					<- c("pdf","eps")
	tb5.var.graphic.type				<- tclVar()  
	tclvalue(tb5.var.graphic.type) 	<- settings$tb5.var.graphic.type	
	comboBox 							<- ttkcombobox(tb5.1,values=tb5.graphic.type,textvariable = tb5.var.graphic.type,width = 17,state = "readonly")

tkgrid(tk2label(tb5.1,text=hz.function.file(function.file,"GT")[2],font = fontHeading,width = label.width ),comboBox,help.button(tb5.1,hz.function.file(function.file,"GT")[2], hz.function.file(function.file,"GT")[3]),padx = pad.val, pady = pad.y,sticky = "we")

####
# only plotting
####
plot.only <- FALSE






    locationFrame 		<- tk2frame(tb5.1)
    path5 				<- tclVar("")
   	buttonRcmdr <- function(..., borderwidth, fg, foreground, relief) ttkbutton(...)

    onBrowse 			<- function(){
        tclvalue(path5) <- tclvalue(tkgetOpenFile(filetypes = "{{Rdata} {.Rdata}} {{All files} *}", initialdir = tclvalue(path2),parent = tt2))
       }
    browseButton <- buttonRcmdr(tb5.1, text=gettext("Browse", domain="R-Rcmdr"), width= 6, command=onBrowse, borderwidth=3)
    locationField <- ttkentry(tb5.1, width=entryWidth, textvariable=path5)
    locationScroll <- ttkscrollbar(tb5.1, orient="horizontal",
        command=function(...) tkxview(tb5.1, ...))
    tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))

tkgrid(tk2label(tb5.1,text= hz.function.file(function.file,"SfR")[2],font = fontHeading,width = label.width ),locationField, browseButton,sticky = "we",padx = pad.val,pady = pad.y,sticky = "we")
######
# tab6 
######

pad.val = 1
tkgrid(tk2label(tb6,text="Paths",font = fontHeading),sticky="w",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(tk2label(tb6,text="selected folders:"),sticky="w",padx = pad.val, pady = pad.y,sticky = "we")

	#path1 <- tclVar(path1)
	path2 <- tclVar(path2)
	path3 <- tclVar(paste(as.character(tclvalue(path1)),"scripts/script.shape.R",sep = ""))

	buttonRcmdr <- function(..., borderwidth, fg, foreground, relief) ttkbutton(...)

  	entryWidth <- 37
  	question1 <- "Browse:"
  	question2 <- "data:"
  	question3 <- "Browse"
  	question4 <- "Browse"

	path.width  <- 16

    locationFrame 		<- tk2labelframe(tb6,text = "cRacker path")
    path1 				<- path1
    directoryFrame 		<- tk2frame(locationFrame)
    onBrowse 			<- function(){
        tclvalue(path1) <- tclvalue(tkchooseDirectory(initialdir = tclvalue(path1),parent = tt2))
       }
    browseButton <- buttonRcmdr(directoryFrame, text=gettext(question1, domain="R-Rcmdr"), width= path.width, command=onBrowse, borderwidth=3)
    locationField <- ttkentry(directoryFrame, width=entryWidth, textvariable=path1)
    test<- function(){tclvalue(path1) <- tclvalue(path1)  }
	test()
    locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
        command=function(...) tkxview(locationField, ...))
    tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))
    
    
#tkgrid(browseButton,locationField, sticky="w",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(directoryFrame, sticky="nw",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(locationFrame, sticky="w",columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")
           

    onBrowse <- function(){
        tclvalue(path2) <- tclvalue(tkchooseDirectory(initialdir = tclvalue(path2),parent = tt2))
        }
    browseButton <- buttonRcmdr(directoryFrame, text=gettext(question2, domain="R-Rcmdr"), width=path.width, command=onBrowse, borderwidth=3)
    locationField <- ttkentry(directoryFrame, width= entryWidth, textvariable=path2)
    test<- function(){tclvalue(path2) <- tclvalue(path2)  }
	test()
     locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
     command=function(...) tkxview(locationField, ...))
     tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))
#tkgrid(browseButton,locationField, sticky="w",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(directoryFrame, sticky="nw",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(locationFrame, sticky="we",padx = pad.val, pady = pad.y,sticky = "we")

####
#checkmark for personal shaping
###
 
#tkgrid(tk2label(tb6,text="       "),sticky="w",columnspan = 1,padx = pad.val, pady = pad.y,sticky = "we")

  
	#locationFrame 	<- tk2labelframe(tb6)
	#directoryFrame 	<- tk2frame(locationFrame)
path.labelframe 		<- tk2labelframe(tb6,text = "cRacker path")
	path.frame	<- tk2frame(path.labelframe)

path1.input <- tclvalue(path1)
path2.input <- tclvalue(path2)


tkgrid(tk2label(path.frame,text = "cracker:"),tk2label(tb6,text = hz.brake.strings(path1.input,70)),sticky = "wn")
tkgrid(tk2label(path.frame,text = "data path:"),tk2label(tb6,text = hz.brake.strings(path2.input,70)),sticky = "wn")

tkgrid(path.frame)
tkgrid(path.labelframe,sticky="ew",columnspan = 3)

#tkgrid(tk2label(tb6,text = "data path:"),tk2label(tb6,text = tclvalue(path2)))
tkgrid(tk2label(tb6,text="Advanced",font = fontHeading),sticky="w",padx = pad.val, pady = pad.y,sticky = "we")


#tkgrid(tk2label(directoryFrame,text="personal shape"),tb6.p,help.button(tb6,"Labeling",tb6.help.shape.script),sticky = "w",columnspan = 3)

#tkgrid(directoryFrame, sticky="nw",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(locationFrame, sticky="w",padx = pad.val, pady = pad.y,sticky = "we")	
  
 	locationFrame <- tk2labelframe(tb6,text = "script peptide filtering")
    directoryFrame <- tk2frame(locationFrame)
    onBrowse <- function(){
          tclvalue(path3) <- tclvalue(tkgetOpenFile(filetypes = "{{R script} {.R}} {{All files} *}",initialdir=tclvalue(path1),parent = tt2))
          }
    browseButton <- buttonRcmdr(directoryFrame, text=gettext(question3, domain="R-Rcmdr"), width=path.width, command=onBrowse, borderwidth=3)
    locationField <- ttkentry(directoryFrame, width= entryWidth, textvariable=path3)
    test<- function(){tclvalue(path3) <- tclvalue(path3)  }
	test()
    locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
    command=function(...) tkxview(locationField, ...))
tkconfigure(locationField, xscrollcommand=function(...) tkset(locationScroll, ...))

  	tb6.help.shape.script 	<- "Use your own shape script, to adjust displacement of peptides that have a certain degree of non measured values. This option replaces the option \"Percentage of allowed NA\"."
	tb6.p 		<- tk2checkbutton(directoryFrame)
	tb6.var.shape.script 	<- tclVar(settings$tb6.var.shape.script)
	tkconfigure(tb6.p,variable= tb6.var.shape.script)

#tkgrid(browseButton,locationField, sticky="we",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(tk2label(directoryFrame,text="personal shape"),tb6.p,help.button(tb6,"Labeling",tb6.help.shape.script),padx = pad.val,pady = pad.y,sticky = "we")


#tkgrid(directoryFrame, sticky="we",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(locationFrame, sticky="we",padx = pad.val, pady = pad.y,sticky = "we")
  
  
#tkgrid(tk2label(tb6,text="       "),sticky="w",columnspan = 1,padx = pad.val, pady = pad.y,sticky = "we")
  
  tkfocus(tb6)
    #locationFrame 		<- tk2frame(tb6)
    path4 				<- tclVar(settings$path4)
    locationFrame 		<- tk2labelframe(tb6,text = "experimental design")
    
    directoryFrame		<- tk2frame(locationFrame)
    
    onBrowse 			<- function(){
        tclvalue(path4) <- tclvalue(tkgetOpenFile(filetypes = "{{tab delimited} {.tab}} {{All files} *}", initialdir = tclvalue(path2),parent = tt2))
       }
    browseButton <- buttonRcmdr(directoryFrame, text=gettext(question4, domain="R-Rcmdr"), width= 10, command=onBrowse)
    locationField <- ttkentry(directoryFrame, textvariable=path4,width = 45)
    test<- function(){tclvalue(path4) <- tclvalue(path4)  }
	test()
    locationScroll <- ttkscrollbar(directoryFrame, orient="horizontal",
        command=function(...) tkxview(locationField, ...))
    tkconfigure(locationField)

####
# group
####

	 temp.fun <- function(){wd <- getwd()
setwd(path2)
hz.show.path()
setwd(wd)
} 


	cb <- tk2checkbutton(directoryFrame)
	tb2.var.groupnorm <- tclVar(settings$tb2.var.groupnorm)#tclVar(settings$tb2.var.groupnorm)
	tkconfigure(cb,variable=tb2.var.groupnorm)
	tkgrid(browseButton,locationField , tk2button(locationFrame,text = "create",width = 6,command = function(){hz.write.design(as.character(tclvalue(path2)),.data);try(temp.fun)}), sticky="w",padx = pad.val, pady = pad.y,sticky = "we")

	tkgrid(tk2label(directoryFrame,text=hz.function.file(function.file,"GSC")[2]),cb,help.button(directoryFrame,hz.function.file(function.file,"GSC")[2] ,hz.function.file(function.file,"GSC")[3]),pady = pad.y,padx = pad.val,sticky = "we")



####
# group filter
####

	cb <- tk2checkbutton(directoryFrame)
	tb2.var.group.filternorm <- tclVar(settings$tb2.var.group.filternorm)#tclVar(settings$tb2.var.group.filternorm)
	tkconfigure(cb,variable=tb2.var.group.filternorm)

	tkgrid(tk2label(directoryFrame,text=hz.function.file(function.file,"GSH")[2]),cb,help.button(directoryFrame,hz.function.file(function.file,"GSH")[2] ,hz.function.file(function.file,"GSH")[3]),pady = pad.y,padx = pad.val,sticky = "we")

tkgrid(directoryFrame,padx = pad.val, pady = pad.y,sticky = "we",columnspan =3)
tkgrid(locationFrame,columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")	
	
#tkgrid(test,sticky ="w",padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(tk2label(tb2,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(tk2label(tb6	,text=spacer),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")
#tkgrid(tk2label(directoryFrame,text="personal shape"),tb6.p,help.button(tb6,"Labeling",tb6.help.shape.script),sticky = "w",columnspan = 3)
	  	
			
######
# save
######    




#tkgrid(tk2label(tt2,text=spacer,width = label.width ),columnspan = 3,padx = pad.val, pady = pad.y,sticky = "we")

	tt2.settings.var 					<- c("current settings","no","load standard settings")
	tt2.settings.method 				<- tclVar()  
	tclvalue(tt2.settings.method) 		<- settings$tt2.settings.var
	comboBox 						<- ttkcombobox(tt2,values=tt2.settings.var,textvariable = tt2.settings.method,width = 17,state = "readonly")

tkgrid(tklabel(tt2,text = hz.function.file(function.file,"SAVE")[2],font = fontHeading,width = label.width ,background = .bg,justify = "center",compound = "center"),comboBox,help.button(tt2,hz.function.file(function.file,"SAVE")[2], hz.function.file(function.file,"SAVE")[3]),pady = 5)  

##############################################################################################################
# End of parameters 
##############################################################################################################

	Cancel.but 	<- tk2button(tt2,text="Stop",command=function() {tclvalue(done)<-1;tkdestroy(tt3)})
	OK.but 		<- tk2button(tt2,text="Start",command=function() {tclvalue(done)<-2;tkdestroy(tt3)})
	
	tkbind(tt2, "<Return>",function(x){tclvalue(done)<-2 ; tkdestroy(tt3)})
	tkbind(tt2, "<Escape>",function(x){tclvalue(done)<-1 ; tkdestroy(tt3)})

	
  
   		

tkgrid(OK.but,Cancel.but,columnspan = 1,pady=5)	

tkdestroy(tk.loading)
tkgrid(tt2)

tkwait.window(tt2)



 	
##############################################################################################################
# End of parameters 
##############################################################################################################
	 tt2.var.expname	<- as.character(tclvalue(tt2.var.expname))
	 settings.old		<- settings



binary.rewrite <- function(x){
	
	x <-tclvalue(x) 
	if(x == 0 | x == "0"){ x <- FALSE}
	if(x == 1 | x == "1"){ x <- TRUE} 
	return(x)
	}
#######
# tab1	 
#######	
	n.correction = FALSE

	  	 settings$tb1.var.protein.intensity <- as.character(tclvalue(tb1.val.protein.intensity))
	 tb1.val.protein.intensity <- as.character(tclvalue(tb1.val.protein.intensity))
     ####
     settings$tb4.val.raw.values <- as.character(tclvalue(tb4.val.raw.values))

     tab.4.var.raw.x <- as.character(tclvalue(tb4.val.raw.values))
     sum.of.total <- "sum"
	if (tab.4.var.raw.x ==tb4.val.raws[3]){tab.4.var.raw.x <- TRUE; n.correction <- FALSE}else{
		
		if(tab.4.var.raw.x ==tb4.val.raws[2]){tab.4.var.raw.x <- FALSE;n.correction <- TRUE}
		if(tab.4.var.raw.x ==tb4.val.raws[1]){tab.4.var.raw.x <- FALSE;n.correction <- FALSE}
		#if(tab.4.var.raw.x ==tb4.val.raws[3]){tab.4.var.raw.x <- FALSE;n.correction <- TRUE; 	sum.of.total <- "sum50"}
		#if(tab.4.var.raw.x ==tb4.val.raws[4]){tab.4.var.raw.x <- FALSE;n.correction <- FALSE; 	sum.of.total <- "sum50"}
	
	}
     tab.4.var.raw.x <- as.character((tab.4.var.raw.x))
	 ######
	 settings$tb1.val.replicates.cp <- as.character(tclvalue(tb1.val.replicates.cp))

	 tb1.val.replicates.cp <- as.character(tclvalue(tb1.val.replicates.cp))
	if (tb1.val.replicates.cp =="1"){tb1.val.replicates.cp <- TRUE}else{tb1.val.replicates.cp <- FALSE}
     tb1.val.replicates.cp <- as.character((tb1.val.replicates.cp)) 
     
     settings$tb1.val.replicates <- as.character(tclvalue(tb1.val.replicates))

	 tb1.val.replicates <- as.character(tclvalue(tb1.val.replicates))
	if (tb1.val.replicates =="1"){tb1.val.replicates <- FALSE}else{tb1.val.replicates <- TRUE}
     tb1.val.replicates <- as.character((tb1.val.replicates)) 
     
      ######
     settings$tb1.duplicates.var = as.character(tclvalue(tb1.val.pep.duplicates)) 
     tb1.val.pep.duplicates = as.character(tclvalue(tb1.val.pep.duplicates)) 
     
     
     #####
     settings$tb1.val.conrev <- as.character(tclvalue(tb1.val.conrev))
     tb1.val.conrev <- as.character(tclvalue(tb1.val.conrev))
	if (tb1.val.conrev =="1"){tb1.val.conrev <- TRUE}else{tb1.val.conrev <- FALSE}
     tb1.val.conrev <- as.character((tb1.val.conrev))
     
#######
# tab2	 
#######	
     tt2.experiment.name  	<- as.character(tclvalue(tt2.experiment.name))
     
	 ##  
	 
	settings$tb2.var.rownorm <- as.character(tclvalue(tb2.var.rownorm))
	 
	tb2.var.rownorm <- as.character(tclvalue(tb2.var.rownorm))
	if (tb2.var.rownorm=="1"){tb2.var.rownorm <- TRUE}else{tb2.var.rownorm = FALSE}
	tb2.var.rownorm <-  as.character(tb2.var.rownorm)


	settings$tb3.var.phospho <- as.character(tclvalue(tb3.var.phospho))	
	tb3.var.phospho <- as.character(tclvalue(tb3.var.phospho))
	if (tb3.var.phospho =="1"){tb3.var.phospho <- TRUE}else{tb3.var.phospho <- FALSE}
    tb3.var.phospho <- as.character(tb3.var.phospho)
    
   
    settings$tb2.val.score 		<- as.numeric(tclvalue(tb2.val.score))
    tb2.val.score 		<- as.numeric(tclvalue(tb2.val.score))
    if(is.na(tb2.val.score)){tb2.val.score <- 0}
    
    
    settings$tb2.var.groupnorm <- tclvalue(tb2.var.groupnorm)
    settings$tb2.var.group.filternorm <- tclvalue((tb2.var.group.filternorm))

    
   	settings$tb2.norm.method	<- tclvalue(tb2.val.norm.method)
    
    
#######
# tab3	 
#######	
	settings$tb3.outlier <- as.character(tclvalue(tb3.var.outlier))
	tb3.var.outlier <- as.character(tclvalue(tb3.var.outlier))
	if(tb3.var.outlier  == tb3.outlier[1] ){tb3.var.outlier ="NA"}
	if(tb3.var.outlier  == tb3.outlier[2] ){tb3.var.outlier ="row"}
	if(tb3.var.outlier  == tb3.outlier[3] ){tb3.var.outlier ="all.below"}
	if(tb3.var.outlier  == tb3.outlier[4] ){tb3.var.outlier ="top.3"}

	
	
		 
	settings$tb2.var.db			<- tclvalue(tb2.var.db)
	db.dat						<- paste("cRackerSequence-",tclvalue(tb2.var.db),sep = "")
	###	if exclusion of redundant = TRUE
	settings$tb2.var.red.pep <- as.character(tclvalue(tb2.var.red.pep))

	tb2.var.red.pep <- as.character(tclvalue(tb2.var.red.pep))
	if (tb2.var.red.pep =="1"){tb2.var.red.pep <- TRUE}else{tb2.var.red.pep <- FALSE}
    tb2.var.red.pep <- as.character((tb2.var.red.pep)) 
    ###	if exclusion of redundant = TRUE
	settings$tb2.var.shape <- as.numeric(tclvalue(tb2.var.shape))
     
     
	## if max qu exclu
	settings$tb2.var.red.pep.mq <- as.character(tclvalue(tb2.var.red.pep.mq))

	tb2.var.red.pep.mq <- as.character(tclvalue(tb2.var.red.pep.mq))
	if (tb2.var.red.pep.mq =="1"){tb2.var.red.pep.mq <- TRUE}else{tb2.var.red.pep.mq <- FALSE}
    tb2.var.red.pep.mq <- as.character(tb2.var.red.pep.mq) 
    ## zero treat
    
   	settings$tb3.val.zero.treat  <-  as.numeric(tclvalue(tb3.val.zero.treat))
   	tb3.val.zero.treat  <-  as.numeric(tclvalue(tb3.val.zero.treat))
   	
   	
   	### build matrix
   		
	## plot.only 

#######
# tab4	 
#######		 

	settings$tb4.quant.method  <- as.character(tclvalue(tb4.quant.method ))
	tb4.quant.method  <- as.character(tclvalue(tb4.quant.method ))
	print("tb4.quant.method")
	if(tb4.quant.method  == "cbn"){tb4.var.cbn <- "1";tb4.val.n15 <- "0"}
	if(tb4.quant.method  == "15n"){tb4.var.cbn <- "0";tb4.val.n15 <- "1"}
	if(tb4.quant.method  == "lf"){tb4.var.cbn <- "0";tb4.val.n15 <- "0"}
print(tb4.quant.method)
	## Reference Protein
	settings$tb4.val.cbn.prot <- as.character(tclvalue(tb4.val.cbn.prot))   
    if(tb4.var.cbn == "1"){tb4.val.cbn.prot <- as.character(tclvalue(tb4.val.cbn.prot)) }else{tb4.val.cbn.prot <- NULL} 
    ## norm ref value	
    tb4.var2.cbn <- as.character(tclvalue(tb4.var2.cbn))
    
    if(tb4.var2.cbn == "1"){tb4.var2.cbn <- TRUE }else{tb4.var2.cbn <- FALSE} 
    
    settings$tb4.val.cbn.shape <- as.numeric(tclvalue(tb4.val.cbn.shape))
    
    
    ## 15n
	
    
    if(tb4.val.n15 == "1"){tb4.val.n15 <- TRUE }else{tb4.val.n15 <- FALSE} 
    ########
    
    settings$tb4.val.n15.log <- as.character(tclvalue(tb4.val.n15.log))
    
    tb4.val.n15.log <- as.character(tclvalue(tb4.val.n15.log))
    
    if(tb4.val.n15.log == "1"){tb4.val.n15.log <- TRUE }else{tb4.val.n15.log <- FALSE}
    #######
        settings$tb4.var.n15.log.correction <- tclvalue(tb4.val.n15.log.correction)#

    tb4.val.n15.log.correction <- tclvalue(tb4.val.n15.log.correction)#
    #######
    settings$tb4.val.ratio <- as.numeric(tclvalue(tb4.val.ratio))
    tb4.val.ratio <- as.numeric(tclvalue(tb4.val.ratio))
    
    ######
    
       
#######
# tab5	 
#######
  
tb5.val.kmeans <- 	tclvalue(tb5.val.kmeans)  

if(tb5.val.kmeans!= "auto"){ 
	tb5.val.kmeans <- as.numeric(tb5.val.kmeans)
		if(is.na(tb5.val.kmeans)){
		tb5.val.kmeans <- "auto"
		}
}
settings$tb5.val.kmeans <- 	tb5.val.kmeans  
 
  	
 	 cracker 			<- as.character(tclvalue(path1))
  	 path.data 			<- as.character(tclvalue(path2)) 
  	 path.script.shape	<- as.character(tclvalue(path3)) 
  	 exp.design			<- as.character(tclvalue(path4)) 


	 settings$tb6.var.shape.script <- as.character(tclvalue(tb6.var.shape.script))
	 tb6.var.shape.script <- as.character(tclvalue(tb6.var.shape.script))
	if (tb6.var.shape.script=="1"){tb6.var.shape.script <- "TRUE"}else{tb6.var.shape.script = "FALSE"}
		tb6.var.shape.script  <-  as.character(tb6.var.shape.script )

settings$tb5.var.p.adjust<- tclvalue(tb5.val.p.adjust)
tb5.val.p.adjust<- tclvalue(tb5.val.p.adjust)
settings$tb5.val.anova.p	<- as.numeric(tclvalue(tb5.val.anova.p	))

tb5.val.anova.p	<- as.numeric(tclvalue(tb5.val.anova.p	))

settings$tb5.var.do.go <- binary.rewrite(tb5.var.do.go)
settings$tb5.go.term.list <- tclvalue(tb5.val.go.term.list)

settings$tb5.val.volcano <- tclvalue(tb5.val.volcano)
settings$tb5.var.color.plots <- tclvalue(tb5.var.color.plots)
settings$tb5.var.graphic.type <- tclvalue(tb5.var.graphic.type)
settings$tb5.val.volcano	<- tclvalue(tb5.val.volcano)
settings$tb5.do.cor <- tclvalue(tb5.val.do.cor)
settings$tb5.var.log2.ratio.thres <- tclvalue(tb5.var.log2.ratio.thres )


tb5.nb.2.fl.bp.frame.val.plot.type<- tclvalue(tb5.nb.2.fl.bp.frame.val.plot.type) 
settings$tb5.nb.2.fl.bp.frame.val.plot.type <- tb5.nb.2.fl.bp.frame.val.plot.type

if(tb5.nb.2.fl.bp.frame.val.plot.type == tb5.1.2.fl.bp.frame.plot.type.var[1] ){
	barpl 	     	= TRUE
	time.grouped	= FALSE
	lineplot.beside	= FALSE

}
if(tb5.nb.2.fl.bp.frame.val.plot.type == tb5.1.2.fl.bp.frame.plot.type.var[2] ){
	barpl 	     	= TRUE
	time.grouped	= TRUE
	lineplot.beside	= FALSE

}
if(tb5.nb.2.fl.bp.frame.val.plot.type == tb5.1.2.fl.bp.frame.plot.type.var[3] ){
	barpl 	     	= FALSE
	time.grouped	= FALSE
	lineplot.beside	= FALSE

}
if(tb5.nb.2.fl.bp.frame.val.plot.type == tb5.1.2.fl.bp.frame.plot.type.var[4] ){
	barpl 	     	= FALSE
	time.grouped	= TRUE
	lineplot.beside	= FALSE

}

if(tb5.nb.2.fl.bp.frame.val.plot.type == tb5.1.2.fl.bp.frame.plot.type.var[5] ){
	barpl 	     	= FALSE
	time.grouped	= TRUE
	lineplot.beside	= TRUE


}

settings$tb5.nb.2.fl.bp.frame.var.xlab <- tclvalue(tb5.1.2.fl.bp.frame.var.xlab)
settings$tb5.val.cluster.method		<- tclvalue(tb5.val.cluster.method)
if(tclvalue(tb5.val.cluster.method) == "hclust"){
	hclust.groups 	     	= TRUE}else{hclust.groups <- FALSE}
settings$tb5.val.anova.log2 <-tclvalue(tb5.val.anova.log2)
settings$val.onetailed.ttest 	<- tclvalue(tb5.val.onetailed.ttest)
####
# tb7	
####







####
# rewrite settings
###
print("save settings")

try(
if(1==1&tclvalue(done) !=3){
tt2.settings.method <- as.character(tclvalue(tt2.settings.method))
if(tt2.settings.method == tt2.settings.var[1]){
try(	save(settings,file= paste(path1.backup,"/settings.Rdata",sep  ="")))
}
if(tt2.settings.method == tt2.settings.var[2]){}
if(tt2.settings.method == tt2.settings.var[3]){
	load(paste(path1.backup,"/settings-standard.Rdata",sep =""))
	try(save(	settings,file = paste(path1.backup,"/settings.Rdata",sep ="")))
tclvalue(done) <- 4

}
})
print("saved settings")



#######
# returning list	 
#######
print(tclvalue(done))
if(tclvalue(done) == 2){     return(ReadAffy(
     				celfile.path 	= list(
     				settings		= settings,

     				cracker 		= cracker,
     				path.data		= path.data,
     				db		 		= db.dat,
					path.shape		= path.script.shape,

					msquant			= tt2.experiment.name,
     				
     				peptidenorm 	= tb1.val.replicates.cp,
     				raw 			= tb1.val.replicates, 
     				peptidemerge 	= tb1.val.protein.intensity,
					dupli.val		= tb1.val.pep.duplicates,
					ex.conrev		= tb1.val.conrev,
					use.raw			= tab.4.var.raw.x, 
					norm.method		= tclvalue(tb2.val.norm.method),
					sum.of.total	= sum.of.total,
					
     				score 			= tb2.val.score,
     				row.norm 		= tb2.var.rownorm,
     				shape 			= as.numeric(tclvalue(tb2.var.shape )),
     				expname			= tt2.var.expname,
     				exclu 			= tb2.var.red.pep,
     				maxq.exclu 		= tb2.var.red.pep.mq,
     				group.norm 		= binary.rewrite(tb2.var.groupnorm),
     				group.filter	= binary.rewrite(tb2.var.group.filternorm),

     				outlier 		= tb3.var.outlier,
     				phospho 		= tb3.var.phospho,
     				zero.treat		= tb3.val.zero.treat,
     				build.matrix	= TRUE,#tb3.var.build.matrix,
					#plot.only		= binary.rewrite(tb3.var.plot.only),
					plot.only 		= tclvalue(path5),

     				cbn.prot		= tb4.val.cbn.prot,
     				shape.prot.norm = as.numeric(tclvalue(tb4.val.cbn.shape)),
     				n.correction	= n.correction,
     				#position = tkwm.geometry(tb1)
     				N15				= tb4.val.n15,
     				row.target.norm = tb4.var2.cbn,
     				n15.log2		= tb4.val.n15.log ,
     				n15.correct.method = tb4.val.n15.log.correction,
     				n15.correct.expect = tb4.val.ratio,

     				
					centers			= tb5.val.kmeans,
					p.adjust.method	= tb5.val.p.adjust,
     				p.value			= tb5.val.anova.p,
     				log2.test		= binary.rewrite(tb5.val.anova.log2),
     				volcano			= binary.rewrite(tb5.val.volcano),
     				onetailed.ttest = binary.rewrite(tb5.val.onetailed.ttest),
					graphic.type	= tclvalue(tb5.var.graphic.type),
					ratio.thres		= tclvalue(tb5.var.log2.ratio.thres),
     				
     				do.GO			= binary.rewrite(tb5.var.do.go),#tb5.do.go.m
     				go.library		= paste("cRackerMapping-",tclvalue(tb5.val.go.term.list),sep = ""),
     				do.cor			= tclvalue(tb5.val.do.cor),
     				do.network		= binary.rewrite(tb5.var.do.network),
     				color.plots		= tclvalue(tb5.var.color.plots),
     				barpl 	     	= barpl, 	     
     				time.grouped  	= time.grouped,
     				lineplot.beside = lineplot.beside,
     				x.xlab			= as.character(tclvalue(tb5.1.2.fl.bp.frame.var.xlab)),
     				exp.design		= exp.design,
     				hclust.groups	= hclust.groups,
     				
     				script.shape 	= tb6.var.shape.script ,
     				quant.method	= tb4.quant.method,
     				calc.empai		= FALSE,
     				empai.sd		= FALSE,
     				empai.from.msms = FALSE,
     				empai.pep.length= 6,
     				empai.reference	= ""
     				
     				     				)))}else{
     				     					
     				     				if(tclvalue(done) == 3){return("switch")}
     									if(tclvalue(done) == 1){return("stopped")}
     				     				
     				     				if(tclvalue(done) == 4){return("reload")}
     				     				}
     				
     				 }
#print(hz.read.parameters(path1 =("/Users/henno/cracker/pkg/data/"),path2 = "",path2.set = list(settings = "",2,3),.data = ""))