help.tk <-
function(title.help,text.help){		
		fontHeading		<- tkfont.create(size=10,weight="bold")
		tth 			<- tktoplevel(bg = .bg)
 	 	tkwm.title(tth,paste("help ",title.help,sep = ""))
 	 	if(nchar(title.help)>1){
		tkgrid(tk2label(tth,text=title.help,font = fontHeading,wraplength = 350,background = .bg))}
		tkgrid(tk2label(tth,text=text.help, wraplength = 400, justify = "center",background = .bg))
		tkgrid(tk2button(tth, text = "Close",
          command = function(){tkdestroy(tth)}))
        }
