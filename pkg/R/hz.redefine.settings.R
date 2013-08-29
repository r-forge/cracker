hz.redefine.settings <- 
function(import.list,template){
require("tcltk2")
	tt <- tktoplevel()
 	tkwm.title(tt,"Warning!")
	done <- 1

	tt.var.import.list 				<- import.list$file.type
	tt.val.import.list				<- tclVar()  
	tclvalue(tt.val.import.list) 	<- template
	comboBox 						<- ttkcombobox(tt,values=tt.var.import.list,textvariable = tt.val.import.list,width = 17,state = "readonly")
	tkgrid(tklabel(tt,text = "Import failed: Please try another template:"))
	tkgrid(comboBox,pady = 5)

    cancel.ok.line <- tkframe(tt)

	Cancel.but 	<- tk2button(cancel.ok.line,text="Stop",command=function() {tclvalue(done)<-3;tkdestroy(tt)})
	OK.but 		<- tk2button(cancel.ok.line,text="Go!",command=function() {tclvalue(done)<-1;tkdestroy(tt)})
	Try.anyway 	<- tk2button(cancel.ok.line,text="Try anyway!",command=function() {tclvalue(done)<-2;tkdestroy(tt)})

	
	tkbind(tt, "<Return>",function(x){tclvalue(done)<-2 ; tkdestroy(tt)})
	tkbind(tt, "<Escape>",function(x){tclvalue(done)<-1 ; tkdestroy(tt)})

tkgrid(OK.but ,Try.anyway,Cancel.but,padx = 1,pady = 1)

tkgrid(cancel.ok.line,columnspan = 2)
	
tkwait.window(tt)
done <- tclvalue(done)
return(list(action = done, type = tclvalue(tt.val.import.list))
)
}