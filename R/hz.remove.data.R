hz.remove.data <-
function(folder = "mapping-library",path1 = path1){



tt<-tktoplevel()

list.path <- paste(path1,folder,sep ="/")

sel.list.path <- list.files(list.path)

wait <- 1

tt.var.import.list 					<- sel.list.path
	tt.val.import.list				<- tclVar()  
	tclvalue(tt.val.import.list) 	<- sel.list.path[1]
	comboBox 							<- ttkcombobox(tt,values=tt.var.import.list,textvariable = tt.val.import.list,width = 17,state = "readonly")
		
tkgrid(tklabel(tt,text = paste("Delete",folder,":")),comboBox)

Cancel.but 	<- tk2button(tt,text="Stop",command=function() {tclvalue(wait)<-1;tkdestroy(tt)})
	OK.but 		<- tk2button(tt,text="Delete file!",command=function() {tclvalue(wait)<-2;tkdestroy(tt)})
	
	tkbind(tt, "<Return>",function(x){tclvalue(wait)<-2 ; tkdestroy(tt)})
	tkbind(tt, "<Escape>",function(x){tclvalue(wait)<-1 ; tkdestroy(tt)})

tkgrid(OK.but,Cancel.but,columnspan = 2,pady=5)	
tkwait.window(tt)

if(tclvalue(wait) == 2 ){
	unlink(paste(path1,"/",folder,"/",tclvalue(tt.val.import.list),sep = ""))
	e
	
}



}
