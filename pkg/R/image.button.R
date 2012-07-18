image.button <-
function(level,title.help,text.help, bgcol = "#008e8d",image = NULL ,fgcol = "white",size = 1,bf = "normal"){
 	tkbutton(level,image = image, compound  = "right",command = function(){help.tk(title.help,text.help)}) 
 	}
