help.button <-
function(level,title.help,text.help, buttontext = " ? " ,fgcol = "white",size = 1,bf = "normal"){
 	tk2button(level,text = buttontext,compound = "right",width = 6,command = function(){help.tk(title.help,text.help)}) 
 	}
