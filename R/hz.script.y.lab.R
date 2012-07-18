hz.script.y.lab <-
function(gui.input = gui.input,...){
if(gui.input$calc.empai){
	y.lab.input <- "emPAI"
}else{
	y.lab.input <- "intensity"
	if(gui.input$n15.log2){
		y.lab.input <- paste("log2",y.lab.input)
	}
	
	
}
if(gui.input$row.norm){
	scale.string <- paste("scaled:")
}else{
	scale.string <- ""
}
	


if(gui.input$quant.method == "lf"){
	ylab.input.new <- bquote(~.("")^.(scale.string)~.(" ")~.(y.lab.input)/Sigma~.(y.lab.input)[.("sample")])
}

if(gui.input$quant.method == "lf" & gui.input$n.correction){
	ylab.input.new <- bquote(~.("")^.(scale.string)~(.(y.lab.input)/Sigma~.(y.lab.input)[.("sample")])~.("*")~(.("n")/max~("samples(n)")))
}

if(!is.null(gui.input$cbn.prot)){

	ylab.input.new <- bquote(~.("")^.(scale.string)~.(" ")~.(y.lab.input)/.("averaged protein ")*.(y.lab.input)[.(gui.input$cpn.prot)])

}

if(!is.null(gui.input$cbn.prot)& gui.input$n.correction){

	ylab.input.new <- bquote(~.("")^.(scale.string)~(.(y.lab.input)/.("averaged protein"))~.("*")~(.("n")/max~("samples(n)")))

}

if(!is.null(gui.input$cbn.prot)){

	ylab.input.new <- bquote(~.("")^.(scale.string)~.(" ")~.(y.lab.input)/.("averaged protein ")*.(y.lab.input)[.(gui.input$cpn.prot)])

}

if(gui.input$quant.method == "15n"& !gui.input$calc.empai){

	ylab.input.new <- bquote(~.("")^.(scale.string)~.(" ")~.(y.lab.input)/.(y.lab.input)[.("labeled")])

}

if(gui.input$quant.method == "15n"& !gui.input$calc.empai & gui.input$n.correction){

	ylab.input.new <- bquote(~.("")^.(scale.string)~.(" ")~.(y.lab.input)/.(y.lab.input)[.("labeled")]~.("*")~(.("n")/max~("samples(n)")))

}

y.lab.input<- ylab.input.new
return(y.lab.input)
}
