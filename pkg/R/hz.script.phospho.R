hz.script.phospho <-
function(.data2,.data,gui.input, hz.exp.des.parse.data2,.col,.design,y.lab.input = hz.script.y.lab.return,prog.max,ratio.prog,pb,ui, plot.loop,path.data= gui.input$path.data, foldername,  colorblind.set, color.blind, hz.cracker.anova.return,plot.type,import.list){
	
	
#hz.cracker.anova.return <- hz.script.return$statistics	
#plot.type <- 1
#.col <- 1
#y.lab.input <- "test"
#colorblind.set <- c(1,2,3,4)
#import.list <- import.list$import.list
#.data2 <- backup
#backup<- .data2 
#.data2 <- .data2  		
#stop()

grep.p <- grep(tolower(import.list$Modifications.identifier),tolower(rownames(.data2$x)))
grep.T <- grep(tolower(import.list$Modifications.identifier),tolower(rownames(.data2$x)),invert = T)
matrix.ref <- .data2
matrix.ref$x <- matrix.ref$x[grep.T,]
matrix.ref$x.sd <- matrix.ref$x.sd[grep.T,]
matrix.ref$prot.n <- matrix.ref$prot.n[grep.T,]

if(length(grep.p) > 0){
.data2$x <- .data2$x[grep.p,]
if(!gui.input$n15.log2){
	.data2$x[is.na(.data2$x)] <- 0
}else{
	.data2$x[is.na(.data2$x)] <- min(.data2$x,na.rm = T)
}
.data2$x.sd <- .data2$x.sd[grep(import.list$Modifications.identifier,rownames(.data2$x.sd)),]
.data2$prot.n <- .data2$prot.n[grep(import.list$Modifications.identifier,rownames(.data2$prot.n)),]

dir.create(.setpath <- paste(gui.input$path.data,foldername,"phosphopeptides-all",sep = "/")) 
setwd(.setpath)	



ui <- cracker.ui.tk;	# Use TclTk for ui
ui$init();				# Init (loads library)
prog.max <- 10000
pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)



assign("hz.exp.des.parse.data2",hz.exp.des.parse.data2,envir = .GlobalEnv)
#error.try <- class(.error<- try(hz.script.plot.main.return <-  hz.script.plot.main(.data2,.data,gui.input, hz.exp.des.parse.data2,.col,.design,y.lab.input = hz.script.y.lab.return,prog.max,ratio.prog,pb,ui, plot.loop = 1,path.data= gui.input$path.data, foldername=foldername, colorblind.set= colorblind.set, color.blind = color.blind)))


	#error.try <- class(.error<- try(hz.script.heatmap.return <- hz.script.heatmap(.data2,gui.input,prog.max,pb,ui,ratio.prog)))
	error.try <- class( .error <- try(hz.script.pca.return <- hz.script.pca(.data2,gui.input, hz.exp.des.parse.data2,prog.max,pb,ui,ratio.prog)))

	error.try <- class(.error<- try(hz.script.heatmap2.return <- hz.script.heatmap2(.data2,gui.input,hz.cracker.anova.return$p.aov, hz.exp.des.parse.data2,.col = 1,colorblind.set,prog.max,pb,ui, plot.type= plot.type,color.blind= color.blind, ratio.prog = ratio.prog)))

	error.try <- class( .error <- try(hz.script.kmeans.return  <- hz.script.kmeans(.data2,gui.input,.design,y.lab.input,colorblind.set,color.blind,hz.script.heatmap2.return$hz.script.hiercl.return$plot.clustering,.col,prog.max,pb,ui)))

	error.try <- class(.error <- try(hz.script.row.plot.space <- hz.script.row.plot(.data2,gui.input,y.lab.input, hz.cracker.anova.return$.aov.new,hz.exp.des.parse.data2,colorblind.set,.col,prog.max,ratio.prog=1000,pb,ui)))

#stop()
print("test")
	#error.try <- class( .error <- try(hz.script.graph(.data2,gui.input)))
setwd("../")


p.split <- strsplit(rownames(.data2$x),"..",fixed = T)

phospho.matrix <- c()
peptidelist.grep <- c()
proteinlist.grep <- c()
sdlist.grep <- c()
for(o in 1:length(rownames(.data2$x))){
	temp.o <- p.split[[o]]
	phospho.matrix <- rbind(phospho.matrix,temp.o)
	
	peptidelist.grep <- c(peptidelist.grep,grep(paste(temp.o[1],temp.o[2],sep = "&"),rownames(.data2$peptidelist)))
	proteinlist.grep <- c(proteinlist.grep,grep(paste(temp.o[1]),rownames(.data2$x[-grep(import.list$Modifications.identifier,rownames(.data2$x)),])))
		
}

#.data2$x[is.na(.data2$x)] <- 0
## norm on protein:

ref.protein <- matrix.ref$x
ref.protein[is.na(ref.protein)] <- 0

ref.proteins <- intersect(phospho.matrix[,1],rownames(ref.protein))

ref.proteins.phospho.matrix <- c()
ref.proteins.phospho.matrix.row <- c()
prots.phospho <- c()
if(length(ref.proteins) !=0){
for(i in ref.proteins){
	temp.i 		<- grep(i, rownames(.data2$x))
	prots.phospho <- c(prots.phospho,temp.i)
	
	ref.proteins.phospho.matrix.row <- c(ref.proteins.phospho.matrix.row, rownames(.data2$x)[temp.i])

	temp.i 		<- as.numeric(.data2$x[temp.i,])
	temp.i[is.na(temp.i)] <- 0
	temp.i.ref <- as.numeric(ref.protein[i,])
	temp.i.ref[is.na(temp.i.ref)] <- 0

	temp.phospho 	<- 	(temp.i/temp.i.ref)
	ref.proteins.phospho.matrix <- rbind(ref.proteins.phospho.matrix,temp.phospho)


}
print("test")
rownames(ref.proteins.phospho.matrix) <- ref.proteins.phospho.matrix.row
colnames(ref.proteins.phospho.matrix) <- colnames(.data2$x)

.data2$x <- ref.proteins.phospho.matrix

inf.matrix <- matrix("",dim(ref.proteins.phospho.matrix)[1],dim(ref.proteins.phospho.matrix)[2])
inf.matrix[is.infinite(ref.proteins.phospho.matrix)& ref.proteins.phospho.matrix > 0] <- "Inf"
inf.matrix[is.infinite(ref.proteins.phospho.matrix)& ref.proteins.phospho.matrix < 0] <- "-Inf"

stop()
.data2$x.sd[is.na(.data2$x.sd)] <- 0
.data2$x.sd[is.na(.data2$x.sd)] <- 0
.data2$x.sd <- apply(.data2$x.sd[prots.phospho,],2,as.numeric) + apply(.data2$x.sd[proteinlist.grep,],2,as.numeric)


dir.create(.setpath <- paste(gui.input$path.data,"/",foldername,"phosphopeptides-protein-reference",sep = "/")) 
setwd(.setpath)	

write.csv(ref.proteins.phospho.matrix,"Phosphopeptides-normalized-on-protein.csv")
write.csv(.data2$x.sd ,"Phosphopeptides-normalized-on-protein-sd.csv")

#.data2$x <- t(test)	
plot.type = 1



.data2$x <- t(test)
write.csv(.data2$x ,"log2-Phosphopeptides-normalized-on-protein.csv")


y.lab.input <- paste("log2( phospho-peptide", y.lab.input,"/","protein", y.lab.input,")")
	error.try <- class(.error <- try(hz.script.row.plot.space <- hz.script.row.plot(.data2,gui.input,y.lab.input, hz.cracker.anova.return$.aov.new,hz.exp.des.parse.data2,colorblind.set,.col,prog.max,ratio.prog,pb,ui)))


try(hz.script.heatmap(.data2,gui.input,prog.max,pb,ui,ratio.prog))
try(hz.script.heatmap2(.data2,gui.input,hz.cracker.anova.return$p.aov, hz.exp.des.parse.data2,.col,colorblind.set,prog.max,pb,ui, plot.type= plot.type,color.blind= color.blind, ratio.prog = ratio.prog))
#try(hz.script.hiercl)
try(hz.script.pca(.data2,gui.input, hz.exp.des.parse.data2,prog.max,pb,ui,ratio.prog))
try(hz.script.kmeans(.data2,gui.input,y.lab.input,colorblind.set,color.blind,.col,prog.max,pb,ui))
try(hz.script.graph(.data2,gui.input,prog.max,pb,ui))

					
					
}else{write.csv("No intersecting proteins available for comparissons.","readme.txt")}
}else{write.csv("No Phosphopeptides detectable.")}
}