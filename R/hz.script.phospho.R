hz.script.phospho <-
function(.data2,.data,gui.input, hz.exp.des.parse.data2,.col,.design,y.lab.input = hz.script.y.lab.return,prog.max,ratio.prog,pb,ui, plot.loop,path.data= gui.input$path.data, foldername,  colorblind.set, color.blind, hz.cracker.anova.return){
.data2 <- .data2  		

.data2$x <- .data2$x[grep(tolower(gui.input$phospho.string),tolower(rownames(.data2$x))),]
.data2$x.sd <- .data2$x.sd[grep(gui.input$phospho.string,rownames(.data2$x.sd)),]
.data2$prot.norm <- .data2$prot.n[grep(gui.input$phospho.string,rownames(.data2$prot.n)),]


dir.create(.setpath <- paste(path2,foldername,"phosphopeptides-all",sep = "/")) 
setwd(.setpath)	

	error.try <- class(.error<- try(hz.script.heatmap.return <- hz.script.heatmap()))
	error.try <- class( .error <- try(hz.script.pca.return <- hz.script.pca(.data2,gui.input)))
	error.try <- class(.error<- try(hz.script.heatmap2(.data2,gui.input,hz.cracker.anova.return$p.aov)))

	error.try <- class( .error <- try(hz.script.kmeans.return  <- hz.script.kmeans(.data2,gui.input,y.lab.input)))
	error.try <- class( .error <- try(hz.script.graph(.data2,gui.input)))
setwd("../")


p.split <- strsplit(rownames(.data2$x),"#")

phospho.matrix <- c()
peptidelist.grep <- c()
proteinlist.grep <- c()
sdlist.grep <- c()
for(o in 1:length(rownames(.data2$x))){
	temp.o <- p.split[[o]]
	phospho.matrix <- rbind(phospho.matrix,temp.o)
	
	peptidelist.grep <- c(peptidelist.grep,grep(paste(temp.o[1],temp.o[2],sep = "&"),rownames(.data2$peptidelist)))
	proteinlist.grep <- c(proteinlist.grep,grep(paste(temp.o[1]),rownames(.data2$x[-grep(gui.input$phospho.string,rownames(.data2$x)),])))
		
}

.data2$x[is.na(.data2$x)] <- 0
## norm on protein:

ref.protein <- .data2$x[-grep(gui.input$phospho.string,rownames(.data2$x)),][proteinlist.grep,]
ref.protein[is.na(ref.protein)] <- 0

ref.proteins <- intersect(phospho.matrix[,1],rownames(ref.protein))

ref.proteins.phospho.matrix <- c()
ref.proteins.phospho.matrix.row <- c()
prots.phospho <- c()
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


rownames(ref.proteins.phospho.matrix) <- ref.proteins.phospho.matrix.row
colnames(ref.proteins.phospho.matrix) <- colnames(.data2$x)

.data2$x <- ref.proteins.phospho.matrix

inf.matrix <- matrix("",dim(ref.proteins.phospho.matrix)[1],dim(ref.proteins.phospho.matrix)[2])
inf.matrix[is.infinite(ref.proteins.phospho.matrix)& ref.proteins.phospho.matrix > 0] <- "Inf"
inf.matrix[is.infinite(ref.proteins.phospho.matrix)& ref.proteins.phospho.matrix < 0] <- "-Inf"

.data2$x.sd[is.na(.data2$x.sd)] <- 0
.data2$x.sd[is.na(.data2$x.sd)] <- 0
.data2$x.sd <- apply(.data2$x.sd[prots.phospho,],2,as.numeric) + apply(.data2$x.sd[proteinlist.grep,],2,as.numeric)


dir.create(.setpath <- paste(path2,foldername,"phosphopeptides-protein-reference",sep = "/")) 
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

					
					
}
