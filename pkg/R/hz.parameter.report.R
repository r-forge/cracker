hz.parameter.report <-
function(gui.input,.data2,.data,.report,.design,foldername){
report.pca.bpca <- .report$report.pca.bpca	
report.heatmap.norm <- .report$report.heatmap.norm
.norm <- ""
if(gui.input$calc.empai){.norm <- paste(.norm,"emPAI\n",sep = "")
}else{
.norm <- paste("ion intensities\n",.norm,sep = "")
}
if(gui.input$quant.method == "lf"){.norm <- paste(.norm,"label free\n",sep = "")}
if(gui.input$quant.method == "cbn"){.norm <- paste(.norm,"reference protein normalization on",gui.input$cbn.prot,"\n",sep = ";")}
if(gui.input$quant.method == "15n"){
	.norm <- paste(.norm,"labeled normalization\n",sep = "")
	
	
	}
if(gui.input$n15.log2){
	.norm <- paste(.norm,"log2 transformation\n",sep = "")
	}
	if(gui.input$n15.correct.method != "none"){
	.norm <- paste(.norm," ",gui.input$n15.correct.method," Correction of ratios to ratio ",gui.input$n15.correct.expect,"\n",sep = "")
	}	
	


if(gui.input$n.correction){.norm <- paste(.norm,"n corrected",sep = ";")}

cat(.norm)
# scaling
.scale <- "none"
if(gui.input$row.norm){
	.scale <- gui.input$norm.method
	.scale <- paste(.scale,"; peptides with NA > ",gui.input$shape*100,"% occurence in samples are excluded ;",sep = "")

}

# general info method
.info <- ""
.info <- paste(.info,"Averaging method, peptides -> protein: ", gui.input$peptidemerge, "\n",sep = "")
.info <- paste(.info,"Duplicated peptide species, have been treated with: ", gui.input$dupli.val,"\n", sep = "")
.info <- paste(.info,"Excluded peptides with score < ", gui.input$score,"\n", sep = "")
if(gui.input$exclu){
.info <- paste(.info,"Excluded redundant peptides\n",sep = "")
}
if(gui.input$phospho){
.info <- paste(.info,"Run phospho peptides separately.\n",sep = "")
}

if(gui.input$outlier!= "NA"){
	if(gui.input$outlier == "top.3"){
		.info <- paste(.info,"Used top 3 intensity peptides.\n",sep = "")
	}
	if(gui.input$outlier != "top.3"){
		.info <- paste(.info,"Used outlier exclusion.\n",sep = "")
	}
	
}

.info <- paste(.info,"Removed 0 values with ",gui.input$zero.treat,"\n",sep = "")

#### statistics
.statistic <- c("kmeans","cluster",gui.input$centers)
.statistic <- rbind(.statistic,c("general","p value",gui.input$p.value))
.statistic <- rbind(.statistic,c("general","p value correction method",gui.input$p.adjust.method))
.statistic <- rbind(.statistic,c("mapping","file",gui.input$go.library))
.statistic <- rbind(.statistic,c("correlation analysis","method",gui.input$do.cor))
.statistic <- rbind(.statistic,c("heatmap","matrix scaling method", report.heatmap.norm))
.statistic <- rbind(.statistic,c("heatmap","data transformation","log2"))
.statistic <- rbind(.statistic,c("heatmap","distance method", "euclidean"))
.statistic <- rbind(.statistic,c("heatmap","agglomeration method (hclust)", "average"))
.statistic <- rbind(.statistic,c("pca", "data transformation", "log2"))
.statistic <- rbind(.statistic,c("pca", "NA replacement",report.pca.bpca ))



.statistic <- rbind(c("Process","option","used"),.statistic)
rownames(.statistic) <- NULL
#.statistic <- as.data.frame(.statistic)
#####
#pca
#heatmap
#hierarchical clustering
#kmeans clustering


parameters.write <- function(){
zz <- file(paste(gui.input$path.data,"/",foldername,"/report.txt",sep = "/"), open="wt")
sink(zz)
sink(zz, type="message")
title.c <- paste("\n# Report cRacker run ",as.character(Sys.time())," #\n",sep = "")
cat("\n",rep("#",nchar(title.c)-2),title.c,rep("#",nchar(title.c)-2 ),"\n",sep = "")
title.c <- paste("\n# sample normalization "," #\n",sep = "")
cat("\n",rep("#",nchar(title.c)-2),title.c,rep("#",nchar(title.c)-2 ),"\n",sep = "")


cat(.norm)

title.c <- paste("\n# general settings/informations "," #\n",sep = "")
cat("\n",rep("#",nchar(title.c)-2),title.c,rep("#",nchar(title.c)-2 ),"\n",sep = "")

cat(.info)

title.c <- paste("\n# statistics settings"," #\n",sep = "")
cat("\n",rep("#",nchar(title.c)-2),title.c,rep("#",nchar(title.c)-2 ),"\n",sep = "")

cat.matrix <- function(.statistic){max.length <-	apply(.statistic,1,function(x){max(nchar(x))})
max.col 	<- apply(.statistic,2,function(x){max(nchar(x))})

for(i in 1:(dim(.statistic)[1])){
	temp <- .statistic[i,]
	
	temp.i <- c()
	for(a in 1:length(temp)){
		
		temp.i <- c(temp.i,paste(as.character(temp[a]),paste(rep(" ",max.col[a]-nchar(temp[a])),collapse = ""),collapse = ""))
		
	}
	
	cat( tempe<- paste("",paste(temp.i,collapse = "\t"),"\n") )
	if(i == 1){
		cat(rep("-",nchar(tempe)),"\n",sep = "")
		
	}
	
	
	
}}

cat.matrix(.statistic)

title.c <- paste("\n# dataset information "," #\n",sep = "")
cat("\n",rep("#",nchar(title.c)-2),title.c,rep("#",nchar(title.c)-2 ),"\n",sep = "")

protein.total <- length(unique(.data$code))
peptide.total <- length(unique(.data$sequence))

analysed.proteins <- dim(.data2$x)[1]
analysed.peptides <- rownames(.data2$used.peptides)


temp.ana<- unique(cbind(.data$code,.data$rawfilename))
proteins.p.sample <- aggregate(temp.ana[,1],list(temp.ana[,2]),length)

temp.ana<- unique(cbind(.data$sequence,.data$rawfilename))
sequences.p.sample <- aggregate(temp.ana[,1],list(temp.ana[,2]),length)

.datar <- .data[as.numeric(.data$Score) > gui.input$score,]
.score.print <- F
if(dim(.datar)[1] != dim(.data)[1]){
	print("score")
.score.print <- T
	
protein.total.score <- length(unique(.datar$code))
peptide.total.score <- length(unique(.datar$sequence))

temp.ana<- unique(cbind(.datar$code,.datar$rawfilename))
proteins.p.sample.score <- aggregate(temp.ana[,1],list(temp.ana[,2]),length)

temp.ana			<- unique(cbind(.datar$sequence,.datar$rawfilename))
sequences.p.sample.score <- aggregate(temp.ana[,1],list(temp.ana[,2]),length)

}

#### statistics
data.info <- c("info:","","")
data.info <- rbind(data.info,c("total number of protein in experiment","", protein.total))
data.info <- rbind(data.info,c("total number of peptide species in experiment","",peptide.total))
if(.score.print){
data.info <- rbind(data.info,c("total number of protein n in experiment","score corrected",protein.total.score))
data.info <- rbind(data.info,c("total number of peptide species  n in experiment","score corrected", peptide.total.score))
colnames(proteins.p.sample.score) <- c("sample","n protein")
colnames(sequences.p.sample.score) <- c("sample","n peptide species")




}

proteins.p.sample <- rbind(c("sample","n protein"),proteins.p.sample)
sequences.p.sample <- rbind(c("sample","n peptides"),sequences.p.sample)



cat("\n\nProteins/peptide species used in experiment\n\n")
cat.matrix(data.info)

cat("\n\nProteins/peptide species used in experiment\n\n")
cat.matrix(proteins.p.sample)

cat("\n\nnumber of proteins/peptides, per sample:\n\n")
cat.matrix(proteins.p.sample)
cat.matrix(sequences.p.sample)
if(.score.print){
cat("\n\nnumber of proteins/peptides above Score threshold, per sample:\n\n")

cat.matrix(proteins.p.sample.score)
cat.matrix(sequences.p.sample.score)	
	
}





sink(type = "message")
sink()
}
parameters.write()

}
