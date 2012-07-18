hz.script.volcano <-
function(.data2,gui.input,extented.info, colorblind.set,color.blind,hz.cracker.anova.return,prog.max,pb,ui){
	
	
	if(!exists("ratio.prog")){ratio.prog <- 1000}
	
anov <- 	hz.cracker.anova.return$p.aov
ttestlist <- hz.cracker.anova.return$ttestlist
data <- 	.data2$x
ratio.thres <- as.numeric(gui.input$ratio.thres)
#info <- read.csv("/Users/cRacker-DEMO/henrik/SRE-sugar-fractions/cracker/ms-analysis-IonIntensity--exp-csv-2011-10-14/INFO-proteinlist.csv")

choose.vec <- 1:dim(data)[2]
comb.vec	<- combn(choose.vec,2)

volcano.output <- c()
ttestlist[,1] <- make.names(ttestlist[,1])
ttestlist[,2] <- make.names(ttestlist[,2])

if(gui.input$graphic.type == "pdf"){
	pdf("volcano-plot.pdf",width = 7)
	par(mai = c(1,1,0.2,0.2),mfrow= c(1,1))

}

.wd <- getwd()

for(i in 1:dim(comb.vec)[2]){

	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb, ratio.prog*7+ratio.prog/dim(comb.vec)[2]*i, label=paste( "Volcanoplot..."))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, ratio.prog/dim(comb.vec)[2]*i, label=paste( "Volcanoplot..."))))
		}
	##############
	
list.1 <- colnames(data)[comb.vec[1,i]]
list.2 <- colnames(data)[comb.vec[2,i]]
list.1 <- gsub(" ","",tolower(list.1))
list.2 <- gsub(" ","",tolower(list.2))

#ttestlist[,1] <- make.names(ttestlist[,1])
#ttestlist[,2] <- make.names(ttestlist[,2])


test <- cbind(c(ttestlist[,1] == list.1|ttestlist[,1] == list.2),c(ttestlist[,2] == list.1|ttestlist[,2] == list.2))
test2 <- apply(test,1,all)
test2[is.na(test2)] <- FALSE

if(!all(!test2)){
	
filtered <- ttestlist[test2,]
if(is.vector(filtered)){
	filtered <- t(as.matrix(filtered))
}

control <- hz.merge.control(gsub(" ","",filtered[,4]),gsub(" ","",rownames(data)))
p.val <- filtered[control,]
#p.val <- cbind(p.val,p.adjust(p.val[,3],method = gui.input$p.adjust.method))

data.log 	<- log2(as.numeric(data[,comb.vec[1,i]])/as.numeric(data[,comb.vec[2,i]]))
error.try <- class(try(test 	<- grep(paste(rownames(data),collapse = "|"), anov[,1])))

if(error.try == "try-error"){
	
	test <- c()
	for(ri in 1:length(rownames(data))){
		test <- c(test,grep(gsub(" ","",rownames(data)[ri]),anov[,1]))
		
	}
	
}

test2 	<- anov[test ,]
if(is.vector(test2)){test2 <- t(as.matrix(test2))}
test2 	<- test2[order(test2[,2]),]

te <- rep(1,length(data.log))

col.set <- c("pink","lightblue","2","4")

if(colorblind.set){
	#col.set <-unlist(color.blind)[c(2,4,5,7)]
	
}
p.val[is.na(p.val)] <- 1


te[data.log >	ratio.thres	& as.numeric(p.val[,3]) < gui.input$p.value] 		<- col.set[1]
te[data.log < 	-ratio.thres	& as.numeric(p.val[,3]) < gui.input$p.value] 	<- col.set[2]
te[data.log >	ratio.thres	& as.numeric(p.val[,5]) < gui.input$p.value] 		<- col.set[3]
te[data.log < 	-ratio.thres	& as.numeric(p.val[,5]) < gui.input$p.value] 	<- col.set[4]

counts <- aggregate(te,list(te),length)
test <- apply(cbind(data.log,-log10(as.numeric(p.val[,3]))),1,function(x){
	return(!any(is.na(x)|is.infinite(x)))
	
})


if(gui.input$graphic.type == "eps"){
	dir.create(.wd.set <- "volcano-eps")
	ps.name <- paste("volcano-plot",colnames(data)[comb.vec[1,i]],colnames(data)[comb.vec[2,i]],".eps",sep = "-")
	 
	postscript(paste(.wd, .wd.set,ps.name,sep = "/"),width = 7,height = 7, paper = "special",onefile = FALSE,horizontal = FALSE)

	par(mai = c(1,1,0.2,0.2),mfrow= c(1,1))
	
}



plot(data.log,-log10(as.numeric(p.val[,3])),col = te,xlab = paste("log2(",colnames(data)[comb.vec[1,i]],"/",colnames(data)[comb.vec[2,i]],")"),ylab = "-log10(p value)",lty = 3,pch = 1,lwd = 2)

abline(h=-log10(gui.input$p.value), v = c(ratio.thres,-ratio.thres),col = "grey",pch = 16)
leg.vec <- c(
				paste("p >",gui.input$p.value			,colnames(data)[comb.vec[1,i]],counts[counts[,1]==col.set[1],2]),
				paste("p >",gui.input$p.value			,colnames(data)[comb.vec[2,i]],counts[counts[,1]==col.set[2],2]),
				paste("adjusted p >",gui.input$p.value	,colnames(data)[comb.vec[1,i]],counts[counts[,1]==col.set[3],2]),
				paste("adjusted p >",gui.input$p.value	,colnames(data)[comb.vec[2,i]],counts[counts[,1]==col.set[4],2])
				

				)
legend("topright",leg.vec,col = col.set,pch = 1,bg = "#ffffff90",bty  = "o",pt.lwd = 2,cex = 0.7,title = paste(length(test[test == TRUE]), "data points"))

temp.data <- cbind(paste(colnames(data)[comb.vec[,i]],collapse ="/"),rownames(data),data.log,p.val[,3],p.val[,5])

}else{
	print("no compar")
}
if(exists("temp.data")){
	volcano.output  <- rbind(volcano.output,temp.data) 
}

}
graphics.off()

#write.csv(ttestlist,"")
if(length(volcano.output) != 0){

colnames(volcano.output)	<-	c("samples","accession","ratio","ttest.p.value",paste("p.value",gui.input$p.adjust.method,"corrected",sep = ".")) 
#test <- merge(signi.double,info[,c(1,4)],by = 1)
if(!exists("extended.info")){
extended.info 		<- .data2$proteinlist.info
volcano.info.add 	<- hz.merge.control(as.character(extended.info[,colnames(extended.info) == "code"]),tolower(volcano.output[,2]))
}
volcano.info.add <- hz.merge.control(as.character(extended.info[,colnames(extended.info) == "code"]),tolower(volcano.output[,2]))


if(length(volcano.info.add) == dim(volcano.output)[1]){
	volcano.output  <- cbind(volcano.output , extended.info[volcano.info.add,-1])
}

write.csv(volcano.output,file = "volcano-data.csv") 
signi.volcano.output <- volcano.output[(as.numeric(volcano.output[,3]) < -ratio.thres | as.numeric(volcano.output[,3]) > ratio.thres ) & as.numeric(volcano.output[,4]) < gui.input$p.value	
,]

na.row.exclude <- apply(signi.volcano.output,1,function(x){all(is.na(x))})
print(dim(signi.volcano.output[!na.row.exclude,]))
write.csv(signi.volcano.output[!na.row.exclude,],"volcano-signi-data-p-value-uncorrected.csv")
unlink("ttest-pvalues.csv")
}

setwd(.wd)
}
