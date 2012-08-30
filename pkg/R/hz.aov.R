hz.aov <-
function(x = aov.export,.data2,gui.input, progressbar = TRUE,prog.max,pb,ui){
	count <- 1
	if(!exists("ratio.prog")){ratio.prog <- 1000}

	#library(gridExtra)	
#p.adjust.method = "BH"
rm(input.all.list)
x 		<- as.data.frame(x, stringsAsFactors = FALSE)

if(!is.na(max(as.numeric(x[,2],na.rm = T)))){
	print("re")
x[,2] <- colnames(.data2$x)[as.numeric(x[,2])]
	
	
}


.acc 	<- x[,1]
uni.acc	<- sort(unique(.acc))
total 	<- length(uni.acc)
p.all	<- c()
pt		<- list()
list.ttest 	<- list()
list.ttest.one.sided <- list()
length.data <- 7#length(unique(x[,2]))
if(length.data>6){
	length.data <- 7 + 0.85*(length.data-6)
	
}

#stop()
#try(pdf("pairwise-ttest-tables.pdf",width = length.data,height = length.data) )## only for ttest

for(i in 1:length(uni.acc)){
#			if(uni.acc[i]== "at1g01470.1 "){stop()}

	#if(i == 20){stop()}
	
if(progressbar ){

	##############	GUI
	ratio.prog2 <- (as.numeric(10000)/8)/as.numeric(total)
	if(exists("ratio.prog.2")){
		if(is.na(ratio.prog.2)){ratio.prog.2 <- 5000}
	}else{ratio.prog.2 <- 5000}
	


	pb.check	<- class(try(ui$setProgressBar(pb, i*ratio.prog2, label=paste(round(i/length(uni.acc)*100),  "% Anova/ttest done"))))
	pb.check	<- class(try(ui$setProgressBar(pb, i*ratio.prog2, label=paste(round(i/length(uni.acc)*100),  "% Anova/ttest done"))))

	while(pb.check == "try-error"){
			print("Warning: User closed window!")
	
				
		pb 			<- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)

		pb.check	<- class(try(ui$setProgressBar(pb, i*ratio.prog2, label=paste(round(i/length(uni.acc)*100),  "% done"))))
	}
	##############
	
	
}		
	temp.i 				<- x[.acc == uni.acc[i],c(2,3)] 
	if(is.vector(temp.i) == FALSE){
		.aov.true <- TRUE
		rownames(temp.i) 	<- 1:dim(temp.i)[1]
	}else{
		.aov.true <- FALSE
		temp.i 				<- t(as.data.frame(temp.i))
		#print(dim(temp.i))
		rownames(temp.i) 	<- 1
	}

	temp.i.aov <- temp.i
	temp.i.aov[,2] <- as.numeric(	temp.i.aov[,2])
	temp.i.aov	<- as.data.frame(temp.i.aov, stringsAsFactors = FALSE)


if(length(unique(temp.i.aov$experiment)) != length(temp.i.aov$experiment)& .aov.true& length(unique(temp.i.aov$experiment)) > 1){
	
	temp.aov 	<- aov(as.numeric(temp.i.aov$intensity)~ as.factor(temp.i.aov$experiment))


	temp.agg <- aggregate(temp.i.aov[,1],list(temp.i.aov[,1]),FUN=length)
			print(temp.agg)

	
	
	temp.agg <- temp.agg[temp.agg[,2] < 2,1]
	oneside.ttest <- c()
	if(length(temp.agg) !=0 ){
		print("oneside-ttest")
	oneside.ttest <- hz.oneside.ttest(temp.i.aov)
	temp.exclu	<- grep(paste(temp.agg,collapse = "|"),temp.i.aov[,1])
	temp.i.aov 		<- temp.i.aov[-temp.exclu,]
	}
	
	

	
	temp.sum 		<- summary(temp.aov)
	p.value 		<- temp.sum[1][[1]][[5]][1]}else{p.value <- 1}
	if(is.na(p.value[1])){p.value <- 1	}
	
	
	p.all 			<- c(p.all,p.value[1])
	
	if(p.value <= .data2$gui.input$p.value | 1 == 1){

#try(	temp.i.aov <- temp.i.aov[-c(1,2,3),])
	#temp.i.aov <- temp[-9,]

	print(temp.i.aov)
	type.test <- aggregate(temp.i.aov[,1],list(temp.i.aov[,1]),length)
	
		assign("temp.type.test",type.test,envir = .GlobalEnv)

	assign("temp",temp.i.aov,envir = .GlobalEnv)
	
	
	temp.pt	<- pairwise.t.test(as.numeric(temp.i.aov[,2]), temp.i.aov[,1],p.adjust.method="none")
	

	
	
	


	list.ttest[[i]]		<- temp.pt
	input				<- temp.pt$p.value
	
	


	input.vec 	<- as.numeric(input)
	input.vec.1 <- rep(rownames(input),dim(input)[2])

	input.vec.2 <- lapply(colnames(input),function(x){rep(x,dim(input)[1])})

	input.list 	<- cbind(input.vec.1,unlist(input.vec.2),input.vec)

assign("input",input,envir= .GlobalEnv)
assign("input.list",input.list,envir= .GlobalEnv)

	input.list 	<- input.list[!is.na(input.list[,3]),]
	include.oneside.ttest<- T
	if(exists("list.ttest.one.sided")& include.oneside.ttest){
		#list.ttest.one.sided[[i]] <- oneside.ttest$test
		
		if(is.matrix(input.list)){rep.factor <- dim(input.list)[1]}else{rep.factor <- 1}
		if(is.matrix(oneside.ttest$one.side.test.results)){
			rep.factor2 <- dim(oneside.ttest$one.side.test.results)[1]
		}else{
			#rep.factor<- t(as.matrix(oneside.ttest$one.side.test.results))
			rep.factor2 <- 1
		}
		
		type.vector <- c(rep("two.sided",rep.factor),rep("one.sided",rep.factor2))

		input.list <- rbind(input.list,oneside.ttest$one.side.test.results[,c(3,2,1)])
		
	}else{
		if(is.matrix(dim(input.list)[1])){rep.factor <- dim(input.list)[1]}else{rep.factor <- 1}
		type.vector <- c(rep("two.sided",rep.factor))
	}
	
	if(is.vector(input.list)){
		input.list <- t(as.matrix(input.list))
		
	}
	input.list	<- cbind(input.list,uni.acc[i])
		#print(is.matrix(input.list))
	if(dim(input.list)[2] != 4){
		#print(input.list)

	}
	if(!exists("input.all.list")){
		type.all.vector <- type.vector
		input.all.list <- input.list
		#print(dim(input.all.list))
	}else{
	
		type.all.vector <- c(type.all.vector,type.vector)
		input.all.list	<- rbind(input.all.list, input.list)
	}

	write.t.tables <- F
	if(write.t.tables){
		if(count == 1){
			hz.print.des <- function(){
				#library(plotrix)
				text.input.title <- "Pairwise t-test"

				text.input1 <- "Selection: Pairwise ttests are applied to proteins, showing a significant difference between mean values using ANOVA. For this selection the uncorrected p values from ANOVA will be used."

				text.input2<- "Pairwise t-test: This test tests protein intensities, calculated from at least 2 peptide intensities. For each protein (not for all protein comparissons together!) the calculated p values are already adjusted according to the set method.The different samples that are tested are given by alias numbers (see next page)."

				plot.new()
				textbox(c(0,1),1, text.input.title)
				textbox(c(0,1),0.9,text.input1,box = FALSE)
				textbox(c(0,1),0.7,text.input2,box = FALSE)
				grid.newpage()
			}
			try(hz.print.des())
			
		if(.data2$gui.input$raw){
			run.type <- 1		
		}else{
			run.type <- 2
		}
		exp.des <- cbind(.data2$exp.design[, run.type],1:length(.data2$exp.design[, run.type]))
		exp.des <- rbind(c("name","alias"),exp.des)
		print(exp.des)
			count <- count+1
			alias <- 1/(dim(exp.des)[1]+2) 
			alias <- seq(from = alias, to = alias*dim(exp.des)[1] , by = alias)
			alias <- rep(alias,(dim(exp.des)[2]))
			grid.text(rev(exp.des),c(rep(0.2,dim(exp.des)[1]),rep(0.4,dim(exp.des)[1])),alias,just = "left")
		
	}
	
#	print(i)
#	print(exp.des)
		hz.write.table(temp.pt,temp.i.aov = temp.i.aov,name = uni.acc[i],acc = .acc,exp.des = exp.des)
		#stop()
	}
#	print(i)
	}else{
	list.ttest[[i]] <- paste(uni.acc[i],": Not above significance threshold in Anova. Pairwise ttest was not applied.")
	}
}


#dev.off()

try(names(list.ttest) <- uni.acc)

p.all.adjust <- p.adjust(p.all,method = gui.input$p.adjust.method)
.return <- cbind(uni.acc,p.all,p.all.adjust)
colnames(.return) <- c( "accession","aov p",paste("aov p",gui.input$p.adjust.method))

if(1 == 0){

ttest.results <- file("pairwise-ttest.txt", open="wt")
sink(ttest.results)
sink(ttest.results,type = "message")
print(list.ttest)
sink(type = "message")
sink()
sink(console,append = TRUE)
sink(console, append = TRUE, type="message")
}
colnames(input.all.list) <- c("sample.1","sample.2","p.value","protein")
write.csv(input.all.list,"ttest-pvalues.csv")
return(list(aov=.return,pt = pt,ttest = list.ttest,ttestlist = input.all.list,type.vector.ttest = type.all.vector))
}
