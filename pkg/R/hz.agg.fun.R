hz.agg.fun <- 
function(
	x,outlier,row.norm,merge.method,ui,pb,prog.max,pb.string
){
									
				temp.i <- 			unlist(x)
				print(length(temp.i))
				if(outlier == "row"& row.norm == TRUE) {
											
						box.temp.t	<- boxplot.stats(as.numeric(temp.i))$out
						for(r in 1:length(box.temp.t)) {
							temp.r <- box.temp.t[r]
							temp.i[temp.i == temp.r] <- NA
						}
						
						if(length(box.temp.t) != 0) {
							box.ex <-1
						} else {
							box.ex <- NA
						}
								
				}
				
				if(outlier == "all.below" &  row.norm == FALSE) {
					temp.t		<- temp.i
					box.temp.t	<- boxplot.stats(as.numeric(temp.t))$out
					box.temp.t	<- box.temp.t[box.temp.t < median(as.numeric(temp.i),na.rm = TRUE)]
					for(r in 1:length(box.temp.t)) {
						temp.r <- box.temp.t[r]
						temp.t[temp.t == temp.r] <- NA
					}
					if(length(box.temp.t) != 0){
						box.ex <- box.ex + 1
					} else {
						box.ex <- NA
					}
					temp.i <- temp.t
				}
				
				if(outlier == "top.3"){
				temp.i.test 	<- temp.i[order(temp.i,decreasing = TRUE)]
				temp.i.test  	<- temp.i.test[!is.na(temp.i.test)]
				temp.i.test 	<- temp.i.test[1:3]

				temp.i[setdiff(1:length(temp.i),grep(paste(temp.i.test,collapse = "|"),temp.i))] <- NA

				}
				
				###### module for bsa test absolute values
				temp.3 			<- temp.i
				temp.3.test 	<- temp.3[order(temp.3,decreasing = TRUE)]
				temp.3.test  	<- temp.3.test[!is.na(temp.3.test)]
				temp.3.test 	<- temp.3.test[1:3]
				temp.3[setdiff(1:length(temp.3),grep(paste(temp.3.test,collapse = "|"),temp.3))] <- NA
				if(merge.method!="sum"){
				if(merge.method == "mean"){
							temp.3.mean <- mean(as.numeric(as.matrix(temp.3)),na.rm = TRUE)}
if(merge.method == "median"){
							temp.3.mean <- median(as.numeric(as.matrix(temp.3)),na.rm = TRUE)}
				}else{
					temp.3.mean <- NA
				}
				
							#print(temp.i)
							aov.data.2  <- paste(as.character(temp.i),collapse = "#")
							temp.o <- temp.i 
							#temp.aov  <- temp.o[!is.na(temp.o)]
							#temp.aov2 <- rep(as.character(x[1]),length(temp.aov))
							#temp.aov3 <- rep(experiment[i],length(temp.aov))
							#aov.data.2 <- cbind(temp.aov2,temp.aov3,temp.aov)
							assign("cRacker.temp.aov.data.2",aov.data.2,envir = .GlobalEnv)
							
							assign("temp.o",temp.o,envir = .GlobalEnv)
#print(temp.o[!is.na(temp.o)])
							if(merge.method == "mean"){
							temp.o.mean <- mean(as.numeric(as.matrix(temp.o)),na.rm = TRUE)}
							
							if(merge.method == "sum"){
							temp.o.mean <- sum(as.numeric(as.matrix(temp.o)),na.rm = TRUE)}
							
							if(merge.method == "median"){
							temp.o.mean <- median(as.numeric(as.matrix(temp.o)),na.rm = TRUE)}
				
							if(merge.method != "sum"){
								temp.o.sd		<- sd(as.numeric(as.matrix(temp.o)),na.rm = TRUE)
							}else{temp.o.sd <- NA}
							
							temp.o.count 	<- length(temp.o[!is.na(temp.o)|temp.o == 0]) 
							
							
							temp.o.sd.rel	<- temp.o.sd/temp.o.mean 
							
							number.o		<- length(as.numeric(temp.o)[!is.na(temp.o)])
							
							data <- c(as.character(x[1]),temp.o.mean,temp.o.sd,temp.o.sd.rel,number.o,temp.3.mean,aov.data.2)
							names(data) <- c("Acc","mean","sd","rel.sd","n","top3","aov.data")
							assign("cracker.counter.temp",cracker.counter.temp+1,envir = .GlobalEnv)
	########## GUI 
	pb.check <- class(try(ui$setProgressBar(pb,cracker.ratio.prog.temp*cracker.counter.temp, label=paste("Averaging Peptides in sample",pb.string[1],"/",pb.string[2]))))

		while(pb.check == "try-error"){
				print("Warning: User closed window!")
				pb <- ui$progressBar(title = "cRacker", min = 0,max = prog.max, width = 300)
				pb.check <- class(try(ui$setProgressBar(pb, cracker.ratio.prog.temp*i, label=paste("Averaging Peptides:",cracker.counter.temp))))
		}
	##############
							#data <- list(data = data,aov= aov.data.2)
							return(data)
}

	
	
			