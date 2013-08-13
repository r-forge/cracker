hz.norm <-
function(	x,					# matrix, should contain numeric values
							margin = 1,		# 1 for rows, 2 for columns
							na.sub = NA,		# Sub for Nas in data set
							norm = "mean",	# "mean" for mean, "median" for median, "sum fpr sum
							saliva = FALSE,
							gender = "W",
							group = NULL,
							template = ""
							){
			

			
			if(is.vector(x)	== TRUE){										col		<- names(x)
							x 		<- t(as.matrix(as.numeric(x)))
							row		<- "t"}else{
			if(dim(x)[1]==1){	row <- rownames(x);
								col <- colnames(x);
								x	<- t(as.matrix(as.vector(x)))}else{	
	row <- rownames(x)
	col <- colnames(x)
	#==============
	# As.numeric, NA substitution
	
	
	x <- as.matrix(x)
	x <- apply(x,2,as.numeric)
	x[sapply(x, is.na)] <- na.sub}}
	x[is.infinite(x)] <- NA
	#==============
	# mean or median normalisation:
	
	if(saliva){
		print("Applying Saliva normalisation!")
		if(gender == "W"){gender = 12}else{gender = 3}
		
			apply(x,1,function(x){
			mean.temp 	<- 	temp.k[c(1:gender)]
	 		mean.temp	<- 	mean.temp[mean.temp != 0]
			mean.rest	<- mean(mean.temp,na.rm=TRUE)
				if(length(mean.rest[!is.na(mean.rest)]) == 0){
					mean.rest <- mean(x,na.rm = TRUE)}
			x	<- x/mean.rest		
					})
		
		}else{
	all.data <- c()		
	
	backup <- x
	
	if(is.null(group)){
		if(margin ==1){
			group <- rep(1,dim(x)[2])

		}
		if(margin ==2){
			group <- rep(1,dim(x)[1])

		}
	}

	for(group.i in 1:length(unique(group))){
		temp.group <- group ==unique(group)[group.i]

		if(margin == 1){
			if(length(group) != dim(backup)[2]){print("ALARM");stop()}
		x <- backup[,temp.group]
		}
		if(margin == 2){
			if(length(group) != dim(backup)[1]){print("ALARM");stop()}

		x <- backup[ temp.group,]
		}

if(is.vector(x)&margin == 1){x <- (as.matrix(x))}
if(is.vector(x)&margin == 2){x <- t(as.matrix(x))}

	if(norm == "sum"){
		
	if(template != ""){
		target.x <- x[-grep(template,row),]
	}else{
		target.x <- x
	}
		

	x.sum <- apply(target.x,margin,function(a){sum(a,na.rm = TRUE)})
		if(margin == 2){
			x <- t(t(x)/x.sum)
		}
		
		if(margin == 1){x <- x/x.sum}	
	}
	
	if(norm == "sum50"){
		
		
		
		
		if(margin == 2){sd.margin <-1}
		if(margin == 1){sd.margin <-2}
		

	if(template != ""){
		target.x <- x[-grep(template,row),]
	}else{
		target.x <- x
	}
	
	row.sum <- apply(target.x,1,function(x){sum(x,na.rm = T)})
	
	target.x <- target.x[!row.sum == 0,]
	
	
	.sd <- 	apply(as.matrix(target.x),sd.margin,function(a){a<- sd(a,na.rm = TRUE);return(a)})

	x.sum <- apply(target.x[.sd <=median(.sd),],margin,function(a){sum(a,na.rm = T)})
	.sum.all <- apply(x,margin,function(a){sum(a,na.rm = T)})
	if(any(.sum.all ==0)){
	print(paste("Warning, the following samples are normalized on sum of total, not enough proteins below 50% sd.\n",names(x.sum)[x.sum == 0]))
	}
		x.sum[x.sum ==0] <- .sum.all[x.sum == 0]

		if(margin == 2){
			x <- t(t(x)/x.sum)
		}
		
		if(margin == 1){x <- x/x.sum}	
	}
	
	if(norm == "equal.peptides"){
		
		if(margin ==1){margin.opp <- 2;.byrow = F}else{margin.opp <- 1;.byrow = T}
		x.num 		<- apply(as.matrix(x),2,as.numeric)
		x.num.logic <- apply(x,margin.opp,function(x){any(is.na(x)|x==0)})
		x.num.sum 	<- x.num[!x.num.logic,]
		dim.x.num.sum 	<- dim(x.num.sum)
		if(!is.vector(x.num.sum)){
			x.num.sum 		<- apply(x.num.sum,margin,sum)
		}
		
		x.num.sum.m <- matrix(x.num.sum,ncol = dim(x)[2],nrow = dim(x)[1],byrow = .byrow)
		x <- x/x.num.sum.m
		
		
	}
	
	if(norm == "median"){
		
		x.median <- apply(x,margin,function(a){median(a,na.rm = TRUE)})
		if(margin == 2){
			x <- t(t(x)/x.median)
		}
		
		if(margin == 1){x <- x/x.median}
		
		}
		
			if(norm == "mean"){
	x.mean <- apply(x,margin,function(a){mean(a,na.rm = TRUE)})
		if(margin == 2){
			x <- t(t(x)/x.mean)
		}
		
		if(margin == 1){x <- x/x.mean}	
	}
	
	if(norm == "z"){
		x.mean 	<- apply(x,margin,function(a){mean(a,na.rm = TRUE)})
		x.sd 			<- apply(x,margin,function(a){sd(a,na.rm = TRUE)})
		x.sd[is.na(x.sd)]			<- 0
		if(margin == 2){
			x <- t(t(x-x.mean)/x.sd)
		}
		
		if(margin == 1){x <- (x-x.mean)/x.sd}	
		
				

		
	}
	
	
		if(norm == "z-pos"){
	x.mean <- apply(x,margin,function(a){mean(a,na.rm = TRUE)})
	x.sd <- apply(x,margin,function(a){sd(a,na.rm = TRUE)})

		if(margin == 2){
			x <- t(t(x-x.mean)/x.sd)
		}
		
		if(margin == 1){x <- (x-x.mean)/x.sd}	
		
				x.min <- min(x,na.rm = T)
				x <- x+ abs(x.min)+abs(x.min)*0.00001
		
	}

	if(margin == 1){
		
	backup[, temp.group] <- x	

	}
	if(margin == 2){
	backup[temp.group,] <- x	
	}
	
	#print(temp.group)
	
	}
	
	}
	x <- backup
	
	colnames(x) <- col
	rownames(x) <- row
	if(!exists("dim.x.num.sum")){
		dim.x.num.sum <- dim(x)
	}
	return(list(x = x,sum = apply(x,margin,function(a){sum(a,na.rm = TRUE)}),method = paste("margin:",margin,"| na.sub:",na.sub,"| Normalisation:",norm),dim.norm.m = dim.x.num.sum))
	#return(apply(x,margin,function(a){sum(a,na.rm = TRUE)}))
	}
