hz.shape <-
function(	x,
						shape = 0.5,
						group.shape = NA
					){
	
x	<- as.matrix(x)
x.shape.all <- c()
n.vec.all <- c()
x.na.true.all <- c()
#shape <- 0.3
if(!all(is.na(!is.na(group.shape)))){
	if(length(group.shape)!= dim(x)[2]){
	group.shape <- rep(1,dim(x)[2])	
	}
}

	x.na  <- x
	x.na[x.na	== "NaN"] <- NA

uni.group.shape <- unique(group.shape)
exclude.matrix <- c()
for(i in 1:length(uni.group.shape)){
	
	if(!is.na(uni.group.shape[i])){
	print(paste("Filtering group",i))
	sub.x <- x.na[,group.shape == uni.group.shape[i]]
	if(is.vector(sub.x)){
		sub.x <- (as.matrix(sub.x))
	}

	sub.x[sub.x== "NaN"] <- NA
	sub.x[sub.x!= "NA"] <- 1
	sub.x[is.na(sub.x)] <- 0

	sub.x.sum <- apply(sub.x,1,sum)
	sub.x.true <- sub.x.sum >=  dim(sub.x)[2]* shape

	
exclude.matrix <- cbind(exclude.matrix, sub.x.true)
	}
}
x.na.true <- apply(exclude.matrix,1,function(x){any(x)})
x.shape <- x[x.na.true,]
print(dim(x.shape))

return(list(shape = x.shape,all = x,n.vec = x.na.true,shape.method= paste("Columns raw:",dim(x)[2],"| Shape factor:",shape)))
	}
