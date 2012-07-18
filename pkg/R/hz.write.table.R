hz.write.table <-
function(x, temp.i.aov = temp.i.aov,name =name,acc = .acc,exp.des = exp.des){


	

#library(gridExtra)
	grid.newpage()
 	uni.acc	<- sort(unique(acc))

	######
	
	test  <- x
	test  <- test$p.value
	names <- uni.acc[i]
	names.vec <- exp.des[,1][-1]
	#print(names.vec)

	# names.vec <- names.vec[exclude.vec<- as.character(unique(c(colnames(test),rownames(test))))]

	
	#print(names.vec)

	#colnames(test) <- names.vec[-dim(test)[2]]	
	#rownames(test) <- exclude.vec[-1]	
	

	
	#	if(length(print(grep("at1g12080",name))) > 0){stop()}

		if(.data2$gui.input$raw){
			run.type <- 1		
		}else{
			run.type <- 2
		}
		exp.des <- cbind(.data2$exp.design[, run.type],1:length(.data2$exp.design[, run.type]))
		exp.des <- rbind(c("name","alias"),exp.des)

	
#print(i)
#print(dim(test))
	if(length(test) == 0){
	test <- matrix(NA,ncol = 2,nrow = 2)
	test[2,2] <- "no comparissons"
	
	}else{
	test <- round(test,digit = 5)
		}
	test <- rbind(colnames(test),test)
	test <- cbind(c(rownames(test)),test)
	#print(i)

	p.col <- test
	p.col[is.na(p.col)] <- "grey"
	p.col[p.col > 0.05] <- "grey"
	p.col[p.col <= 0.05] <- 2
	
	
	
	x1 <- 1/(dim(test)[1]+2) 
	y1 <- 1/(dim(test)[2]+2)
	
	x1 <- seq(from = x1, to = x1*dim(test)[1] , by = x1)
	x1 <- rep(x1,(dim(test)[2]))
	
	y1 <- seq(from = y1, to = y1*dim(test)[1] , by = y1)
	y1 <- rep(y1,each = dim(test)[2])

	test[1,1]	<- name
	p.col[1,]	<- 1
	p.col[,1]	<- 1
	p.col[1,1] 	<- 4
	
	rot <- matrix(1,dim(test),dim(test))
	rot[,1] <- 90
	just <- matrix("left",dim(test),dim(test))


	
	rownames(test) 	<- 1:dim(test)[1]
	colnames(test)	<- 1:dim(test)[2]
	
		grid.text(t(test),x = (x1),y = rev(y1),draw = TRUE,just = c(0,1),vjust = "left",gp=gpar(col = t(p.col),fontsize = 9),rot = rot, check.overlap = TRUE)
	
	return(x)
}
