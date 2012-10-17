hz.change.nodePar <-
function(x,sclus,temp.col,temp.lwd,col.temp){
	if(is.leaf(x)){
		att.n <- attributes(x)
		it<- grep(att.n$label,sclus$labels)
		it	<<- it+1
				if(!any(is.na(temp.col))){

		attr(x,"nodePar") <- c(	att.n$nodePar,list(	lab.col = temp.col[it],
											pch = "",
											lab.cex = 1.2
											
											 )
											)
											}
		if(any(!is.na(col.temp))){
		attr(x,"edgePar") <- c(	att.n$edgePar,list(	lwd = temp.lwd	,
											col = col.temp[it]
																			 )
											)
		}



		}
		return(x)
	
	}
	
	
