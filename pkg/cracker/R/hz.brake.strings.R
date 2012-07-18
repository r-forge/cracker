hz.brake.strings <-
function(x,n = 20){
	x <- as.character(x)
	n <- as.numeric(n)
	if(is.na(n)){
		n = 20
	}
	
	print(n)
	print(nchar(x))
	if(nchar(x) >n){

		for(i in 1:ceiling(nchar(x)/n)){
			
			x.i	<- substring(x, (n*i-n)+1,n*i)
			if(i == 1){x.new <- x.i}else{
			x.new <- paste(x.new,x.i,sep = "\n")
			}
		}
		x <- x.new
	}
	return(x)
}
