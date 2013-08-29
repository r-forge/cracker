hz.oneside.ttest <- 
function(temp.i.aov){
	
t.test.one.side.agg <- aggregate(temp.i.aov[,1],list(temp.i.aov[,1]),length)

t.test.one.side.1		<- t.test.one.side.agg[t.test.one.side.agg[,2]==1,1]
t.test.one.side.2		<- t.test.one.side.agg[t.test.one.side.agg[,2]>1,1]

one.side.test.results <- c()
for(i.one.side in 1:length(t.test.one.side.1)){
		combinations 		<- rbind(t.test.one.side.1[i.one.side],t.test.one.side.2)
	for(i.one.side.2 in 1:dim(combinations)[2]){
		combinations.i  <- combinations[,i.one.side.2 ]
	 	temp.order 		<- c(grep(t.test.one.side.1[i.one.side],combinations.i),grep(t.test.one.side.1[i.one.side],combinations.i,invert = T))

		temp.i.one.side.test <- temp.i.aov[,2][temp.i.aov[,1] == combinations.i[temp.order[2]]]
		temp.i.one.side.mu	 <- temp.i.aov[,2][temp.i.aov[,1] == combinations.i[temp.order[1]]]

		if(mean(temp.i.one.side.test,na.rm = T) > temp.i.one.side.mu){
		alternative.test <- "greater"	
		}else{
		alternative.test <- "less"				
		}	
		if(!all(temp.i.one.side.test == temp.i.one.side.mu)){
error.try<- class(try(		one.side.test <- t.test(x=temp.i.one.side.test,mu = temp.i.one.side.mu,alternative = alternative.test)
)
)		}else{
			one.side.test <- list(p.value = 1)
			
		}
		if(error.try == "error-try"){
						one.side.test <- list(p.value = 1)

		}
		
		one.side.test.results <- rbind(one.side.test.results,c(one.side.test$p.value,combinations.i,temp.i.aov[,1][temp.i.aov[,1] == combinations.i[temp.order[1]]]))
	}
	


}

.matrix <- as.matrix(one.side.test.results[,1])
rownames(.matrix ) <- one.side.test.results[,2]
rownames(.matrix ) <- one.side.test.results[,3]

return(list(one.side.test.results=one.side.test.results,test = one.side.test))

}