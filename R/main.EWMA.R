require(parallel)

######################################################################################################
head.addr <- 'C:/Users/bolus/OneDrive/文档/GitHub/EWMA-Chart-in-Nonstandard-Situations/R/head.R'
######################################################################################################


ARL0.vec <- c(370, 500)

lambda.vec <- c(0.2, 0.3, 0.4)

m.vec <- c(100, 200)

in.vec <- expand.grid(ARL0.vec, lambda.vec, m.vec)

len <- dim(in.vec)[1]

cl <- makeCluster(detectCores() - 1)

clusterExport(cl, c('in.vec'))
clusterEvalQ(cl, source('C:/Users/bolus/OneDrive/文档/GitHub/EWMA-Chart-in-Nonstandard-Situations/R/head.R'))

cc.vec <- unlist(parLapply(
	cl, 
	X = 1:len, 
	function(X) EWMA.get.cc.MC(ARL0 = in.vec[X, 1], interval = c(2.3, 3.5), xmin = 0, xmax = 1, 
		ymin = 0, ymax = 1, lambda = in.vec[X, 2], mm = in.vec[X, 3], ss = Inf, tt = 10, reltol = 1e-6)
))
		
stopCluster(cl)

result <- cbind(in.vec, cc.vec)

save(result, file = 'C:/Users/bolus/OneDrive/文档/GitHub/EWMA-Chart-in-Nonstandard-Situations/R/Rcc.result.Rdata')