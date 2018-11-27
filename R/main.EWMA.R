source('https://raw.githubusercontent.com/bolus123/EWMA-Chart-in-Nonstandard-Situations/master/R/head.R')

a <- EWMA.get.cc.MC(ARL0 = 370, interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.1, mm = 20, ss = Inf, tt = 50, reltol = 1e-6, tol = 1e-6)
    
b <- EWMA.get.cc.Conditional.MC(p0 = 0.1, interval = c(1, 6.5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = 0.1, eplison = 0.1, ARL0 = 370, mm = 20, ss = Inf, tt = 50, reltol = 1e-6, tol = 1e-6)   

######################################################################################################

require(parallel)

cl <- makeCluster(detectCores() - 1)

clusterEvalQ(cl, source('https://raw.githubusercontent.com/bolus123/EWMA-Chart-in-Nonstandard-Situations/master/R/head.R'))    
    
    
######################################################################################################
#head.addr <- '/home/yyao17/SPC/EWMA/ChartingConstants/Code/head.R'
######################################################################################################

start.time <- Sys.time()

ARL0.vec <- c(370, 500)

lambda.vec <- c(0.2, 0.3, 0.4)

m.vec <- c(15, 20, 25, 30, 50, 100, 150, 200)

in.vec <- expand.grid(lambda.vec, ARL0.vec, m.vec)

len <- dim(in.vec)[1]

clusterExport(cl, c('in.vec'))

cc.cuc.vec <- unlist(parLapplyLB(
	cl, 
	X = 1:len, 
	function(X) {EWMA.get.cc.MC(ARL0 = in.vec[X, 2], interval = c(1, 5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = in.vec[X, 1], mm = in.vec[X, 3], ss = Inf, tt = 100, reltol = 1e-6, tol = 1e-6)}
))

end.time <- Sys.time()

result.cuc <- list(time = end.time - start.time, result = cbind(in.vec, as.matrix(cc.cuc.vec, ncol = 1)))

save(result.cuc, file = '/home/yuhuiyao/Documents/EWMA/code/result.cuc.Rdata')

end.time - start.time

######################################################################################################

start.time <- Sys.time()

ARL0.vec <- c(370, 500)#c(370, 500)

p0.vec <- c(0.05, 0.1)

lambda.vec <- c(0.2, 0.3)

eplison.vec <- c(0.1, 0.2)

m.vec <- c(30, 50, 100)

in.vec <- expand.grid(p0.vec, lambda.vec, eplison.vec, ARL0.vec, m.vec)

len <- dim(in.vec)[1]

clusterExport(cl, c('in.vec'))

cc.con.vec <- unlist(parLapplyLB(
	cl, 
	X = 1:len, 
	function(X) {EWMA.get.cc.Conditional.MC(p0 = in.vec[X, 1], interval = c(1, 6.5), xmin = 0, xmax = 1, 
	ymin = 0, ymax = 1, lambda = in.vec[X, 2], eplison = in.vec[X, 3], ARL0 = in.vec[X, 4], mm = in.vec[X, 5], ss = Inf, tt = 100, reltol = 1e-6, tol = 1e-6)}
))
		
	
end.time <- Sys.time()

result.con <- list(time = end.time - start.time, result = cbind(in.vec, as.matrix(cc.con.vec, ncol = 1)))

save(result.con, file = '/home/yuhuiyao/Documents/EWMA/code/result.con.Rdata')

end.time - start.time

######################################################################################################

stopCluster(cl)
