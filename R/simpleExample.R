source('https://raw.githubusercontent.com/bolus123/EWMA-Chart-in-Nonstandard-Situations/master/R/head.R')

a <- EWMA.get.cc.MC(ARL0 = 370, interval = c(1, 5), xmin = 0, xmax = 1, 
                    ymin = 0, ymax = 1, lambda = 0.1, mm = 20, ss = Inf, tt = 50, reltol = 1e-6, tol = 1e-6)

b <- EWMA.get.cc.Conditional.MC(p0 = 0.1, interval = c(1, 6.5), xmin = 0, xmax = 1, 
                                ymin = 0, ymax = 1, lambda = 0.1, eplison = 0.1, ARL0 = 370, mm = 20, ss = Inf, tt = 50, reltol = 1e-6, tol = 1e-6)  