####################################################################################################################################################
    # Phase I Xbar Chart Based on Champ and Jones(2004); see also Yao et al. (2017), and Yao and Chakraborti (2018)
####################################################################################################################################################

#use too many functions from package mvtnorm
require(mvtnorm)
#require(adehabitatLT)
#dchi <- adehabitatLT::dchi

#focus on these 3 functions from package mvtnorm


#pmvt <- mvtnorm::pmvt
#qmvt <- mvtnorm::qmvt
#pmvnorm <- mvtnorm::pmvnorm


####################################################################################################################################################

#source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/MKLswitch.R')

####################################################################################################################################################
    #parts of getting the charting constant l
####################################################################################################################################################

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function

corr.par.f <- function(nu, c4.option) {
        if (c4.option == TRUE) {
            corr.par <- c4.f(nu)
        } else {
            corr.par <- 1
        }
}

PH1.corr.f <- function(k, off.diag = - 1 / (k - 1)){
                                                                                   #correlation matrix
    crr <- diag(k)
    crr[which(crr == 0)] <- off.diag

    crr

}

#first.der.c4.f <- function(nu){
#
#    beta.part <- (beta(nu / 2, 1 / 2)) ^ (-1)
#
#    digamma1 <- digamma(nu / 2)
#    digamma2 <- digamma(nu / 2 + 1 / 2)
#
#    - sqrt(2 * pi) / 2 * beta.part * ((digamma1 - digamma2) / sqrt(nu) + nu ^ (- 3 /2))
#
#}
#
#first.der.cons.f <- function(nu, tau){
#
#    1 / 2 * (1 - 1 / (nu - tau)) ^ (- 1 / 2) * (nu - tau) ^ - 2
#
#}
#
#cons.f <- function(nu, tau){
#
#    sqrt((nu - tau - 1) / (nu - tau))
#
#}


####################################################################################################################################################
    #get l using the multivariate t cdf
####################################################################################################################################################

PH1.get.cc.mvt <- function(
                 k
                 ,nu
                 ,FAP = 0.1
                 #,Phase1 = TRUE
                 ,off.diag = -1/(k - 1)
                 #,alternative = '2-sided'
                 ,c4.option = TRUE
                 ,maxiter = 10000

){

    alternative = '2-sided'                                                   #turn off the alternative
                                                                              #The purpose of this function is
    #MCMC <- FALSE                                                             #to obtain l and k based on
                                                                              #multivariate t.
                                                                              #MCMC part is not available now.

    #if (off.diag == NULL) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))

    corr.P <- PH1.corr.f(k = k, off.diag = off.diag)                              #get correlation matrix with equal correlations

    pu <- 1 - FAP

    corr.par <- corr.par.f(nu, c4.option)

    #if (MCMC == TRUE) {
        #MVN.Q.Gibbs.Sampling(
        #    pu,
        #    MCMC.maxsim,
        #    corr.P,
        #    tails = alternative,
        #    burn = MCMC.burn,
        #    search.interval = MCMC.search.interval
        #)



    #} else {

        L <- ifelse(
                alternative == '2-sided',
                qmvt(pu, df = nu, sigma = corr.P, tail = 'both.tails', maxiter = maxiter)$quantile,
                qmvt(pu, df = nu, sigma = corr.P, maxiter = maxiter)$quantile
            )
                                                      #get L by multivariate T



    #}



    c.i <- L * corr.par * sqrt((k - 1) / k)             #get K

    list(l = L, c.i = c.i)

}

####################################################################################################################################################
    #get L by multivariate Normal
####################################################################################################################################################

#Old code is using chi distribution from package adehabitatLT
#
#PH1.joint.pdf.mvn.chisq <- function(Y, K, m, nu, sigma, alternative = '2-sided') {
#
#    s <- length(Y)
#
#    L <- K / sqrt((m - 1) / m * nu) * Y / c4.f(nu)
#
#    dpp <- lapply(
#            1:s,
#            function(i){
#
#                LL <- rep(L[i], m)
#
#                ifelse(
#                    alternative == '2-sided',
#                    pmvnorm(lower = -LL, upper = LL, sigma = sigma),
#                    pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
#                )
#
#            }
#
#    )
#
#    dpp <- unlist(dpp)
#
#    dpp * dchi(Y, df = nu)
#
#
#}

PH1.joint.pdf.mvn.chisq <- function(
                                Y
                                ,c.i
                                ,k
                                ,nu
                                ,sigma
                                #,alternative = '2-sided'
                                ,c4.option = TRUE)
{

    alternative = '2-sided'                                                   #turn off the alternative

    s <- length(Y)

    corr.par <- corr.par.f(nu, c4.option)

    L <- c.i / sqrt((k - 1) / k * nu) * sqrt(Y) / corr.par

    dpp <- lapply(
            1:s,
            function(i){

                LL <- rep(L[i], k)

                ifelse(
                    alternative == '2-sided',
                    pmvnorm(lower = -LL, upper = LL, sigma = sigma),
                    pmvnorm(lower = rep(-Inf, k), upper = LL, sigma = sigma)
                )

            }

    )

    dpp <- unlist(dpp)

    dpp * dchisq(Y, df = nu)


}

PH1.root.mvn.F <- function(
                    c.i
                    , k
                    , nu
                    , sigma
                    , pu
                    #, alternative = '2-sided'
                    , c4.option = TRUE
                    , subdivisions = 2000
                    , rel.tol = 1e-2)
{
    alternative = '2-sided'                                                   #turn off the alternative
                                                            #The purpose of this function is
    #s <- length(Y)                                          #to obtain appropriate L
    #                                                        #by multivariate normal
    #L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)          #
    #
    #pp <- lapply(
    #        1:s,
    #        function(i){
    #
    #            LL <- rep(L[i], m)
    #            ifelse(
    #                alternative == '2-sided',
    #                pmvnorm(lower = -LL, upper = LL, sigma = sigma),
    #                pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
    #            )
    #
    #        }
    #
    #)
    #
    #pp <- mean(unlist(pp))


    pp <- integrate(
            PH1.joint.pdf.mvn.chisq,
            lower = 0,
            upper = Inf,
            c.i = c.i,
            k = k,
            nu = nu,
            sigma = sigma,
            #alternative = alternative,
            c4.option = c4.option,
            subdivisions = subdivisions,
            rel.tol = rel.tol
        )$value

    pu - pp


}


PH1.get.cc.mvn <- function(
                 k
                 ,nu
                 ,FAP = 0.1
                 #,Phase1 = TRUE
                 ,off.diag = -1/(k - 1)
                 #,alternative = '2-sided'
                 ,c4.option = TRUE
                 ,interval = c(1, 7)
                 ,maxiter = 10000
                 ,subdivisions = 2000
                 ,tol = 1e-2

){
    alternative = '2-sided'                                                   #turn off the alternative
                                                         #The purpose of this function is
                                                            #to obtain L and K based on
    #MCMC <- FALSE                                           #multivariate normal.
                                                            #MCMC part is not available now.
    #if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))

    corr.par <- corr.par.f(nu, c4.option)

    corr.P <- PH1.corr.f(k = k, off.diag = off.diag)

    pu <- 1 - FAP

    #Y <- sqrt(rchisq(maxsim, df = nu))

    #if (MCMC == TRUE) {

        #MVN.Q.Gibbs.Sampling(
        #    pu,
        #    MCMC.maxsim,
        #    corr.P,
        #    tails = alternative,
        #    burn = MCMC.burn,
        #    search.interval = MCMC.search.interval
        #)



    #} else {


        c.i <- uniroot(
                PH1.root.mvn.F,
                interval = interval,
                k = k,
                nu = nu,
                #Y = Y,
                sigma = corr.P,
                pu = pu,
                #alternative = alternative,
                c4.option = c4.option,
                subdivisions = subdivisions,
                tol = tol,
                rel.tol = tol,
                maxiter = maxiter
        )$root

    #}

    L <- c.i / corr.par * sqrt(k / (k - 1))

    list(l = L, c.i = c.i)


}

#get.cc.mvn(20, 80)

####################################################################################################################################################

PH1.get.cc <- function(
            k
            ,nu
            ,FAP = 0.1
            ,off.diag = -1/(k - 1)
            #,alternative = '2-sided'
            ,c4.option = TRUE
            ,maxiter = 10000
            ,method = 'direct'
            ,indirect.interval = c(1, 7)
            ,indirect.subdivisions = 100L
            ,indirect.tol = .Machine$double.eps^0.25


){
    alternative = '2-sided'                                                   #turn off the alternative
                                                  #The purpose of this function is to obtain L and K
                                                    #by multivariate T or multivariate normal.
                                                    #Multivariate normal is so time-consuming
                                                    #that I do not recommend.

    #Phase1 <- TRUE                                  #need to delete if need to use Phase 2

	#method <- 'direct'								#need to delete if need to use indirect

    #if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))

    is.int <- ifelse(nu == round(nu), 1, 0)

	#if (is.int == 0) cat('Nu is not an integrer. Please check.', '\n') #need to delete if need to use indirect

    #if (method == 'direct' & is.int == 1) {                       #using multivariate T to obtain L and K

	if (is.int == 1 && method == 'direct') {
                                               #need to delete if need to use indirect
        PH1.get.cc.mvt(
            k = k
            ,nu = nu
            ,FAP = FAP
            #,Phase1 = Phase1
            ,off.diag = off.diag
            #,alternative = alternative
            ,c4.option = c4.option
            ,maxiter = maxiter
        )

    } else if (is.int == 1 && method == 'indirect'){

        cat('Nu is an integer. Using the indirect method may slow the computation process down.', '\n')

        PH1.get.cc.mvn(
            k = k
            ,nu = nu
            ,FAP = FAP
            #,Phase1 = Phase1
            ,off.diag = off.diag
            #,alternative = alternative
            ,c4.option = c4.option
            ,interval = indirect.interval
            #,maxsim = indirect.maxsim
            ,subdivisions = indirect.subdivisions
            ,maxiter = maxiter
            ,tol = indirect.tol
        )

    }

    else if (is.int == 0 && method == 'direct'){   #using multivariate normal to obtain L and K


        stop('Nu is not an integer. Please use the indirect method instead.')


    } else if (is.int == 0 && method == 'indirect') {

        PH1.get.cc.mvn(
            k = k
            ,nu = nu
            ,FAP = FAP
            #,Phase1 = Phase1
            ,off.diag = off.diag
            #,alternative = alternative
            ,c4.option = c4.option
            ,interval = indirect.interval
            #,maxsim = indirect.maxsim
            ,subdivisions = indirect.subdivisions
            ,maxiter = maxiter
            ,tol = indirect.tol
        )

    } else {

        stop('Unknown method. Please select the direct method or the indirect method.')

    }


}

####################################################################################################################################################
    #reverse the process
####################################################################################################################################################
#only support the direct method
#PH1.get.FAP0 <- function(
#            c.i
#            ,k
#            ,nu
#            ,off.diag = -1/(k - 1)
#            ,alternative = '2-sided'
#            ,c4.option = TRUE
#            ,maxiter = 10000
#            ,indirect.subdivisions = 100L
#            ,indirect.tol = .Machine$double.eps^0.25
#
#){
#
#    #method <- 'direct'
#
#    #Phase1 <- TRUE
#
#    is.int <- ifelse(nu == round(nu), 1, 0)
#
#    #if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
#
#    corr.P <- PH1.corr.f(k = k, off.diag = off.diag)
#
#    corr.par <- corr.par.f(nu, c4.option)
#
#    L <- c.i / corr.par * sqrt(k / (k - 1))
#
#
#    if (method == 'direct' & is.int == 1) {                       #using multivariate T to obtain L and K
#
#        ll <- rep(L, k)
#
#        ifelse(
#            alternative == '2-sided',
#            1 - pmvt(lower = -ll, upper = ll, df = nu, corr = corr.P, algorithm = TVPACK, abseps = 1e-12),
#            1 - pmvt(lower = -Inf, upper = ll, df = nu, corr = corr.P, algorithm = TVPACK, abseps = 1e-12)
#        )
#
#    }
#
#
#
#
#}

####################################################################################################################################################
    #Build Control Chart
####################################################################################################################################################

PH1XBAR <- function(
			X
            ,c.i = NULL
			,FAP = 0.1
			,off.diag = -1/(k - 1)
			#,alternative = '2-sided'
            ,model = 'ANOVA-based'
            ,c4.option = TRUE
			,plot.option = TRUE
			,maxiter = 10000
			,method = 'direct'
			,indirect.interval = c(1, 7)
			,indirect.subdivisions = 100L
			,indirect.tol = .Machine$double.eps^0.25
) {

    alternative = '2-sided'                                                   #turn off the alternative


	k <- dim(X)[1]
	n <- dim(X)[2]

	X.bar <- rowMeans(X)

	X.bar.bar <- mean(X)

    if (model == 'ANOVA-based') {

        nu <- k - 1

        corr.par <- corr.par.f(nu, c4.option)

        sigma.v <- sqrt(var(X.bar)) / corr.par

    } else if (model == 'basic') {

        nu <- k * (n - 1)

        corr.par <- corr.par.f(nu, c4.option)

        sigma.v <- sqrt(sum(apply(X, 1, var)) / k) / corr.par / sqrt(n)

    } else {

        stop("need to specify whether it is based on the ANOVA-based model or others")

    }

    if (is.null(c.i)) {
        c.i <- PH1.get.cc(
                k = k
                ,nu = nu
                ,FAP = FAP
                ,off.diag = off.diag
                #,alternative = alternative
                ,c4.option = c4.option
                ,maxiter = maxiter
                ,method = method
                ,indirect.interval = indirect.interval
                ,indirect.subdivisions = indirect.subdivisions
                ,indirect.tol = indirect.tol
         )$c.i
    } else {

            c.i <- c.i

    }

	LCL <- X.bar.bar - c.i * sigma.v
	UCL <- X.bar.bar + c.i * sigma.v

	if (plot.option == TRUE){

		plot(c(1, k), c(LCL, UCL), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Mean', type = 'n')

		axis(side = 1, at = 1:k)

		points(1:k, X.bar, type = 'o')
		points(c(-1, k + 2), c(LCL, LCL), type = 'l', col = 'red')
		points(c(-1, k + 2), c(UCL, UCL), type = 'l', col = 'red')
		points(c(-1, k + 2), c(X.bar.bar, X.bar.bar), type = 'l', col = 'blue')
		text(round(k * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
		text(round(k * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)


	}

	list(CL = X.bar.bar, sigma = sigma.v, PH1.cc = c.i, k = k, nu = nu, LCL = LCL, UCL = UCL, CS = X.bar)

}


####################################################################################################################################################
    #Example about table 6.3 (Montgomery, 2013)
####################################################################################################################################################
PH1XBAR.data <- function(){
	matrix(c(
		74.03	,	74.002	,	74.019	,	73.992	,	74.008	,
		73.995	,	73.992	,	74.001	,	74.011	,	74.004	,
		73.988	,	74.024	,	74.021	,	74.005	,	74.002	,
		74.002	,	73.996	,	73.993	,	74.015	,	74.009	,
		73.992	,	74.007	,	74.015	,	73.989	,	74.014	,
		74.009	,	73.994	,	73.997	,	73.985	,	73.993	,
		73.995	,	74.006	,	73.994	,	74	,	74.005	,
		73.985	,	74.003	,	73.993	,	74.015	,	73.988	,
		74.008	,	73.995	,	74.009	,	74.005	,	74.004	,
		73.998	,	74	,	73.99	,	74.007	,	73.995	,
		73.994	,	73.998	,	73.994	,	73.995	,	73.99	,
		74.004	,	74	,	74.007	,	74	,	73.996	,
		73.983	,	74.002	,	73.998	,	73.997	,	74.012	,
		74.006	,	73.967	,	73.994	,	74	,	73.984	,
		74.012	,	74.014	,	73.998	,	73.999	,	74.007	,
		74	,	73.984	,	74.005	,	73.998	,	73.996	,
		73.994	,	74.012	,	73.986	,	74.005	,	74.007	,
		74.006	,	74.01	,	74.018	,	74.003	,	74	,
		73.984	,	74.002	,	74.003	,	74.005	,	73.997	,
		74	,	74.01	,	74.013	,	74.02	,	74.003	,
		73.982	,	74.001	,	74.015	,	74.005	,	73.996	,
		74.004	,	73.999	,	73.99	,	74.006	,	74.009	,
		74.01	,	73.989	,	73.99	,	74.009	,	74.014	,
		74.015	,	74.008	,	73.993	,	74	,	74.01	,
		73.982	,	73.984	,	73.995	,	74.017	,	74.013

	), ncol = 5, byrow = T)
}





####################################################################################################################################################
    # Phase II Xbar chart. Please see Yao and Chakraborti (2018)
####################################################################################################################################################

PH2.inner.normal <- function(Z, Y, c.ii, k, nu = k - 1, c4.option = TRUE){

    corr.par <- corr.par.f(nu, c4.option)

    qn <- Z / sqrt(k) + c.ii / corr.par / sqrt(nu) * sqrt(Y)

    pnorm(qn)

}

PH2.CFAR.intgrand <- function(u, v, c.ii, k, nu = k - 1, c4.option = TRUE){

    Z <- qnorm(u)
    Y <- qchisq(v, nu)

    p1 <- PH2.inner.normal(Z, Y, c.ii, k, nu, c4.option)
    p2 <- PH2.inner.normal(Z, Y, -c.ii, k, nu, c4.option)

    1 - p1 + p2
}

PH2.CARL.intgrand <- function(u, v, c.ii, k, nu = k - 1, c4.option = TRUE) 1 / PH2.CFAR.intgrand(u, v, c.ii, k, nu, c4.option)


####################################################################################################################################################


PH2.get.cc.uc <- function(ARL0, k, nu = k - 1, c4.option = TRUE, interval = c(1, 10), u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

        PH2.root.finding.uc <- function(c.ii, k, nu = k - 1, c4.option = TRUE, ARL0, u, v) {

            ARL0 - mean(PH2.CARL.intgrand(u, v, c.ii, k, nu, c4.option))

    }

    k <- k

    rt <- uniroot(PH2.root.finding.uc, interval = interval, k = k, nu = nu, c4.option = c4.option, ARL0 = ARL0, u = u, v = v, tol = tol)$root

    rt


}

####################################################################################################################################################


PH2.root.finding.EPC <- function(p, k, nu = k - 1, c.ii, ARLb, c4.option = TRUE, u = runif(100000), v = runif(100000)){

    1 - p - mean(PH2.CARL.intgrand(u, v, c.ii, k, nu, c4.option) >= ARLb)

}

PH2.get.cc.EPC <- function(p, k, nu = k - 1, eps = 0.1, ARL0 = 370, c4.option = TRUE, interval = c(1, 10), u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

    ARLb <- (1 - eps) * ARL0

    rt <- uniroot(
            PH2.root.finding.EPC,
            interval = interval,
            p = p,
            k = k,
            nu = nu,
            ARLb = ARLb,
            c4.option = c4.option,
            u = u,
            v = v,
            tol = tol
        )$root

    rt

}

#rnum <- 100000
#u <- runif(rnum)
#v <- runif(rnum)

#debug(PH2.get.c.ii.EPC)
#PH2.get.c.ii.EPC(p = 0.05, k = 100, eps = 0.1, ARL0 = 500, c4.option = TRUE, interval = c(1, 5), u = u, v = v)


PH2.get.ARLb.EPC <- function(p, k, nu = k - 1, c.ii, ARL0 = NULL, c4.option = TRUE, interval = c(1, 1000), u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

    rt <- uniroot(
            PH2.root.finding.EPC,
            interval = interval,
            p = p,
            k = k,
            nu = nu,
            c.ii = c.ii,
            c4.option = c4.option,
            u = u,
            v = v,
            tol = tol
        )$root

    if (is.null(ARL0)) {
        list(ARLb = rt, eps = NULL, ARL0 = NULL)
    } else {
        list(ARLb = rt, eps = 1 - rt / ARL0, ARL0 = ARL0)
    }


}
#PH2.get.ARLb.EPC(p = 0.1, k = 100, c.ii = 3, c4.option = TRUE, u = u, v = v)
#PH2.get.ARLb.EPC(p = 0.1, k = 100, c.ii = 3, ARL0 = 370, c4.option = TRUE, u = u, v = v)

PH2.get.k.EPC <- function(p, c.ii, eps = 0.1, ARL0 = 370, c4.option = TRUE, interval = c(1000, 3000), model = 'ANOVA-based', n = 10, u = runif(100000), v = runif(100000), tol = .Machine$double.eps^0.25){

    if (model == 'ANOVA-based') {

        PH2.root.finding.k.EPC <- function(p, k, c.ii, ARLb, c4.option = TRUE, n = 10, u, v){

            nu <- k - 1

            PH2.root.finding.EPC(p, k, nu, c.ii, ARLb, c4.option, u, v)

        }


    } else if (model == 'basic') {

        PH2.root.finding.k.EPC <- function(p, k, c.ii, ARLb, c4.option = TRUE, n = 10, u, v){

            nu <- (n - 1) * k

            PH2.root.finding.EPC(p, k, nu, c.ii, ARLb, c4.option, u, v)

        }

    } else {

        stop("need to specify whether it is based on the ANOVA-based model or others")
    }


    ARLb <- (1 - eps) * ARL0

    rt <- uniroot(
            PH2.root.finding.k.EPC,
            interval = interval,
            p = p,
            c.ii = c.ii,
            ARLb = ARLb,
            c4.option = c4.option,
            n = n,
            u = u,
            v = v,
            tol = tol
        )$root

    list(k = rt, PH1.sample.size = rt * n)

}

####################################################################################################################################################


PH2XBAR <- function(
            X
            ,PH1.info = list(X = NULL, mu = NULL, sigma = NULL, k = NULL, n = NULL, model = "ANOVA-based")
            ,c.ii.info = list(c.ii = NULL, method = 'UC', ARL0 = 370, p = 0.05, eps = 0.1, interval.c.ii.UC = c(1, 3.2), interval.c.ii.EPC = c(1, 10), UC.tol = .Machine$double.eps^0.25, EPC.tol = .Machine$double.eps^0.25)
            #,alternative = '2-sided'
            ,c4.option = TRUE
            ,plot.option = TRUE
            ,maxsim = 100000
) {

    alternative = '2-sided'                                                   #turn off the alternative

    k.ii <- dim(X)[1]

    X.bar <- rowMeans(X)


    if (is.null(PH1.info$X)){

        k <- PH1.info$k

        if (PH1.info$model == 'ANOVA-based') {

            nu <- k - 1

        } else if (PH1.info$model == 'basic') {

            if (is.null(PH1.info$n)) {

                stop('Subgroup size needs to be defined')

            } else {

                n <- PH1.info$n
                nu <- k * (n - 1)

            }

        } else {

            stop('The model in Phase 1 needs to be defined')

        }

        corr.par <- corr.par.f(nu, c4.option)

        X.bar.bar <- PH1.info$mu
        sigma.v <- PH1.info$sigma / corr.par

    } else {

        PH1.chart <- PH1XBAR(X = PH1.info$X, c.i = 1, model = PH1.info$model, plot.option = FALSE)

        k <- PH1.chart$k
        nu <- PH1.chart$nu

        X.bar.bar <- PH1.chart$CL
        sigma.v <- PH1.chart$sigma

    }

    u <- runif(maxsim)
    v <- runif(maxsim)

    UC <- NULL
    EPC <- NULL

    n.c.ii <- 0

    UC.flag <- 0
    EPC.flag <- 0

    if (is.null(c.ii.info$c.ii)) {

        if (c.ii.info$method == 'both') {

            UC.flag <- 1
            EPC.flag <- 1

        } else if (c.ii.info$method == 'UC') {

            UC.flag <- 1

        } else if (c.ii.info$method == 'EPC'){

            EPC.flag <- 1

        } else {

            stop("Please check the method of obtaining charting constants in Phase 2")
        }

        if (UC.flag == 1) {

            n.c.ii <- n.c.ii + 1

            UC <- PH2.get.cc.uc(
                    ARL0 = c.ii.info$ARL0
                    , k = k
                    , nu = nu
                    , c4.option = c4.option
                    , interval = c.ii.info$interval.c.ii.UC
                    , u = u
                    , v = v
                    , tol = c.ii.info$UC.tol
                )

        }

        if (EPC.flag == 1){

            n.c.ii <- n.c.ii + 1

            EPC <- PH2.get.cc.EPC(
                    p = c.ii.info$p
                    , k = k
                    , nu = nu
                    , eps = c.ii.info$eps
                    , ARL0 = c.ii.info$ARL0
                    , c4.option = c4.option
                    , interval = c.ii.info$interval.c.ii.EPC
                    , u = u
                    , v = v
                    , tol = c.ii.info$EPC.tol
                )

        }

        c.ii <- list(UC = UC, EPC = EPC)

    } else {

            n.c.ii <- length(c.ii.info$c.ii)

            c.ii <- c.ii.info$c.ii

    }

    LCL <- X.bar.bar - unlist(c.ii) * sigma.v
    UCL <- X.bar.bar + unlist(c.ii) * sigma.v

    if (plot.option == TRUE){

        plot(c(1, k.ii), c(min(c(LCL, X.bar)), max(c(UCL, X.bar))), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Mean', type = 'n')

        points(1:k.ii, X.bar, type = 'o')

        axis(side = 1, at = 1:k.ii)
        points(c(-1, k.ii + 2), c(X.bar.bar, X.bar.bar), type = 'l', col = 'blue')

        for (ii in 1:n.c.ii){

            points(c(-1, k.ii + 2), c(LCL[ii], LCL[ii]), type = 'l', col = 'red')
            points(c(-1, k.ii + 2), c(UCL[ii], UCL[ii]), type = 'l', col = 'red')

            if (is.null(c.ii.info$c.ii)) {

                lab.in <- paste('(', names(c.ii)[ii], ')', sep = '')

            } else {

                lab.in <- NULL

            }

            ulab <- paste('UCL', lab.in, ' = ', round(UCL[ii], 4), sep = '')
            llab <- paste('LCL', lab.in, ' = ', round(LCL[ii], 4), sep = '')

            text(round(k.ii * 0.8), UCL[ii], ulab, pos = 1)
            text(round(k.ii * 0.8), LCL[ii], llab, pos = 3)

        }




    }

    list(CL = X.bar.bar, sigma = sigma.v, k = k, nu = nu, PH2.cc = unlist(c.ii), LCL = LCL, UCL = UCL, CS = X.bar)

}

####################################################################################################################################################



PH2XBAR.data <- function(){

    matrix(
        c(
            240, 243, 250, 253, 248,
            238, 242, 245, 251, 247,
            239, 242, 246, 250, 248,
            235, 237, 246, 249, 246,
            240, 241, 246, 247, 249,
            240, 243, 244, 248, 245,
            240, 243, 244, 249, 246,
            245, 250, 250, 247, 248,
            238, 240, 245, 248, 246,
            240, 242, 246, 249, 248,
            240, 243, 246, 250, 248,
            241, 245, 243, 247, 245,
            247, 245, 255, 250, 249,
            237, 239, 243, 247, 246,
            242, 244, 245, 248, 245,
            237, 239, 242, 247, 245,
            242, 244, 246, 251, 248,
            243, 245, 247, 252, 249,
            243, 245, 248, 251, 250,
            244, 246, 246, 250, 246,
            241, 239, 244, 250, 246,
            242, 245, 248, 251, 249,
            242, 245, 248, 243, 246,
            241, 244, 245, 249, 247,
            236, 239, 241, 246, 242,
            243, 246, 247, 252, 247,
            241, 243, 245, 248, 246,
            239, 240, 242, 243, 244,
            239, 240, 250, 252, 250,
            241, 243, 249, 255, 253
        )
        ,nrow = 30
        ,ncol = 5
        ,byrow = TRUE
    )

}

#PH2.get.k.EPC(p = 0.05, c.ii = 3, eps = 0.2, ARL0 = 370, c4.option = TRUE, interval = c(100, 400), type = 'ANOVA-based', subgroup.size = 10, u = u, v = v)


#
#m <- dim(x)[1]
#n <- dim(x)[2]
#
#x.bar.bar <- mean(x)
#
#sigma.v <- sqrt(sum(apply(x, 1, var)) / m) / c4.f(m * (n - 1))
#
#k <- get.cc(m, m * (n - 1), 0.05)$K
#
#LCL <- x.bar.bar - k * sigma.v / sqrt(n)
#UCL <- x.bar.bar + k * sigma.v / sqrt(n)
#
#
#x.bar <- rowMeans(x)
#
#plot(c(1, m), c(LCL * 49997/50000, UCL * 50003/50000), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Means', type = 'n')
#
#axis(side = 1, at = 1:25)
#
#points(1:m, x.bar, type = 'o')
#points(c(-1, m + 2), c(LCL, LCL), type = 'l', col = 'red')
#points(c(-1, m + 2), c(UCL, UCL), type = 'l', col = 'red')
#points(c(-1, m + 2), c(x.bar.bar, x.bar.bar), type = 'l', col = 'blue')
#text(20, UCL * 50001/50000, paste('UCL = ', round(UCL, 4)))
#text(20, LCL * 49999/50000, paste('LCL = ', round(LCL, 4)))
#
#pos <- c(rep(c(3, 1), 5), 1, 3, 3, 1, 3, 1, 1, 3, 1, 3, 1, 3, 1, 3, 1) #rep(c(1, 3), 8))[-26]
#
#text(1:m, x.bar, round(x.bar, 4), cex = 0.7, pos = pos)
#
#
####################################################################################################################################################
    #Example about how to get K and L by FAP0
####################################################################################################################################################
#get.cc(
#            20
#            ,80
#            ,FAP = 0.1
#            ,off.diag = NULL
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'direct'
#)
#get.cc(
#            20
#            ,80
#            ,FAP = 0.1
#            ,off.diag = NULL
#            ,alternative = '2-sided'
#            ,maxiter = 10000
#            ,method = 'indirect'
#)

####################################################################################################################################################
    #Example about how to get FAP by K
####################################################################################################################################################
#
#k <- c(
#22.84794
#,10.52163
#,7.386245
#,6.057107
#,5.340985
#,4.898076
#,4.598652
#,4.383316
#,3.842428
#,3.619871
#,3.498807
#,3.422752
#,3.280788
#,3.182391
#,3.151007
#,3.135567
#,3.126384
#,3.120294
#
#
#
#
#
#)
#
#m.seq <- c(
#    3
#    ,4
#    ,5
#    ,6
#    ,7
#    ,8
#    ,9
#    ,10
#    ,15
#    ,20
#    ,25
#    ,30
#    ,50
#    ,100
#    ,150
#    ,200
#    ,250
#    ,300
#
#)
#
#
#i <- 0
#
#record <- rep(NA, 18)
#
#for (m in m.seq){
#
#    i <- i + 1
#
#    nu <- m - 1
#
#    record[i] <- get.FAP0(k[i], m, nu)
#
#}
#
#as.matrix(record)




