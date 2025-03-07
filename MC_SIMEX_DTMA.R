#  Copyright 2025 Annamaria Guolo (University of Padova) 
#  Permission to use, copy, modify and distribute this software and
#  its documentation, for any purpose and without fee, is hereby granted,
#  provided that:
#  1) this copyright notice appears in all copies and
#  2) the source is acknowledged through a citation to the paper
#     Guolo A. and Caruso E. (2025). Misclassification SIMEX in meta-analysis of accuracy of diagnostic tests. Submitted.
#  The Authors make no representation about the suitability of this software
#  for any purpose.  It is provided "as is", without express or implied warranty

## parameter vector theta=c(mu.eta, mu.xi, sigma2.eta, sigma2.xi, rho)
## vector of information xx = eta.obs, xi.obs, var.etaobs, var.xiobs)

library(mvtnorm)
library(far)
library(simex)

mcsimex.results <- function(TP, FN, TN, FP, P=NULL, N=NULL,
                            lambda=c(0.0, 0.5, 1.0, 1.5, 2.0),
                            B=100, p=5,
                            extr.method='quadratic2'){
    ans <- list()
    P <- TP + FN
    N <- TN + FP
    data.bin <- data.frame(TP, P, FP, N, TN, FN)
    n <- NROW(data.bin)
    M <- length(lambda)
    sp.obs <- se.obs <- rep(NA, n)
    for(i in 1:n){
        sp.obs[i] <- data.bin$TN[i]/data.bin$N[i]  
        se.obs[i] <- data.bin$TP[i]/data.bin$P[i]
        if(sp.obs[i]==1 | sp.obs[i]==0)
            sp.obs[i] <-  (data.bin$TN[i]+0.5)/(data.bin$N[i]+1)
        if(se.obs[i]==1 | se.obs[i]==0)
            se.obs[i] <- (data.bin$TP[i]+0.5)/(data.bin$P[i]+1)
    }
    
    prepare.data <- function(data.all.lambda){
        data.ready <- list()  
        for(b in 1:B){
            data.b <- data.all.lambda[b,,] 
            data.b <- data.frame(data.b,
                                 'P'=data.b[,'TP'] + data.b[,'FN'],
                                 'N'=data.b[,'FP'] + data.b[,'TN'])
            this.sp.obs <- this.se.obs <- rep(NA, n)
            for(i in 1:n){
                this.sp.obs[i] <- data.b$TN[i]/data.b$N[i]
                this.se.obs[i] <- data.b$TP[i]/data.b$P[i]
                if(this.sp.obs[i]==1 | this.sp.obs[i]==0){
                    this.sp.obs[i] <- (data.b$TN[i] + 0.5)/(data.b$N[i] + 1)
                }
                if(this.se.obs[i]==1 | this.se.obs[i]==0){
                    this.se.obs[i] <- (data.b$TP[i] + 0.5)/(data.b$P[i] + 1)
                }
            }
            data.ready[[b]] <- cbind(this.se.obs, this.sp.obs)
        }
        return(data.ready)
    }
    
    est.to.function <- function(theta, x){  ##single b
        eta.obs <- x[,1]
        xi.obs <- x[,2]
        if(theta[3] < 0 | theta[4] < 0 | theta[5] > 1 | theta[5] < -1)
            return(NA)
        else{
            val <-  sum(log(dnorm(xi.obs, theta[2], sqrt(theta[4])) *
                            dnorm(eta.obs, theta[1]+theta[5]*sqrt(theta[3]/theta[4])*(xi.obs-theta[2]), sqrt(theta[3]*(1-theta[5]*theta[5])))))
            return(val)
        }
    }
    
    est.function <- function(x){ ## x: columns of data simulation in SIM.step
        x[,1] <- qlogis(x[,1])
        x[,2] <- qlogis(x[,2])
        theta <- c(mean(x[,1]),
                   mean(x[,2]),
                   var(x[,1]),
                   var(x[,2]),
                   cor(x[,1], x[,2]))
        val <- try(optim(theta, est.to.function, x=x, control=list(fnscale=-1), hessian=TRUE) , silent=TRUE)
        if(class(val)=='try-error')
            val <- try(optim(theta, est.to.function, x=x, control=list(fnscale=-1), hessian=TRUE, method='BFGS') , silent=TRUE)
        if(class(val)=='try-error'){
            xx <- cbind(se.obs, sp.obs)
            xx[,1] <- qlogis(xx[,1])
            xx[,2] <- qlogis(xx[,2])
            theta <- c(mean(xx[,1]),
                       mean(xx[,2]),
                       var(xx[,1]),
                       var(xx[,2]),
                       cor(xx[,1], xx[,2]))
            val <- try(optim(theta, est.to.function, x=xx, control=list(fnscale=-1), hessian=TRUE) , silent=TRUE)
            if(class(val)=='try-error')
                val <- try(optim(theta, est.to.function, x=xx, control=list(fnscale=-1), hessian=TRUE, method='BFGS') , silent=TRUE)
        }
        if(class(val)!='try-error'){
            theta.original <- val$par
            variance.theta <- solve(-val$hessian)
            hessian.theta <- val$hessian        
        }
        else {
            theta.original <- rep(NA, length(theta))
            variance.theta <- hessian.theta <- diag(NA, length(theta))
        }
        names(theta.original) <- colnames(variance.theta) <- c('eta', 'xi', 'vareta', 'varxi', 'rho')
        return(list(theta.original=theta.original,
                    variance.theta=variance.theta,
                    hessian.theta = hessian.theta))
    }
    
    SIM.step <- function(lambda.notzero){ ## lambda != 0
        values <- list()
        data.tmp <- array(NA, dim=c(B, n, 4), dimnames=list(1:B, 1:n, c('TP','FN','FP','TN')))
        vec.gradient <- matrix(0.0, nrow=n, ncol=p)
        hess.matrix <- matrix(0.0, nrow=n, ncol=p^2)
        for(m in 1:length(lambda.notzero)){
            for(i in 1:n){
                ## simulation
                mc.simex <- matrix(c(sp.obs[i], 1-se.obs[i], 1-sp.obs[i], se.obs[i]), ncol=2, byrow=TRUE)
                xx <- factor(c(rep(0, data.bin[i,'TN'] + data.bin[i,'FN']),
                               rep(1, data.bin[i,'FP'] + data.bin[i,'TP'])))
                x <- data.frame(x1=xx)
                colnames(mc.simex) <- levels(xx)
                mc.matrix <- list(x1 = mc.simex)
                x.mc <- try(replicate(B,
                                      misclass(data.org = x, mc.matrix = mc.matrix, k = lambda.notzero[m]),
                                      simplify='array'),
                            silent=TRUE)
                if(class(x.mc)=='try-error')
                    print(paste('SE + SP -1 smaller than 0 in study ', i, '; lambda equal to ', lambda.notzero[m]))
                x.true <- c(rep(0, data.bin[i,'N']), rep(1, data.bin[i,'P']))
                tab.mc <- lapply(x.mc, function(w) table(w, x.true))
                data.tmp[,i,'TP'] <- unlist(lapply(tab.mc, function(w) w[2,2]))
                data.tmp[,i,'FN'] <- unlist(lapply(tab.mc, function(w) w[1,2]))
                data.tmp[,i,'TN'] <- unlist(lapply(tab.mc, function(w) w[1,1]))               
                data.tmp[,i,'FP'] <- unlist(lapply(tab.mc, function(w) w[2,1]))
            }
            ## estimation for fixed lambda
            data.ready.lambda <- prepare.data(data.tmp)
            all.theta.b.lambda <- lapply(data.ready.lambda, est.function)
            theta.b <- do.call(rbind,lapply(all.theta.b.lambda, function(x) rbind(x[[1]])))
            all.data <- mapply(list, data.b=data.ready.lambda, theta.b=all.theta.b.lambda, SIMPLIFY = FALSE)
            var.theta.b <- do.call(rbind, lapply(all.theta.b.lambda, function(x) c(x[[2]])))
            index <- which(!is.na(apply(var.theta.b, 1, sum)))
            theta.b <- theta.b[index,]
            var.theta.b <- var.theta.b[index,]
            theta.lambda <- apply(theta.b, 2, mean)
            var.theta.b <- matrix(apply(var.theta.b, 2, mean), ncol=p)
            var.theta.b2 <- var(theta.b)           
            values[[m]] <- list(theta.lambda=theta.lambda,
                                var.lambda=var.theta.b - var.theta.b2)
        }
        return(values)  
    }
    
    
    ## EXtrapolation step
    ## theta.hat : SE, SP
    theta.sx <- var.sx <- matrix(0.0, nrow=length(lambda), ncol=p)
    colnames(theta.sx) <- colnames(var.sx) <- c('eta', 'xi', 'vareta', 'varxi', 'rho')
    ris0 <- est.function( cbind(se.obs, sp.obs))
    theta.sx.lambda0 <- ris0$theta
    var.sx.complete.lambda0 <- ris0$variance.theta
    var.sx.lambda0 <- diag(ris0$variance.theta)
    ris <- SIM.step(lambda[-1])    
    theta.sx <- do.call(rbind, lapply(ris, function(x) rbind(x[[1]])))
    var.sx <- do.call(rbind, lapply(ris, function(x) rbind(diag(x[[2]]))))
    var.sx.complete.tmp <- list2array(lapply(ris, function(x)(x[[2]])))
    theta.sx <- rbind(theta.sx.lambda0, theta.sx)
    var.sx <- rbind(var.sx.lambda0, var.sx)
    var.sx.complete <- array(NA, dim=c(p, p, length(lambda))) 
    var.sx.complete[, , 1] <- var.sx.complete.lambda0
    var.sx.complete[, , -1] <- var.sx.complete.tmp
    
    ## extrapolation
    extr.theta <- matrix(0.0, ncol=p, nrow=1)
    extr.var <- matrix(0.0, ncol=p, nrow=1)
    extr.var.complete <- matrix(0.0, ncol=p, nrow=p)
    ## columns: p parameters, rows: m for extrapolation
    
    if(extr.method=='quadratic'){
        extrapolation.theta <- lm(theta.sx ~ lambda + I(lambda^2))
        extr.theta <- predict(extrapolation.theta, newdata = data.frame(lambda = -1))
        extrapolation.var <- lm(var.sx ~ lambda + I(lambda^2))
        extr.var <- predict(extrapolation.var, newdata = data.frame(lambda = -1))
    }
    if(extr.method=='quadratic2'){
        extrapolation.theta <- lm(theta.sx ~ I(lambda^2))
        extr.theta <- predict(extrapolation.theta, newdata = data.frame(lambda = -1))
        extrapolation.var <- lm(var.sx ~ I(lambda^2))
        extr.var <- predict(extrapolation.var, newdata = data.frame(lambda = -1))
    }
    diag(extr.var.complete) <- extr.var
    if(extr.method=='quadratic'){
        for(j in 1:(p-1)){
            values <- var.sx.complete[,j,j+1]
            extrapolation.var.complete <- lm(values ~ lambda + I(lambda^2))
            extr.var.complete[j,j+1] <- extr.var.complete[j+1,j] <- predict(extrapolation.var.complete, newdata = data.frame(lambda = -1))
        }
    }
    if(extr.method=='quadratic2'){
        for(j in 1:(p-1)){
            values <- var.sx.complete[,j,j+1]
            extrapolation.var.complete <- lm(values ~ I(lambda^2))
            extr.var.complete[j,j+1] <- extr.var.complete[j+1,j] <- predict(extrapolation.var.complete, newdata = data.frame(lambda = -1))
        }
    }
    
    colnames(extr.theta) <- c('mu.eta', 'mu.xi', 'var.eta', 'var.xi', 'rho')
    ans$theta <- extr.theta
    colnames(extr.var) <- c('mu.eta', 'mu.xi', 'var.eta', 'var.xi', 'rho')
    ans$se <- sqrt(extr.var)
    rownames(extr.var.complete) <- colnames(extr.var.complete) <- c('mu.eta', 'mu.xi', 'var.eta', 'var.xi', 'rho')
    ans$var.matrix <- extr.var.complete
    rownames(theta.sx) <- lambda
    ans$all.theta <- theta.sx
    ans$lambda <- lambda
    ans$extr.method <- extr.method
    class(ans) <- "mcsimex.results"
    return(ans)
}

