"coef.lars" <-
function(object, ...)
{
	predict(object, type = "coefficient", ...)$coef
}
"cv.folds" <-
function(n, folds = 10)
{
	split(sample(1:n), rep(1:folds, length = n))
}
"cv.lars" <-
function(x, y, K = 10, fraction = seq(from = 0, to = 1, length = 100), 
           trace = FALSE, plot.it = TRUE, se = TRUE, ...)
{
  all.folds <- cv.folds(length(y), K)
  residmat <- matrix(0, length(fraction), K)
  for(i in seq(K)) {
    omit <- all.folds[[i]]
    fit <- lars(x[ - omit,  ], y[ - omit], trace = trace, ...)
    fit <- predict(fit, x[omit,  ,drop=FALSE], mode = "fraction", s = fraction
                   )$fit
    if(length(omit)==1)fit<-matrix(fit,nrow=1)
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    if(trace)
      cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(fraction = fraction, cv = cv, cv.error = cv.error)
  if(plot.it) plotCVLars(object,se=se)
  invisible(object)
}
"downdateR" <-
function(R, k = p)
{
	p <- dim(R)[1]
	if(p == 1)
		return(NULL)
	R <- delcol(R, rep(1, p), k)[[1]][ - p,  , drop = FALSE]
	attr(R, "rank") <- p - 1
	R	# Built-in Splus utility
}
"error.bars" <-
function(x, upper, lower, width = 0.02, ...)
{
	xlim <- range(x)
	barw <- diff(xlim) * width
	segments(x, upper, x, lower, ...)
	segments(x - barw, upper, x + barw, upper, ...)
	segments(x - barw, lower, x + barw, lower, ...)
	range(upper, lower)
}
"lars" <-
  function(x, y, type = c("lasso", "lar", "forward.stagewise"), trace = FALSE, Gram, 
           eps = .Machine$double.eps,  max.steps, use.Gram = TRUE)
{
### program automatically centers and standardizes predictors.
###
### Original program by Brad Efron September 2001
### Recoded by Trevor Hastie November 2001
### Computational efficiency December 22, 2001
### Bug fixes and singularities February 2003
### Conversion to R April 2003
### Copyright Brad Efron and Trevor Hastie
  call <- match.call()
  type <- match.arg(type)
  TYPE <- switch(type,
                 lasso = "LASSO",
                 lar = "LAR",
                 forward.stagewise = "Forward Stagewise")
  if(trace)
    cat(paste(TYPE, "sequence\n"))
  nm <- dim(x)
  n <- nm[1]
  m <- nm[2]
  im <- inactive <- seq(m)
  one <- rep(1, n)
  vn <- dimnames(x)[[2]]	
### Center x and y, and scale x, and save the means and sds
  meanx <- drop(one %*% x)/n
  x <- scale(x, meanx, FALSE)	# centers x
  normx <- sqrt(drop(one %*% (x^2)))
  nosignal<-normx/sqrt(n) < eps
  if(any(nosignal))# ignore variables with too small a variance
    {
    ignores<-im[nosignal]
    inactive<-im[-ignores]
    normx[nosignal]<-eps*sqrt(n)
    if(trace)
      cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < \eps; dropped for good\n")	#
  }
  else ignores <- NULL #singularities; augmented later as well
  names(normx) <- NULL
  x <- scale(x, FALSE, normx)	# scales x
  if(use.Gram & missing(Gram)) {
    if(m > 500 && n < m)
      cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n"
          )
    if(trace)
      cat("Computing X'X .....\n")
    Gram <- t(x) %*% x	#Time saving
  }
  mu <- mean(y)
  y <- drop(y - mu)
  Cvec <- drop(t(y) %*% x)
  ssy <- sum(y^2)	### Some initializations
  residuals <- y
  if(missing(max.steps))
    max.steps <- 8*min(m, n-1)
  beta <- matrix(0, max.steps + 1, m)	# beta starts at 0
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL	# maintains active set
  actions <- as.list(seq(max.steps))	
                                        # a signed index list to show what comes in and out
  drops <- FALSE	# to do with type=="lasso" or "forward.stagewise"
  Sign <- NULL	# Keeps the sign of the terms in the model
  R <- NULL	###
### Now the main loop over moves
###
  k <- 0
  while((k < max.steps) & (length(active) < min(m - length(ignores),n-1)) )
    {
      action <- NULL
      k <- k + 1
      C <- Cvec[inactive]	#
### identify the largest nonactive gradient
      Cmax <- max(abs(C))	### Check if we are in a DROP situation
      if(!any(drops)) {
        new <- abs(C) >= Cmax - eps
        C <- C[!new]	# for later
        new <- inactive[new]	# Get index numbers
### We keep the choleski R  of X[,active] (in the order they enter)
        for(inew in new) {
          if(use.Gram) {
            R <- updateR(Gram[inew, inew], R, drop(Gram[
                                                        inew, active]), Gram = TRUE,eps=eps)
          }
          else {
            R <- updateR(x[, inew], R, x[, active], Gram
                         = FALSE,eps=eps)
          }
          if(attr(R, "rank") == length(active)) {
            ##singularity; back out
            nR <- seq(length(active))
            R <- R[nR, nR, drop = FALSE]
            attr(R, "rank") <- length(active)
            ignores <- c(ignores, inew)
            action <- c(action,  - inew)
            if(trace)
              cat("LARS Step", k, ":\t Variable", inew, 
                  "\tcollinear; dropped for good\n")	#
          }
          else {
            if(first.in[inew] == 0)
              first.in[inew] <- k
            active <- c(active, inew)
            Sign <- c(Sign, sign(Cvec[inew]))
            action <- c(action, inew)
            if(trace)
              cat("LARS Step", k, ":\t Variable", inew, 
                  "\tadded\n")	#
          }
        }
      }
      else action <-  - dropid
      Gi1 <- backsolve(R, backsolvet(R, Sign))	
### Now we have to do the forward.stagewise dance
### This is equivalent to NNLS
      dropouts<-NULL
      if(type == "forward.stagewise") {
        directions <- Gi1 * Sign
        if(!all(directions > 0)) {
          if(use.Gram) {
            nnls.object <- nnls.lars(active, Sign, R, 
                                     directions, Gram[active, active], trace = 
                                     trace, use.Gram = TRUE,eps=eps)
          }
          else {
            nnls.object <- nnls.lars(active, Sign, R, 
                                     directions, x[, active], trace = trace, 
                                     use.Gram = FALSE,eps=eps)
          }
          positive <- nnls.object$positive
          dropouts <-active[-positive]
          action <- c(action, -dropouts)
          active <- nnls.object$active
          Sign <- Sign[positive]
          Gi1 <- nnls.object$beta[positive] * Sign
          R <- nnls.object$R
          C <- Cvec[ - c(active, ignores)]
        }
      }
      A <- 1/sqrt(sum(Gi1 * Sign))
      w <- A * Gi1	# note that w has the right signs
      if(!use.Gram) u <- drop(x[, active, drop = FALSE] %*% w)	###
### Now we see how far we go along this direction before the
### next competitor arrives. There are several cases
###
### If the active set is all of x, go all the way
      if(length(active) >=  min(n-1, m - length(ignores) ) ) {
        gamhat <- Cmax/A
      }
      else {
        if(use.Gram) {
          a <- drop(w %*% Gram[active,  - c(active,ignores), drop = FALSE])
        }
        else {
          a <- drop(u %*% x[,  - c(active, ignores), drop=FALSE])
        }
        gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))	
### Any dropouts will have gam=0, which are ignored here
        gamhat <- min(gam[gam > eps], Cmax/A)	
      }
      if(type == "lasso") {
        dropid <- NULL
        b1 <- beta[k, active]	# beta starts at 0
        z1 <-  - b1/w
        zmin <- min(z1[z1 > eps], gamhat)
        if(zmin < gamhat) {
          gamhat <- zmin
          drops <- z1 == zmin
        }
        else drops <- FALSE
      }
      beta[k + 1,  ] <- beta[k,  ]
      beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
      if(use.Gram) {
        Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
      }
      else {
        residuals <- residuals - gamhat * u
        Cvec <- drop(t(residuals) %*% x)
      }
      Gamrat <- c(Gamrat, gamhat/(Cmax/A))
      arc.length <- c(arc.length, gamhat)	
### Check if we have to drop any guys
      if(type == "lasso" && any(drops)) {
        dropid <- seq(drops)[drops]	
                                        #turns the TRUE, FALSE vector into numbers
        for(id in rev(dropid)) {
          if(trace)
            cat("Lasso Step", k+1, ":\t Variable", active[
                                                        id], "\tdropped\n")
          R <- downdateR(R, id)
        }
        dropid <- active[drops]	# indices from 1:m
        beta[k+1,dropid]<-0  # added to make sure dropped coef is zero
        active <- active[!drops]
        Sign <- Sign[!drops]
      }
      if(!is.null(vn))
        names(action) <- vn[abs(action)]
      actions[[k]] <- action
      inactive <- im[ - c(active, ignores)]
    }
  beta <- beta[seq(k + 1),  ]	#
  dimnames(beta) <- list(paste(0:k), vn)	### Now compute RSS and R2
  if(trace)
    cat("Computing residuals, RSS etc .....\n")
  residuals <- y - x %*% t(beta)
  beta <- scale(beta, FALSE, normx)
  RSS <- apply(residuals^2, 2, sum)
  R2 <- 1 - RSS/RSS[1]
  Cp <- ((n - k - 1) * RSS)/rev(RSS)[1] - n + 2 * seq(k + 1)
  object <- list(call = call, type = TYPE, R2 = R2, RSS = RSS, Cp = Cp, 
                 actions = actions[seq(k)], entry = first.in, Gamrat = Gamrat, 
                 arc.length = arc.length, Gram = if(use.Gram) Gram else NULL, 
                 beta = beta, mu = mu, normx = normx, meanx = meanx)
  class(object) <- "lars"
  object
}
"nnls.lars" <-
  function(active, Sign, R, beta, Gram, eps = 1e-10, trace = FALSE, use.Gram = TRUE)
{
### Modified 05/15/03 to allow for more than one addition to the set
### Lawson and Hanson page 161
### Go back to the first positive coefficent vector; can assume its in order
### Note that X'y is constant for all these guys;
### we assume WOLOG this constant is 1
### We also assume we have come into this because we have a negative coeff  
### If use.Gram is FALSE, then Gram comes in as x
  if(!use.Gram) x <- Gram	# to avoid confusion
  M<-m <- length(active)
  im <- seq(m)
  positive <- im
  zero <- NULL
  ### Get to the stage where beta.old is all positive
  while(m>1) {
    zero.old<-c(m,zero)
    R.old <- downdateR(R, m)
    beta0 <- backsolve(R.old, backsolvet(R.old, Sign[ - zero.old]))*Sign[-zero.old]
    beta.old <- c(beta0,rep(0,length(zero.old)))
    if(all(beta0 >0))break
    m <-m-1
    zero<-zero.old
    positive<-im[-zero]
    R<-R.old
    beta<-beta.old
  }
### Now we do the NNLS backtrack dance
  while(TRUE) {
    while(!all(beta[positive] > 0)) {
      alpha0 <- beta.old/(beta.old - beta)
      alpha <- min(alpha0[positive][(beta <= 0)[positive]])
      beta.old <- beta.old + alpha * (beta - beta.old)
      dropouts<-match(alpha,alpha0[positive],0)
### Used to have the following line, but this failed occasionally
###   dropouts <- seq(positive)[abs(beta.old[positive]) < eps]
      for(i in rev(dropouts)) R <- downdateR(R, i)
      positive <- positive[ - dropouts]	
                                        # there is an order in R
      zero <- im[ - positive]
      beta0 <- backsolve(R, backsolvet(R, Sign[positive])) * 
        Sign[positive]
      beta <- beta.old * 0
      beta[positive] <- beta0
    }
### Now all those in have a positive coefficient
    if(use.Gram) {
      w <- 1 - Sign * drop(Gram %*% (Sign * beta))	
                                        #should be zero for some
    }
    else {
      jw <- x %*% (Sign * beta)
      w <- 1 - Sign * drop(t(jw) %*% x)
    }
    if((length(zero) == 0) || all(w[zero] <= 0))
      break
    add <- order(w)[M]
    if(use.Gram) {
      R <- updateR(Gram[add, add], R, drop(Gram[add, 
                                                positive]), Gram = TRUE,eps=eps)
    }
    else {
      R <- updateR(x[, add], R, x[, positive], Gram = FALSE,eps=eps)
    }
    positive <- c(positive, add)
    zero <- setdiff(zero, add)
    beta0 <- backsolve(R, backsolvet(R, Sign[positive])) * Sign[
                                                        positive]
    beta[positive] <- beta0
  }
  if(trace)
    {
      dropouts<-active[-positive]
      for(i in dropouts){
          cat("NNLS Step:\t Variable", i, "\tdropped\n")
        }
    }
  list(active = active[positive], R = R, beta = beta, positive = positive
       )
}
"plotCVLars" <-
function(cv.lars.object,se=TRUE){
  attach(cv.lars.object)
      plot(fraction, cv, type = "b", ylim = range(cv, cv + cv.error, 
                                     cv - cv.error))
    if(se)
      error.bars(fraction, cv + cv.error, cv - cv.error, 
                 width = 1/length(fraction))
  detach(cv.lars.object)
  
invisible()
}
"plot.lars" <-
  function(x, xvar=c("norm","df","arc.length"), breaks = TRUE, plottype = c("coefficients", "Cp"), 
           omit.zeros = TRUE, eps = 1e-10, ...)
{
  object <- x
  plottype <- match.arg(plottype)
  xvar <- match.arg(xvar)
  coef1 <- object$beta	### Get rid of many zero coefficients
  coef1 <- scale(coef1, FALSE, 1/object$normx)
  if(omit.zeros) {
    c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))
    nonzeros <- c1 > eps
    cnums <- seq(nonzeros)[nonzeros]
    coef1 <- coef1[, nonzeros]
  }
  else cnums <- seq(ncol(coef1))
  s1<-switch(xvar,
             norm={
               s1 <- apply(abs(coef1), 1, sum)
               s1/max(s1)
             },
             df=seq(length(object$arc.length)+1),
             arc.length=cumsum(c(0,object$arc.length))
             )
  xname<-switch(xvar,
                norm="|beta|/max|beta|",
                df="Df",
                arc.length="Arc Length"
                )
                
  if(plottype == "Cp") {
    Cp <- object$Cp
    plot(s1, Cp, type = "b", xlab=xname,main = object$type, ...)
  }
  else {
      matplot(s1, coef1, xlab = xname, ..., type = "b", 
              pch = "*", ylab = "Standardized Coefficients")
      title(object$type,line=2.5)
      abline(h = 0, lty = 3)
      axis(4, at = coef1[nrow(coef1),  ], label = paste(cnums
                                            ), cex = 0.80000000000000004, adj = 0)
      if(breaks) {
        axis(3, at = s1, labels = paste(seq(s1)-1),cex=.8)
        abline(v = s1)
      }

  }
  invisible()
}


"predict.lars" <-
  function(object, newx, s, type = c("fit", "coefficients"), mode = c("step", 
                                                               "fraction", "norm"), ...)
{
  mode <- match.arg(mode)
  type <- match.arg(type)
  if(missing(newx) & type == "fit") {
    warning("Type=fit with no newx argument; type switched to coefficients"
            )
    type <- "coefficients"
  }
  betas <- object$beta
  sbetas <- scale(betas, FALSE, 1/object$normx)	#scaled for unit norm x
  kp <- dim(betas)
  k <- kp[1]
  p <- kp[2]
  steps <- seq(k)
  if(missing(s)) {
    s <- steps
    mode <- "step"
  }
  sbeta <- switch(mode,
                  step = {
                    if(any(s < 0) | any(s > k))
                      stop("Argument s out of range")
                    steps
                  }
                  ,
                  fraction = {
                    if(any(s > 1) | any(s < 0))
                      stop("Argument s out of range")
                    nbeta <- drop(abs(sbetas) %*% rep(1, p))
                    nbeta/nbeta[k]
                  }
                  ,
                  norm = {
                    nbeta <- drop(abs(sbetas) %*% rep(1, p))
                    if(any(s > nbeta[k]) | any(s < 0))
                      stop("Argument s out of range")
                    nbeta
                  }
                  )
  sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
  sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
  usbeta<-unique(sbeta)
  useq<-match(usbeta,sbeta)
  sbeta<-sbeta[useq]
  betas<-betas[useq,]
  coord <- approx(sbeta, seq(sbeta), sfrac)$y
  left <- floor(coord)
  right <- ceiling(coord)
  newbetas <- ((sbeta[right] - sfrac) * betas[left,  , drop = FALSE] + (sfrac -
                                                         sbeta[left]) * betas[right,  , drop = FALSE])/(sbeta[right] - sbeta[
                                                                                          left])
  newbetas[left == right,  ] <- betas[left[left == right],  ]
  robject <- switch(type,
                    coefficients = list(s = s, fraction = sfrac, mode = mode, 
                      coefficients = drop(newbetas)),
                    fit = list(s = s, fraction = sfrac, mode = mode, fit = drop(
                                                                       scale(newx, object$meanx, FALSE) %*% t(newbetas)) + object$
                      mu))
  robject
}
"print.lars" <-
function(x, ...)
{
	cat("\nCall:\n")
	dput(x$call)
	cat("R-squared:", format(round(rev(x$R2)[1], 3)), "\n")
	actions <- x$actions
	jactions <- unlist(actions)
	jsteps <- rep(seq(along = actions), sapply(actions, length))
	actmat <- rbind(jsteps, jactions)
	vn <- names(jactions)
	if(is.null(vn))
		vn <- rep("", length(jactions))
	dimnames(actmat) <- list(c("Step", "Var"), vn)
	cat(paste("Sequence of", x$type, "moves:\n"))
	print(actmat[2:1,  ])
	invisible(x)
}
"updateR" <-
  function(xnew, R = NULL, xold, eps = .Machine$double.eps, Gram = FALSE)
{
###Gram argument determines the nature of xnew and xold
  xtx <- if(Gram) xnew else sum(xnew^2)
  norm.xnew <- sqrt(xtx)
  if(is.null(R)) {
    R <- matrix(norm.xnew, 1, 1)
    attr(R, "rank") <- 1
    return(R)
  }
  Xtx <- if(Gram) xold else drop(t(xnew) %*% xold)
  r <- backsolvet(R, Xtx)
  rpp <- norm.xnew^2 - sum(r^2)
  rank <- attr(R, "rank")	### check if R is machine singular
  if(rpp <= eps)
    rpp <- eps
  else {
    rpp <- sqrt(rpp)
    rank <- rank + 1
  }
  R <- cbind(rbind(R, 0), c(r, rpp))
  attr(R, "rank") <- rank
  R
}
"backsolvet"<-
function(r, x, k=ncol(r))
{
  backsolve(r,x,k,transpose=TRUE)
}
"delcol" <-
function(r, z, k = p)
{
	p <- dim(r)[1]
	r <- r[,  - k, drop = FALSE]
	z <- as.matrix(z)
	dz <- dim(z)
	storage.mode(r) <- storage.mode(z) <- "double"
	.Fortran("delcol",
		r,
		as.integer(p),
		as.integer(k),
		z,
		as.integer(dz[1]),
		as.integer(dz[2]),
                PACKAGE="lars")[c(1, 4)]
}
".First.lib" <-
function (lib, pkg) 
  library.dynam("lars", pkg, lib)
