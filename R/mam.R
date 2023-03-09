### Marginal Additive Model ###

#' Fit a Marginal Additive Model
#'
#' Fit a Marginal Additive Model. This function implements Algorithms ... in ...
#' See details.
#'
#' @param smooth A list of objects of class \code{xx.smooth.spec} created using \code{mgcv::s()}.
#' Any smooth that is compatible with \code{gamm4::gamm()} is compatible; see the documentation for that
#' function. To our knowledge the main (only?) incompatible smooth is the tensor product \code{bs='te'},
#' and you can use \code{bs='t2'} in that case.
#' @param re A **two-sided** random effects formula of the form \code{y~(1|id)} compatible with use in \code{lme4::lmer}.
#' See the documentation for that function, or \code{lme4::lFormula}. This formula should not contain fixed effects,
#' only grouping terms and the response. **This is where the response is specified**. If you have no group-level random effects,
#' then you don't need \code{mam}, use \code{mgcv::gam}.
#' @param fe A **one-sided** fixed effects formula of the form \code{~x} to be passed to \code{model.matrix()}.
#' @param dat The data to be used to fit the conditional model. A \code{data.frame} containing all of the
#' variables in \code{smooth}, \code{re}, and \code{fe}, including the response (specified through \code{re}).
#' @param margdat A \code{data.frame} whose rows are the points where the approximate marginal means will be computed.
#' Defaults to \code{dat}, and we have so far found no reason not to do this, but you can do something different if you want.
#' @param preddat Data at which to compute fitted values. Also defaults to \code{dat}, for now.
#' @param control Output of \code{mam::mam_control()}. See documentation for that function for options.
#' @param ... not used.
#'
#' @returns Object of classes \code{c('mam','gam','glm','lm')} representing the marginal additive model.
#'
#' @export
#'
#' @examples
#' ## Not run
#'
mam <- function(smooth,re,fe,dat,margdat=dat,preddat=dat,est.fun = FALSE,control=mam_control(),...) {
  ## SETUP CONTROL ##
  method <- control$method[1]
  varmethod <- control$varmethod[1]
  absorbcons <- TRUE # Do NOT change this.
  verbose <- control$verbose
  centered <- control$centered
  k <- control$k
  retcond <- control$retcond
  ## END SETUP CONTROL ##

  ## MODEL SETUP ##
  # 1: SMOOTHS #
  Xr <- Xf <- newsparsemat(nrow(dat),0)
  Xrmarg <- Xfmarg <- newsparsemat(nrow(margdat),0)
  Xrpred <- Xfpred <- newsparsemat(nrow(preddat),0)
  numsmooth <- 0
  r <- 0
  if (!is.null(smooth)) {
    tm <- Sys.time()
    if (verbose) cat("Constructing smooths... ")
    if (!inherits(smooth,'list')) smooth <- list(smooth)
    # Conditional model
    SS <- lapply(lapply(smooth,mgcv::smoothCon,data=dat,absorb.cons = absorbcons),'[[',1)
    XXmarg <- Reduce(cbind,lapply(SS,mgcv::PredictMat,data = margdat))
    XXpred <- Reduce(cbind,lapply(SS,mgcv::PredictMat,data = preddat))
    numsmooth <- length(smooth) # Number of smooth terms

    EE <- lapply(lapply(lapply(SS,'[[','S'),'[[',1),eigen)

    p <- sapply(lapply(EE,'[[','vectors'),ncol)
    r <- sapply(lapply(EE,'[[','values'),function(x) sum(x>.Machine$double.eps))
    m <- p-r
    URlist <- mapply(function(x,y) x[ ,1:y],lapply(EE,'[[','vectors'),r,SIMPLIFY = FALSE)
    UFlist <- mapply(
      function(x,y,z) {
        if (y<z) return(x[ ,(1+y):z])
        newsparsemat(z,z)
      },lapply(EE,'[[','vectors'),r,p,SIMPLIFY = FALSE)
    URlist <- lapply(URlist,cbind) # Ensure they stay matrices
    UFlist <- lapply(UFlist,cbind) # Ensure they stay matrices

    UR <- Matrix::bdiag(URlist)
    UF <- Matrix::bdiag(UFlist)
    # if m=1 UF gets coerced to numeric
    if (!is.matrix(UF)) UF <- cbind(UF)

    Dpi <- Matrix::Diagonal(sum(r),1 / sqrt(Reduce(c,lapply(lapply(EE,'[[','values'),function(x) x[x>.Machine$double.eps]))))

    Xlist <- lapply(SS,'[[','X')
    X <- Reduce(cbind,Xlist)

    Xr <- as.matrix(X %*% UR %*% Dpi)
    Xf <- as.matrix(X %*% UF)
    dups <- !base::duplicated(t(Xf)) & apply(Xf,2,function(x) !all(x==0)) # Remove the duplicated intercepts
    if (length(dups) > 1) Xf <- Xf[ ,which(dups)]

    # Marginal and prediction
    # Use the SAME Ur and Uf
    Xrmarg <- as.matrix(XXmarg %*% UR %*% Dpi)
    Xrpred <- as.matrix(XXpred %*% UR %*% Dpi)
    Xfmarg <- as.matrix(XXmarg %*% UF)
    Xfpred <- as.matrix(XXpred %*% UF)
    # Same dups
    if (length(dups) > 1) {
      Xfmarg <- Xfmarg[ ,which(dups)]
      Xfpred <- Xfpred[ ,which(dups)]
    }

    dt <- difftime(Sys.time(),tm,units = 'secs')
    if (verbose) cat("finished, took",round(dt,4),"seconds.\n")
  }
  # END SMOOTHS #

  # 2: FIXED EFFECTS #
  if(!is.null(fe)){
    tm <- Sys.time()
    if (verbose) cat("Constructing linear effects... ")
    Xlin <- stats::model.matrix(fe,data=dat)[,-1,drop=F]
    Xf <- cbind(Xlin,Xf)
    Xfmarg <- cbind(stats::model.matrix(fe,data=margdat)[,-1],Xfmarg) # same for marg
    Xfpred <- cbind(stats::model.matrix(fe,data=preddat)[,-1],Xfpred) # same for pred
    dt <- difftime(Sys.time(),tm,units = 'secs')
    if (verbose) cat("finished, took",round(dt,4),"seconds.\n")
  }

  # Finally, add an intercept if not already present
  if (!any(apply(Xf,2,function(x) all(x==1)))) {
    Xf <- cbind(1,Xf)
    Xfmarg <- cbind(1,Xfmarg)
    Xfpred <- cbind(1,Xfpred)
  }
  # END FIXED EFFECTS #

  ## END MODEL SETUP ##

  ## CONDITIONAL MODEL ##

  # Fit with TMB
  tm <- Sys.time()
  if (verbose) cat("Fitting conditional model... ")
  reform <- lme4::lFormula(re,data=dat)
  redim <- Reduce(sum,lapply(reform$reTrms$cnms,length))
  LamforCond <- methods::as(Matrix::t(reform$reTrms$Lambdat)[1:redim,1:redim],'TsparseMatrix')
  tmbdat <- list(
    XF = as.matrix(Xf),
    XR = as.matrix(Xr),
    # A = methods::as(Matrix::t(reform$reTrms$Zt),'dgTMatrix'),
    # Lam = methods::as(Matrix::t(reform$reTrms$Lambdat),'dgTMatrix'),
    A = methods::as(Matrix::t(reform$reTrms$Zt),'TsparseMatrix'),
    Lam = LamforCond,

    Lind = as.integer(reform$reTrms$Lind-1)[1:length(LamforCond@x)],
    diagind = as.integer(reform$reTrms$theta), # Diagonals initialized to 1, off-diags to 0.
    y = stats::model.frame(re,dat)[,1], # Response
    p = as.integer(numsmooth), # Number of smooth terms
    r = r # Rank of each smooth
  )
  tmbparams <- with(tmbdat,list(
    betaF = rep(0,ncol(XF)),
    bR = rep(0,ncol(XR)),
    theta = rep(0,length(diagind)),
    logsmoothing = rep(0,p),
    U = rep(0,ncol(A))
  ))
  rand <- c('bR','U','betaF') # "REML", or whatever it is, integrates over betaF too

  template <- TMB::MakeADFun(
    data = c(model = "conditional_model",tmbdat),
    parameters = tmbparams,
    random = rand,
    silent = TRUE,
    DLL = 'mam_TMBExports'
  )
  X <- with(tmbdat,cbind(XF,XR))
  Xmarg <- cbind(Xfmarg,Xrmarg)
  Xpred <- cbind(Xfpred,Xrpred)
  if(centered) Xpred <- sweep(Xpred,2,colMeans(Xpred),'-') # Return centered smooths.

  if(est.fun) Xpred[,apply(Xpred, MARGIN = 2, sd) < 1e-3] <- 0

  if (method == 'BFGS') {
    opt <- with(template,stats::optim(par,fn,gr,method='BFGS',hessian=TRUE))
    # Rename output to match trustOptim
    opt$solution <- opt$par
  } else if (method == 'trust') {
    utils::capture.output(opt <- with(template,trustOptim::trust.optim(
      x = par,
      fn = fn,
      gr = gr,
      hs = function(x) methods::as(methods::as(numDeriv::jacobian(gr,x),'generalMatrix'),'CsparseMatrix'),
      method = 'Sparse'
    )))
  } else {
    stop(paste0("Unrecognized optimization method:",method,"\n"))
  }

  # Point estimates
  tmbcoefs <- with(template$env,last.par[random])
  tmbbetaF <- tmbcoefs[names(tmbcoefs)=='betaF']
  tmbbR <- tmbcoefs[names(tmbcoefs)=='bR']
  tmbU <- tmbcoefs[names(tmbcoefs)=='U']

  condcoefs <- c(tmbbetaF,tmbbR)
  condest <- Xpred %*% condcoefs
  # separate linear and smooth coefficients
  if(!is.null(fe)){
    condlincoefs <- condcoefs[1:ncol(Xlin)]
    condsmoothcoefs <- condcoefs[(ncol(Xlin)+1):length(condcoefs)]
  }else{
    condlincoefs <- NULL
    condsmoothcoefs <- condcoefs
  }
  # compute variances & SEs
  if(retcond==TRUE){### requires potentially slow inversion (can speed up later and remove if statement)

    condprecmat <- TMB::sdreport(template,getJointPrecision = TRUE)$jointPrecision
    condidx <- which(rownames(condprecmat) %in% c('betaF','bR'))

    # OLD:
    # condvarmat <- Matrix::solve(condprecmat)[condidx,condidx] ## saving in order to output ## this is potentially slow
    # condvar2 <- diag(as.matrix(Xpred %*% condvarmat %*% t(Xpred))) ##old: diag(Xpred %*% solve(condprecmat)[condidx,condidx] %*% t(Xpred))

    # NEW:
    Cchol <- tryCatch(Matrix::Cholesky(condprecmat,LDL = FALSE,perm = FALSE),error = function(e) e,warning = function(w) w)
    if (inherits(Cchol,'condition')) {
      condvar <- NULL
      retcond <- FALSE # for later
    } else {
      newX <- newsparsemat(ncol(condprecmat),nrow(Xpred))
      newX[condidx, ] <- t(Xpred)
      condvar <- Matrix::colSums(Matrix::solve(Cchol,newX,system="L")^2)
    }

    }else{
    condvar <- condvarmat <- NULL
  }

  dt <- difftime(Sys.time(),tm,units = 'secs')
  if (verbose) cat("finished, took",round(dt,4),"seconds.\n")
  ## END CONDITIONAL MODEL ##

  #############################################################################
  ## STEP 2: Compute the marginal means ##

  tm <- Sys.time()
  if (verbose) cat("Approximating marginal means... ")
  ## Get the marginal means now ##
  thetaest <- opt$solution[names(template$par) == 'theta'] # Trust region
  loglambdaest <- opt$solution[names(template$par) == 'logsmoothing'] # Trust region
  gpix <- numeric(nrow(Xfmarg))
  paramlist <- list(
    betaF = tmbbetaF,bR = tmbbR,
    theta = thetaest,loglambda = loglambdaest,
    U_hat = tmbcoefs[names(tmbcoefs)=='U']
  )

  paramvec <- Reduce(c,paramlist)
  # J <- newsparsemat(length(gpix),length(paramvec))

  # GHQ grid
  renum <- length(lapply(reform$reTrms$cnms,length))

  ghqgrid <- mvQuad::createNIGrid(dim = redim,type = 'GHe',level = k)
  ghqnodes <- as.matrix(mvQuad::getNodes(ghqgrid)); class(ghqnodes) <- 'matrix'
  ghqweights <- as.numeric(mvQuad::getWeights(ghqgrid))
  if (ncol(ghqnodes) == 1) {
    ghqweights <- ghqweights * as.numeric(stats::dnorm(ghqnodes))
  } else {
    ghqweights <- ghqweights * apply(apply(ghqnodes,1,stats::dnorm),2,prod)
  }

  # Get marginal RE design matrix
  reformmarg <- lme4::lFormula(re,data=margdat)
  Zmarg <- Matrix::t(reformmarg$reTrms$Zt)
  # TODO: unequal group sizes
  Zflat <- flatten_bdmat(Zmarg,unique(table(reformmarg$reTrms$flist$id)),redim)

  # Get a single subject's RE lambda and Lind and diagind
  LamforGHQ <- tmbdat$Lam[1:redim,1:redim]
  if (is.numeric(LamforGHQ)) {
    # Guard against 1d case where it converts to numeric
    LamforGHQ <- as.matrix(LamforGHQ)
    # LamforGHQ <- methods::as(LamforGHQ,'dgTMatrix')
    LamforGHQ <- methods::as(methods::as(methods::as(Matrix::Matrix(LamforGHQ), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  }
  diagindMarg <- tmbdat$diagind # This one is the same
  # The Lind is NOT the same. TODO: unit test this for correlated and uncorrelated int/slope
  LindMarg <- seq(0,max(tmbdat$Lind))
  # templatelist <- vector(mode='list',length=length(gpix))

  # NEW: single template
  tmbdatamarg <- list(
    XF=Xfmarg,XR=Xrmarg,Z=Zflat,
    Lam = LamforGHQ,Lind = LindMarg,diagind = diagindMarg,m=redim,
    Q = ghqnodes,w = ghqweights
  )
  template_marg <- TMB::MakeADFun(
    data = c(model = "marginal_mean2",tmbdatamarg),
    parameters = paramlist,
    silent = TRUE,
    DLL = "mam_TMBExports",
    ADreport = TRUE
  )
  gpix <- template_marg$fn(paramvec)
  J <- template_marg$gr(paramvec)

  # for (i in 1:length(gpix)) {
  #   # Get the vector of random effects COVARIATES
  #   # NOTE: in the single group case, this should be (1,z) where z is the variable that there is a
  #   # random slope for.
  #   # Problem: when z = 0, lme4 doesn't store it as a zero, it stores it as an empty space
  #   # So the dimension of Zmarg[i, ] is less than redim, and this causes an error
  #   # in TMB due to mismatched dimensions
  #   # So check the dimension of v = Zmarg[i, ]. If it's less than redim, pad it
  #   # with zeros according to renum.
  #   # THIS HAS NOT BEEN THOROUGHLY TESTED, and PROBABLY DOES NOT WORK for multiple
  #   # grouping factors
  #   # TODO: test this for multiple grouping factors
  #   # TODO: implement a hack where you just add a small number to data zeroes avoid this.
  #   if (verbose) {
  #     if ((i-1)%%10==0) {
  #       cat(paste0("\nMarginal mean: i = ",i))
  #     } else {
  #       cat(",",i)
  #     }
  #   }
  #   v <- as.numeric(methods::as(Zmarg[i, ],'sparseVector')@x)
  #   if (length(v) < redim) {
  #     # Add a zero in the second place
  #     v <- splice(v,0,2)
  #   }
  #   tmbdatamarg <- list(
  #     # link = 1L,
  #     xf = Xfmarg[i, ],xr = Xrmarg[i, ],
  #     v = v,
  #     Lam = LamforGHQ,Lind = LindMarg,diagind = diagindMarg,
  #     Q = ghqnodes,w = ghqweights
  #   )
  #   tmp <- TMB::MakeADFun(
  #     data = c(model = "marginal_mean",tmbdatamarg),
  #     parameters = paramlist,
  #     silent = TRUE,
  #     DLL = "mam_TMBExports"
  #   )
  #   gpix[i] <- tmp$fn(paramvec) # NOTE: logit marg prob
  #   J[i, ] <- methods::as(tmp$gr(paramvec),'sparseVector')
  # }
  dt <- difftime(Sys.time(),tm,units = 'secs')
  if (verbose) cat("finished, took",round(dt,4),"seconds.\n")
  tm <- Sys.time()
  if (verbose) cat("Computing variance factor... ")
  # Joint precision from the conditional model, ordered to match the order of the params in the marginal TMB template
  sdr <- TMB::sdreport(template,getJointPrecision = TRUE)
  H <- sdr$jointPrecision
  # Compute the variance one of three ways
  if (varmethod == 0) {
    # 0 == "auto", so choose based on whether any dd < 0
    L <- Matrix::Cholesky(H,perm=TRUE,LDL=TRUE)
    D <- Matrix::solve(L,system="D")
    dd <- D@x
    varmethod <- 2 # Full Hessian
    if (any(dd <0) ) varmethod <- 1 # Kass & Steffey
  }
  if (varmethod == 1) {
    # Kass & Steffy. Compute the two variance factors
    fixparam <- names(template$par)
    fullnames <- colnames(H)
    condidx <- (1:ncol(H))[-which(fullnames %in% fixparam)]
    # Compute the "conditional" variance factor
    Hcond <- H[condidx,condidx]
    Lcond <- Matrix::Cholesky(Hcond,perm=TRUE,LDL=TRUE)
    Dcond <- Matrix::solve(Lcond,system="D")
    Jcond <- J[ ,condidx]
    Vt <- sqrt(Dcond) %*% Matrix::solve(Lcond,Matrix::solve(Lcond,t(Jcond),system = "P"),system = "L") ## v transpose

    # Compute the "marginal" variance factor
    margidx <- (1:ncol(H))[which(fullnames %in% fixparam)]
    Hmarg <- methods::as(methods::as(Matrix::Matrix(opt$hessian),'generalMatrix'),'CsparseMatrix')
    Lmarg <- Matrix::Cholesky(Hmarg,perm=TRUE,LDL=TRUE)
    Dmarg <- Matrix::solve(Lmarg,system="D")
    Jmarg <- J[ ,margidx]
    Vtmarg <- sqrt(Dmarg) %*% Matrix::solve(Lmarg,Matrix::solve(Lmarg,t(Jmarg),system = "P"),system = "L") ## v transpose
  } else if (varmethod == 2) {
    # Full Hessian
    L <- Matrix::Cholesky(H,perm=TRUE,LDL=TRUE)
    D <- Matrix::solve(L,system="D")
    Vt <- sqrt(D) %*% Matrix::solve(L,Matrix::solve(L,t(J),system = "P"),system = "L") ## v transpose
    # Get a 0 Vtmarg, for later
    margdim <- length(template$par)
    Vtmarg <- newsparsemat(margdim,nrow(J))
  }
  dt <- difftime(Sys.time(),tm,units = 'secs')
  if (verbose) cat("finished, took",round(dt,4),"seconds.\n")
  ## END MARGINAL MEANS ##

  ## MARGINAL MODEL ##

  tm <- Sys.time()
  if (verbose) cat("Fitting MAM... ")
  ## Fit the MAM ##
  # Coefficients
  XtXmarg <- Matrix::crossprod(Xmarg)
  XtXmarginv <- solve(XtXmarg)
  mamcoef <- Matrix::solve(XtXmarg,Matrix::crossprod(Xmarg,gpix))

  # Fitted values
  mamest <- Xpred %*% mamcoef

  # Remove zero rows; don't affect the colsums and vastly improves speed
  mamvarfactor_cond <- Vt[Matrix::rowSums(Vt)!=0, ] %*% Xmarg %*% XtXmarginv
  mamvarfactor_marg <- Vtmarg[Matrix::rowSums(Vtmarg)!=0, ] %*% Xmarg %*% XtXmarginv

  mamcoefvar <- Matrix::colSums((mamvarfactor_cond)^2) + Matrix::colSums((mamvarfactor_marg)^2)
  mamestvar <- Matrix::colSums((mamvarfactor_cond %*% t(Xpred))^2) + Matrix::colSums((mamvarfactor_marg %*% t(Xpred))^2)

  ## separate linear and smooth coefficients
  if(!is.null(fe)){
    mamlincoefs <- mamcoef[1:ncol(Xlin)]
    names(mamlincoefs) <- colnames(Xlin)
    mamsmoothcoefs <- mamcoef[(ncol(Xlin)+1):length(mamcoef)]
  }else{
    mamlincoefs <- NULL
    mamsmoothcoefs <- mamcoef
  }
  # Get names for smooth coefs
  if (!is.null(smooth)) {
    smoothnames <- Reduce(c,lapply(lapply(SS,'[[','term'),paste0,collapse='+'))
    smoothdims <- Reduce(c,lapply(SS,'[[','bs.dim'))
  }

  mamvarmat <- NULL

  ## exponentiate diagonal terms of theta
  thetaest[as.integer(reform$reTrms$theta)==1] <- exp(thetaest[as.integer(reform$reTrms$theta)==1])
  LamforGHQ@x[] = thetaest[LindMarg+1] # zero indexing
  Sig <- LamforGHQ%*%Matrix::t(LamforGHQ)

  dt <- difftime(Sys.time(),tm,units = 'secs')
  if (verbose) cat("finished, took",round(dt,4),"seconds.\n")

  ## END MARGINAL MODEL ##

  ## Return results
  mam <- list(fitted = mamest,fitted_se = sqrt(mamestvar),coef_se = sqrt(mamcoefvar),
              coeflin = mamlincoefs, coefsmooth = mamsmoothcoefs,logsmoothing = loglambdaest,Xpred = Xpred)
  variance <- list(sigma = Sig,
                   theta = thetaest,
                   lambda=exp(loglambdaest),
                   mamvarfactor_cond=mamvarfactor_cond,
                   mamvarfactor_marg=mamvarfactor_marg,
                   optresults = opt)
  conditional <- list(predU = matrix(tmbU,nrow=length(unique(dat$id)),byrow=TRUE))
  marginal <- list(prob = ilogit(gpix))
  if (retcond) {
    conditional <- c(conditional,list(fitted = as.numeric(condest),fitted_se = sqrt(as.numeric(condvar)),
                        coeflin = condlincoefs, coefsmooth = condsmoothcoefs,
                        # var = condvarmat,
                        theta=Sig, lambda=exp(loglambdaest),
                        predU = matrix(tmbU,nrow=length(unique(dat$id)),byrow=TRUE) ))
  }
  list(mam = mam,conditional = conditional,marginal=marginal,variance = variance)
}
