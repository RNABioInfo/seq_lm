# Stable NB profile likelihood, including extreme fractional low counts.
# Reparameterized intercept avoids underflow in edgeR's one-way null fits.
softplus <- function(z) pmax(z,0)+log1p(exp(-abs(z)))
stable_null_deviance <- function(y,off,beta,target,phi) {
    eta <- off+outer(beta,target);eta <- eta-apply(eta,1,max)
    total <- rowSums(y);n <- nrow(y);q <- ncol(y)
    poisson <- rep_len(phi,n)==0
    phi <- pmax(rep_len(phi,n),1e-12)
    inv <- 1/phi;logphi <- log(phi)
    lo <- log(pmax(total,1e-300))-log(q)-40
    hi <- log(pmax(total,1e-300))-apply(eta,1,min)+40
    alpha <- log(pmax(total,1e-300))-log(rowSums(exp(eta)))
    for(k in seq_len(60)) {
        z <- eta+alpha+logphi
        prob <- plogis(z);comp <- plogis(-z)
        score <- rowSums(y*comp-prob*inv)
        deriv <- -rowSums((y+inv)*prob*comp)
        lo[score>0] <- alpha[score>0];hi[score<=0] <- alpha[score<=0]
        nxt <- alpha-score/deriv
        bad <- !is.finite(nxt)|nxt<=lo|nxt>=hi
        nxt[bad] <- (lo[bad]+hi[bad])/2
        converged <- abs(score)<=1e-10*pmax(1,total)
        converged[poisson] <- TRUE
        nxt[converged] <- alpha[converged]
        alpha <- nxt
        if(all(converged)) break
    }
    z <- eta+alpha+logphi
    first <- y*(-log1p(inv/y)+softplus(-z));first[y==0] <- 0
    result <- 2*rowSums(first+inv*(-log1p(phi*y)+softplus(z)))
    if(any(poisson)) {
        lm <- (eta+alpha)[poisson,,drop=FALSE]; yy <- y[poisson,,drop=FALSE]
        term <- yy*(log(yy)-lm);term[yy==0]<-0
        result[poisson] <- 2*rowSums(term-yy+exp(lm))
    }
    pmax(0,result)
}

# Invert the installed edgeR QL single-contrast test, conditional on its fitted
# offsets, effective NB dispersion, moderated s2.post and denominator df.
# Specialized and validated for a two-group, no-covariate design.
profile_intervals <- function(fit, target, level=.95, iterations=30L) {
    stopifnot(length(target)==ncol(fit$counts), all(target %in% c(0,1)),
              ncol(fit$design)==2L, length(unique(target))==2L)
    y <- fit$counts; n <- nrow(y)
    phi <- rep_len(fit$dispersion,n)
    if (!is.null(fit$average.ql.dispersion)) phi <- phi/fit$average.ql.dispersion
    df <- if(is.null(fit$df.residual.zeros)) fit$df.residual.adj else fit$df.residual.zeros
    df <- pmin(rep_len(fit$df.prior,n)+df,sum(fit$df.residual))
    critical <- rep_len(fit$s2.post,n)*qf(level,1,df)
    off <- as.matrix(fit$offset)+matrix(0,n,ncol(y))
    base <- fit$deviance
    reference <- which(fit$design[which(target==0)[1],]==1)
    treatment <- which(fit$design[which(target==1)[1],]==1)
    stopifnot(length(reference)==1,length(treatment)==1)
    co <- fit$unshrunk.coefficients
    if(is.null(co)) co <- fit$coefficients
    center <- drop(co[,treatment]-co[,reference])
    zero0 <- rowSums(y[,target==0,drop=FALSE])==0
    zero1 <- rowSums(y[,target==1,drop=FALSE])==0
    boundary <- zero0|zero1
    # Anchor close to the boundary maximum without using divergent coefficients.
    anchor <- center; anchor[zero0] <- 40; anchor[zero1] <- -40
    deviance_at <- function(beta) {
        pmax(0, stable_null_deviance(y,off,beta,target,phi)-base)
    }
    mu <- fit$fitted.values
    w <- mu/(1+phi*mu)
    if(!is.null(fit$weights)) w <- w*fit$weights
    v <- 1/rowSums(w[,target==0,drop=FALSE])+1/rowSums(w[,target==1,drop=FALSE])
    startstep <- pmax(.05,pmin(10,sqrt(v*critical)))
    startstep[!is.finite(startstep)] <- 1
    answers <- list(); errors <- list()
    for (direction in c(-1,1)) {
        unbounded <- if(direction<0) zero1 else zero0
        low <- rep(0,n); high <- startstep
        for(k in seq_len(32)) {
            value <- deviance_at(anchor+direction*high)
            expand <- value<critical & !unbounded
            if(!any(expand)) break
            high[expand] <- high[expand]*2
        }
        bracketed <- deviance_at(anchor+direction*high)>=critical
        for(k in seq_len(iterations)) {
            mid <- (low+high)/2
            value <- deviance_at(anchor+direction*mid)
            inside <- value<critical
            low[inside] <- mid[inside];high[!inside] <- mid[!inside]
        }
        endpoint <- anchor+direction*(low+high)/2
        residual <- abs(deviance_at(endpoint)-critical)/pmax(critical,1e-12)
        endpoint[unbounded] <- direction*Inf
        endpoint[!bracketed&!unbounded] <- NA_real_
        residual[unbounded|!bracketed] <- NA_real_
        answers[[length(answers)+1]] <- endpoint/log(2)
        errors[[length(errors)+1]] <- residual
    }
    low <- answers[[1]]; high <- answers[[2]];center <- center/log(2)
    radius <- pmax(center-low,high-center)
    radius[boundary] <- Inf
    radius[!is.finite(critical)|df<=0|critical<=0] <- NA_real_
    list(lower=low,upper=high,center=center,R=radius,df_total=df,
         critical_deviance=critical,effective_phi=phi,
         curvature_se=sqrt(v*rep_len(fit$s2.post,n))/log(2),
         lower_root_relative_error=errors[[1]],upper_root_relative_error=errors[[2]],
         boundary=boundary)
}
