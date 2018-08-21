## squarem2 function, from the SQUAREM R package, by RAVI VARADHAN and CHRISTOPHE ROLAND Scandinavian Journal of Statistics, Vol. 35: 335–353, 2008
## link to cran webpage for SQUAREM R package:
## https://cran.r-project.org/web/packages/SQUAREM/index.html

squarem2<-function (par, fixptfn, ..., control = list()) 
{
    control.default <- list(K = 1, square = TRUE, method = 3, 
        step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1, 
        tol = 1e-07, maxiter = 1500, trace = FALSE)
    namc <- names(control)
    if (!all(namc %in% names(control.default))) 
        stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
    ctrl <- modifyList(control.default, control)
    method <- ctrl$method
    maxiter <- ctrl$maxiter
    tol <- ctrl$tol
    kr <- ctrl$kr
    step.min <- ctrl$step.min0
    step.max0 <- ctrl$step.max0
    step.max <- ctrl$step.max0
    mstep <- ctrl$mstep
    trace <- ctrl$trace
    if (trace) 
        cat("Squarem-2 \n")
    iter <- 1
    feval <- 0
    kount <- 0
    conv <- TRUE
    while (feval < maxiter) {
        extrap <- TRUE
        p1 <- try(fixptfn(par, ...), silent = TRUE)
        feval <- feval + 1
        if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) 
            break
        q1 <- p1 - par
        sr2 <- crossprod(q1)
        if (sqrt(sr2) < tol) 
            break
        p2 <- try(fixptfn(p1, ...), silent = TRUE)
        feval <- feval + 1
        if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) 
            break
        q2 <- p2 - p1
        sq2 <- sqrt(crossprod(q2))
        res <- sq2
        if (sq2 < tol) 
            break
        sv2 <- crossprod(q2 - q1)
        srv <- crossprod(q1, q2 - q1)
        alpha <- switch(method, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
        alpha <- max(step.min, min(step.max, alpha))
        p.new <- par + 2 * alpha * q1 + alpha^2 * (q2 - q1)
        if (abs(alpha - 1) > 0.01) {
            ptmp <- try(fixptfn(p.new, ...), silent = TRUE)
            feval <- feval + 1
            if (class(ptmp) == "try-error" | any(is.nan(unlist(ptmp)))) {
                p.new <- p2
                if (alpha == step.max) 
                  step.max <- max(step.max0, step.max/mstep)
                alpha <- 1
                extrap <- FALSE
            }
            else {
                res <- sqrt(crossprod(ptmp - p.new))
                parnorm <- sqrt(crossprod(p2)/length(p2))
                kres <- kr * (1 + parnorm) + sq2
                p.new <- if (res <= kres) 
                  ptmp
                else p2
                if (res > kres) {
                  if (alpha == step.max) 
                    step.max <- max(step.max0, step.max/mstep)
                  alpha <- 1
                  extrap <- FALSE
                }
            }
        }
        if (alpha == step.max) 
            step.max <- mstep * step.max
        if (step.min < 0 & alpha == step.min) 
            step.min <- mstep * step.min
        if (trace) 
            cat("Residual: ", res, "  Extrapolation: ", extrap, 
                "  Steplength: ", alpha, "\n")
        par <- p.new
        iter <- iter + 1
    }
    if (feval >= maxiter) 
        conv <- FALSE
    return(list(par = par, value.objfn = NA, iter = iter, fpevals = feval, 
        objfevals = 0, convergence = conv))
}
