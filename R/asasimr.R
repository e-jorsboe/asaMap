## First source asasimr.R

#' Simulate ancestry specific assocation
#'
#' @param qvec individual admixture proportions from population 1 (n vec)
#' @param pvec allele frequencies in the two ancestral populations (2 vec)
#' @param bvec population specific effects (2 vec for add and 3 vec for rec) or (4 long if you want to simulate effect for un-counted allele)
#' @param model either "add" or "rec"
#' @param cvec intercept and effects of covariates if any (p vec)
#' @param cov design matrix with intercept (p vec)
#' @param quant simulate quantitative trait, logistic regression model if FALSE
#' @return list of genotype vector (gt[n]), phenotype vector (y[n]), ancestral state (as[2,n]), allelic genotypes (ag[2,n]) and design for simulated population specific effects (x [2,n] for additive and [3,n] for recessive)
#' @export
#' @examples
#' set.seed(5)
#'
#' N <- 2500
#' qdist <- c(0,0.25,0.5,0.75,1)
#' Q <- rep(qdist,each=N/length(qdist))
#' mafs <- c(0.2,0.2)
#'
#' # Simulate additive model, only effect in pop1:
#' effs <- c(0.5,0) # length 2 for additive
#' simA <- asaSim(qvec=Q,pvec=mafs,bvec=effs,model="add",cvec=0)
#' simAL <- asaSim(qvec=Q,pvec=mafs,bvec=effs,model="add",cvec=0,quant=FALSE)
#'
#' # Simulate recessive model, only effect in pop2:
#' effs <- c(0.5,0,0) # length 3 for recessive
#' simR <- asaSim(qvec=Q,pvec=mafs,bvec=effs,model="rec",cvec=0)
#' simRL <- asaSim(qvec=Q,pvec=mafs,bvec=effs,model="rec",cvec=0,quant=FALSE)
#'
#' # Simulate additive model, effect in pop1 of alternative allele:
#' effs1 <- c(0,0,0.5,0) # length 4 for additive
#' simA1 <- asaSim(qvec=Q,pvec=mafs,bvec=effs1,model="add",cvec=0)
#' simAL1 <- asaSim(qvec=Q,pvec=mafs,bvec=effs1,model="add",cvec=0,quant=FALSE)
#' 
asaSim <- function(qvec, pvec, bvec, model, cvec, cov, quant=TRUE,wrongAllele=F,setSeed=NULL){
    
    if(!is.null(setSeed)){
        set.seed(setSeed)
    }
    
    nInd <- length(qvec)
    ## add intersect
    if(missing(cov)){
        cov <- matrix(1, nrow=nInd, ncol=1)
    }
    
    qmat <- rbind(qvec,1-qvec)
    ## Allelic state (1=pop1, 2=pop2)
    as<-apply(qmat, 2, function(x)sample(1:2,2,replace=TRUE,prob=x)) ## 2 x nInd
    ## Sample allelic genotype (1=risk allele)
    ag <- matrix(rbinom(2*nInd,1,pvec[as]),nrow=2,byrow=FALSE) ## 2 x nInd

    ## checking that bvec is either 2 or 4 long for the additive model
    if(model=="add" & length(bvec)!=2 & length(bvec)!=4){
        stop("bvec must have length of 2 or 4")
    }

    ## checking that bvec is either 3 or 5 long for the recessive model
    if(model=="rec" & length(bvec)!=3 & length(bvec)!=5){
        stop("bvec must have length of 3 or 5")
    }
    
    ## simulating states for the additive model
    if(model=="add") {
        x1 <- as.integer(ag[1,]*(as[1,]==1)+ag[2,]*(as[2,]==1))
        x2 <- as.integer(ag[1,]*(as[1,]==2)+ag[2,]*(as[2,]==2))
        x3 <- as.integer((1-ag[1,])*(as[1,]==1)+(1-ag[2,])*(as[2,]==1))  
        x <- cbind(x1,x2,x3)
        ## simulating - if bvec is longer than 2 also simulates alternate alleles    
        if(length(bvec)>3){
            x4 <- as.integer((1-ag[1,])*(as[1,]==2)+(1-ag[2,])*(as[2,]==2))
            x<-cbind(x,x4)
        }
    ## simulating states for the recessive model
    } else if(model=="rec") {
        x1 <- as.integer((colSums(as) == 2) * (colSums(ag) == 2))
        x2 <- as.integer((colSums(as) == 4) * (colSums(ag) == 2))
        xM <- as.integer((colSums(as) == 3) * (colSums(ag) == 2))
        x <- cbind(x1,x2,xM)
        if(length(bvec)>3){
          x3 <- as.integer((colSums(as) == 2) * (colSums(ag) == 0))
          x4 <- as.integer((colSums(as) == 4) * (colSums(ag) == 0))
          x <- cbind(x,x3,x4)
        }
        
    } else {
        stop('Model must be either "add" or "rec", please.')
    }

    
    if(model=="add"){
        if(length(bvec)==4){
            ## also incorporating effects from both alternative allele (pop 1 and 2)
            evec <- x %*% bvec + cov %*% cvec
        } else{
            ## only incorporating effect from primary allele (x1 and x2)
            evec <- x[,1:2] %*% bvec + cov %*% cvec
        }
    } else{
        if(length(bvec)==5){
            ## also incorporating effects from both recessive alternative allele (pop 1 and 2)
            evec <- x %*% bvec + cov %*% cvec
        } else{
            ## only incorporating effect from recessive primary allele (x1, x2 and xM)
            evec <- x[,1:3] %*% bvec + cov %*% cvec
        }
    }
    ## Simulated responce
    ## note that standard deviation is 1
    if(quant){
        yvec<-rnorm(n = nInd, mean = evec)
    }else{
        yvec <- rbinom(n = nInd, size = 1, prob = plogis(evec))
    }
    gt<- colSums(ag)
    return(list(gt = gt, y = yvec, as = as, ag = ag, x = x))
}

#' xs are true design as defined in simulation function
asamap_true <- function(xs, ys, qvec, mafs, # mandatory arguments
                        covs, start, # can be missing
                        model='add', qu=TRUE) {

  model_list <- get_models(mod = model)
  test_list <- get_tests(mod = model)
  models <- list()
  tests <- list()

  if(missing(covs)){
    covs<-matrix(1,nrow=length(ys),ncol=1)
  }

  # Design matrix, one row per observation
  design <- cbind(xs, covs)

  if(missing(start)){
    start <- runif(ncol(design),-1,1)
    if(qu){
      start <- c(start,sd(ys))
    }
  }

  for(mdl in model_list){
    submodel <- get_submodel(f = mdl, s = start, d = design, m = model)
    # Re-using dosage_fit function for fitting model and getting mll
    models[[mdl]] <- dosage_fit(phenos = ys, design = submodel$design, qu = qu)
  }

  for (tst in test_list) {
    tst_mdls <- unlist(strsplit(tst, "v", fixed=TRUE))
    full <- models[[tst_mdls[[1]]]]
    null <- models[[tst_mdls[[2]]]]
    dfs <- length(full$par) - length(null$par)
    m2lq <- 2*(null$value-full$value)
    pval <- pchisq(q = m2lq, df = dfs,lower.tail = FALSE)
    tests[[tst]] <- c(df = dfs, m2lq = m2lq, pval = pval)
  }

  return(list(models=models, tests=tests))
}

sim_fit_all <- function(N=2500, effs, qdist=c(0,0.25,0.5,0.75,1), mafs=c(0.2,0.2), gen_model="add", qu=TRUE, cov){
  Q <- rep(qdist,each=N/length(qdist))

  if(missing(cov)){
    cov <- matrix(1, nrow=N, ncol=1)
    cvec <- 0
  }

  # Simulate data
  asa_sim <- asaSim(qvec=Q, pvec=mafs, bvec=effs, model=gen_model, cvec=cvec, cov=cov, quant=qu)
  # Fit ancestry specific estimates using EM
  asa <- asamap(gs = asa_sim$gt, ys = asa_sim$y, qvec = Q, mafs = mafs, cov=cov, model=gen_model, qu=qu, se=TRUE)
  # Fit ancestry specific estimates using dosage
  asa_dose <- asamap_dosage(gs = asa_sim$gt, ys = asa_sim$y, qvec = Q, mafs = mafs, covs = cov, model = gen_model, qu = qu)
  # Fit using true ancestry
  asa_true <- asamap_true(xs = asa_sim$x, ys = asa_sim$y, qvec = Q, mafs = mafs, covs = cov, model = gen_model, qu = qu)

  return(list(asa=asa, dose=asa_dose, true=asa_true))
}

