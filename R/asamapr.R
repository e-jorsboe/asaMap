#' Make augmented design matrix for additive model from genotypes and covariates
#'
#' @param gs genotypes (vector of length n)
#' @param covs matrix of covariates including intecept vector (matrix, n x p)
#' @return augmented design matrix for full model (matrix 4n x p+3)
make_design_add <- function(gs, covs) {

    add_x1Cg <- matrix(c(0,0,0,0,1,1,0,0,2,1,1,0), nr=4)
    add_x2Cg <- matrix(c(0,0,0,0,0,0,1,1,0,1,1,2), nr=4)
    ## design matrix for alternative allele for pop 1 
    add_x3Cg <- matrix(c(2,1,1,0,1,0,1,0,0,0,0,0), nr=4)
    
    n <- length(gs)
    design <- cbind(x1=c(add_x1Cg[,gs+1]),
                    x2=c(add_x2Cg[,gs+1]),
                    x3=c(add_x3Cg[,gs+1]),
                    covs[rep(1:n,each=4),])
    return(design)
}




#' Make augmented design matrix from recessive model from genotypes and covariates
#'
#' @param gs genotypes (vector of length n)
#' @param covs matrix of covariates including intecept vector (matrix, n x p)
#' @return augmented design matrix for full model (natrix 4n x p+5)
make_design_rec <- function(gs, covs, model) {
    
    rec_x1Cg <- matrix(c(0,0,0,0,0,0,0,0,1,0,0,0), nr=4)
    rec_x2Cg <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,1), nr=4)
    rec_xmCg <- matrix(c(0,0,0,0,0,0,0,0,0,1,1,0), nr=4)
    ## design matrix for being homozygous for alternative allele for pop 1  
    rec_x3Cg <- matrix(c(1,0,0,0,0,0,0,0,0,0,0,0), nr=4)
    ## design matrix for being homozygous for alternative allele for pop 2
    rec_x4Cg <- matrix(c(0,0,0,1,0,0,0,0,0,0,0,0), nr=4)
    
    n <- length(gs)
    design <- cbind(x1=c(rec_x1Cg[,gs+1]),
                    x2=c(rec_x2Cg[,gs+1]),
                    xm=c(rec_xmCg[,gs+1]),
                    x3=c(rec_x3Cg[,gs+1]),
                    x4=c(rec_x4Cg[,gs+1]),
                    covs[rep(1:n,each=4),])
    return(design)
}


#' Make augmented design matrix from genetic model, genotypes and covariates
#'
#' @param gs genotypes (vector of length n)
#' @param covs matrix of covariates including intecept vector (matrix n x p)
#' @param model genetic model eiter "add" or "rec"
#' @return augmented design matrix for full model (4n x p+2 for additive and  4n x p+3 for recessive)
make_design <- function(gs, covs, model) {
  if (model == "add") {
    design <- make_design_add(gs, covs)
  } else if (model=="rec") {
    design <- make_design_rec(gs, covs)
  } else {
    stop("Model must be either 'add' or 'rec', please.")
  }
  return(design)
}


#' Conditional distribution of ancestral state "s" given genotype "g", admixture proportion from pop 1 "q", and allele frequencies "f".
#'
#' @param g genotype (length 1 in {0,1,2})
#' @param q admixture prop from pop 1 (length 1 in [0,1])
#' @param f allele frequencies in ancestral population (length 2 in [0,1]^2)
#' @return prob dist on ancestral state space {1,1},{1,2}.{2,1},{2,2} (length 4 and sums to 1)
prop_states_given_geno <- function(g,q,f){
  ## remember that this is the probability of the state given the genotype, NOT the alleleic genotype, why formulas differ
  t<-c(as.numeric(g>0),as.numeric(g>1))
  un_11 <- q^2*f[1]^(t[1]+t[2])*(1-f[1])^(2-t[1]-t[2])
  un_12 <- q*(1-q)*f[1]^t[1]*(1-f[1])^(1-t[1])*f[2]^t[2]*(1-f[2])^(1-t[2])
  un_21 <- (1-q)*q*f[2]^t[1]*(1-f[2])^(1-t[1])*f[1]^t[2]*(1-f[1])^(1-t[2])
  un_22 <- (1-q)^2*f[2]^(t[1]+t[2])*(1-f[2])^(2-t[1]-t[2])
  c(un_11,un_12,un_21,un_22)/sum(c(un_11,un_12,un_21,un_22))
}

#' Making prior distribution across states
#'
#' @param geno genotype vector (length n)
#' @param admix admixture proportion from pop 1 (vector length n)
#' @param freq allele frequencies in pop 1 and 2 (vector length 2)
#' @return Probability of ancestral state conditional on genotype (vector length 4n with one entry per possible ancestral configuration \{1,1\}, \{1,2\}, \{2,1\}, \{2,2\} for each individual)
make_mixture_props <- function(geno, admix, freq) {
  prior <- c(apply(
    X=cbind(geno,admix),
    MARGIN=1,
    FUN=function(x){prop_states_given_geno(g=x[1],q=x[2],f=freq)}
  ))
  return(prior)
}

#' Making prior distribution across states
#'
#' @param geno genotype vector (length n)
#' @param ancprobs ancestry probabilities (matrix n x 3 OR matrix n x 4)
#' @return Probability of ancestral state conditional on genotype (vector length 4n with one entry per possible ancestral configuration \{1,1\}, \{1,2\}, \{2,1\}, \{2,2\} for each individual)
make_mixture_propsV2 <- function(geno, ancprobs) {
  if(ncol(ancprobs)==4){
    prior<-as.vector(t(ancprobs))
  } else if(ncol(ancprobs)==3){
    ancprobs2<-cbind(ancprobs[,1],ancprobs[,2]/2,ancprobs[,2]/2,ancprobs[,3])
    prior<-as.vector(t(ancprobs2))
  } else{
    q()
    print("Ancestry probabilities (ancprobs) do not have either 3 or 4 columns!")
  }
  return(prior)
}


#' Calculate joint probabilities of phenotypes and ancestral states for quantitative traits
#'
#' @param params Parameters, first genetic effects, then intercept and additional covariate effects and last sd (length 2/3 + p + 1)
#' @param phenos Individual phenotypes (augmented to length 4n, i.e. each repeated 4 times)
#' @param design Augmented design matrix ( matrix 4n x length(params)-1 )
#' @param prior Probability of ancestral state conditional on genotype (vector length 4n) - see ?make_mixture_props
#' @return Joint probability of phenotype and state given parameters (matrix 4 x n)
make_joint_lm <- function(params, phenos, design, prior) {
  probs<-dnorm(x = phenos,
               mean = design%*%(params[-length(params)]),
               sd = params[length(params)])

  ## emil
  print("probs:")
  print(head(probs))
  
  matrix(probs * prior,nr=4)
}

#' Calculate log of joint probabilities of phenotypes and ancestral states for quantitative traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with log of joint probability of phenotype and state given parameters
make_joint_lm_log <- function(params, phenos, design, prior) {
  probs <- dnorm(x = phenos,
                 mean = design%*%(params[-length(params)]),
                 sd = params[length(params)], log=T)

  ## emil

  print('design')
  print(design[5:8,])
  print(params)
  
  print("mean:")
  print(head(design%*%(params[-length(params)]),16))
  print(length(design%*%(params[-length(params)])))
  
  print("sd:")
  print(head(params[length(params)]))

  print("dnorm:")
  print(head(exp(probs)))

  
  matrix(probs + log(prior), nr=4)
}


#' Calculate joint probabilities of phenotypes and ancestral states for dichotomuous traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with joint probability of phenotype and state given parameters
make_joint_bin <- function(params, phenos, design, prior) {
  probs <- dbinom(x=phenos,
                  size=1,
                  prob=plogis(design%*%params))
  matrix(probs * prior,nr=4)
}


#' EM update step for quantitative traits
#'
#' @return list of two: pars_new are new parameter estimates (effects and sd) and mll_new is the updated minus log likelihood
update_em_lm <- function(pars, phen, dsgn, prior) {
  ## joint pheno and state posterior prob
  joint_pheno_state <- make_joint_lm(params = pars, phenos = phen, design = dsgn, prior = prior)
  ## Marginal probs of y (state summed out)
  prob_pheno <- colSums(joint_pheno_state)
  ## minus log likelihood (for old parameters actually)
  mll_new <- -sum(log(prob_pheno))
  ## degrees of freedom for calculatind new sd estimate
  df <- length(prob_pheno) - ncol(dsgn)
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  ## New parameter guess
  fit <- fit_lm(design = dsgn, pheno = phen, weights = wgts)
  #fit <- lm(phen~.-1,weights = wgts,data = data.frame(phen,dsgn))
  sigm <- sqrt(sum((fit$residuals^2) * wgts) / df)
  pars_new <- c(fit$coefficients, sigm)
  return(list(pars_new = pars_new, mll_new = mll_new))
}

#' colSums in log-space
colSumsLog <- function(x){
  k <- apply(x, 2, max)
  log(rowSums(exp(t(x) - k))) + k
}

#' EM update step for quantitative traits, calculations done in log space
#'
#' @param pars Model parameters (length 2/3 + p +1 )
#' @param phen Individual phenotypes augmented to length 4 x N
#' @param dsgn Augmented design matrix (matrix 4n x length(params)-1)
#' @param prior Probability of ancestral state conditional on genotype (vector length 4n) - see ?make_mixture_props
#' @return list of two: pars_new are new parameter estimates (effects and sd) and mll_new is the updated minus log likelihood
update_em_lm_log<-function(pars, phen, dsgn, prior) {
  ## joint pheno and state posterior prob
  joint_pheno_state_log <- make_joint_lm_log(params = pars, phenos = phen,
                                             design = dsgn, prior = prior)


    ## emil
    print("prior:")
    print(head(prior))

    print("joint_pheno_state:")
    print(fed(make_joint_lm(params = pars, phenos = phen,design = dsgn, prior = prior)))
    
    print("joint_pheno_state_log:")
    print(fed(joint_pheno_state_log))


    

  prob_pheno_log <- colSumsLog(joint_pheno_state_log)
  mll_new<- -sum(prob_pheno_log)
  df <- length(phen)/4-ncol(dsgn)
  wgts_log<- c(t(t(joint_pheno_state_log) - prob_pheno_log))
  wgts <- exp(wgts_log)
  ## New parameter guess
  fit <- fit_lm(design = dsgn, pheno = phen, weights = wgts)
  #fit <- lm(phen~.-1,weights = wgts,data = data.frame(phen,dsgn))
  sigm <- sqrt(sum((fit$residuals^2) * wgts) / df)
  pars_new <- c(fit$coefficients, sigm)

  fit <- fit_lm(design = dsgn, pheno = phen, weights = wgts)
  sigm <- sqrt(sum((fit$residuals^2) * wgts) / df)
  pars_new <- c(fit$coefficients, sigm)

  return(list(pars_new = pars_new, mll_new = mll_new))
}



#' EM update step for dichotomuous traits
#'
#' @return list of two: pars_new are new parameter estimates (effects) and mll_new is the updated minus log likelihood
update_em_bin <- function(pars, phen, dsgn, prior) {
  ## joint pheno and state posterior prob
  joint_pheno_state <- make_joint_bin(params = pars, phenos = phen, design = dsgn, prior = prior)
  ## Marginal probs of y (state summed out)
  prob_pheno <- colSums(joint_pheno_state)
  ## minus log likelihood (for old parameters actually)
  mll_new <- -sum(log(prob_pheno))
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  ## New parameter guess
  fit <- fit_bin(design = dsgn, pheno = phen, weights = wgts)
  pars_new <- fit$coefficients
  return(list(pars_new = pars_new, mll_new = mll_new))
}

#' EM update step wrapper for choosing quantitative or dichotomuous traits
update_em <- function(pars, phen, dsgn, prior, quant) {
  if (quant == TRUE) {
    #update_em_lm(pars, phen, dsgn, prior)
    update_em_lm_log(pars, phen, dsgn, prior)
  } else {
    update_em_bin(pars, phen, dsgn, prior)
  }
}

#' Minus log likelihood for full model for quantitative traits (only used by optim checks)
#'
#' @param params Numerical vector of coefficients (and sigma for quantitative)
#' @param desi Augmented design matrix (4 entries per observation) (4n x p)
#' @param phen Augmented phenotypes length 4n
#' @param prior Numerical vector 4 entries per obs of conditional prob of state given genotypes
#' @param quant Logical: is trait quantitiative?
#' @return Minus log of likelihood
calc_mll_lm <- function(params, desi, phen, prior){
  joint_phen_state<-make_joint_lm(params = params, phenos = phen, design = desi, prior = prior)
  return(-sum(log(colSums(joint_phen_state))))
}

#' Minus log likelihood for full model for dicho traits (only used by optim checks)
#'
#' @param params Numerical vector of coefficients (and sigma for quantitative)
#' @param desi Augmented design matrix (4 entries per observation) (4n x p)
#' @param phen Augmented phenotypes length 4n
#' @param prior Numerical vector 4 entries per obs of conditional prob of state given genotypes
#' @param quant Logical: is trait quantitiative?
#' @return Minus log of likelihood
calc_mll_bin <- function(params, desi, phen, prior){
  joint_phen_state <- make_joint_bin(params = params, phenos = phen, design = desi, prior = prior)
  return(-sum(log(colSums(joint_phen_state))))
}

#' EM algorithm controller
#'
#' @param initial Start values for optimization (if quantitative length is ncol(desi)+1 else length is ncol(desi))
#' @param maxI Max number iterations of EM algo
#' @param phe Observed phenotypes 4n long
#' @param desi Sesign matrix 4n times 2+nCov
#' @param pri Prior dist over states, 4n long
#' @param qua Is trait quantitative? true/false
#' @param tole Convergence tolerence
#' @return list of par (estimates), value (minus log likelihood), counts (iterations), convergence and about (how did algo terminate)
control_em <- function(initial, maxI, phe, desi, pri, qua, tole){
  pars_old <- initial
  mll_old <- Inf
  for (iter in 1:maxI) {
    
    update <- update_em(pars = pars_old, phen = phe, dsgn = desi, prior = pri, quant=qua)
    pars_new <- update$pars_new
    mll_new <- update$mll_new
    if (mll_new > mll_old + 10e-6) {
      print(mll_new)
      print(mll_old)
      em_out <- list(par=pars_new, value = mll_new, counts = c(iter, iter), convergence = 1, about = "mll_new > mll_old + 10e-6")
      print("EM step in wrong direction.")
      break
    } else if (mll_old - mll_new < tole) {
      em_out <- list(par = pars_new, value = mll_new, counts = c(iter, iter), convergence = 0, about = "mll_old - mll_new < tol")
      break
    } else if ( iter == maxI){
      em_out <- list(par = pars_new, value = mll_new, counts = c(iter, iter), convergence = 1, about="iter == maxI")
    } else {
      pars_old <- pars_new
      mll_old <- mll_new
    }
  }
  return(em_out)
}


#' Get start parmeters and design matrix for chosen additive (sub)model
#'
#' @param f chosen (sub)model
#' @param s start parameters
#' @param d augmented design matrix
get_submodel_add <- function(f,s,d){
    ## full model with also the alternative allele for pop 1
    if(f=="M0"){
        start <- s
        design <- d        
    }else if(f=="M1"){
        start <- s[-3]
        design <- d[,-3]
    }else if(f=="M2"){
        start <- s[c(-1,-3)]
        design <- d[,c(-1,-3)]
    }else if(f=="M3"){
        start <- s[c(-2,-3)]
        design <- d[,c(-2,-3)]
    }else if(f=="M4"){
        start <- s[c(-1,-3)]
        design <- cbind(d[,1]+d[,2],d[,-c(1,2,3)])
    }else if(f=="M5"){
        start <- s[-c(1,2,3)]
        design <- d[,-c(1,2,3)]
    }else{stop("(Sub)model must be in M1-M5")}
    if(!is.matrix(design)){design <- matrix(design)}
    return(list(start=start,design=design))
}

#' Get start parmeters and design matrix for chosen recessiven (sub)model
#'
#' @param f chosen (sub)model
#' @param s start parameters
#' @param d augmented design matrix
get_submodel_rec <- function(f,s,d){
    ## full model with also the alternative allele for pop 1
    if(f=="R0"){
        start <- s
        design <- d
    }else if(f=="R1"){
        start <- s[-c(4,5)]
        design <- d[,-c(4,5)]
    }else if(f=="R2"){
        start <- s[-c(3,4,5)]
        design <- cbind(d[,1],d[,2]+d[,3],d[,-c(1,2,3,4,5)])
    }else if(f=="R3"){
        start <- s[-c(3,4,5)]
        design <- cbind(d[,1]+d[3],d[,2],d[,-c(1,2,3,4,5)])
    }else if(f=="R4"){
        start <- s[-c(2,3,4,5)]
        design <- d[,-c(2,3,4,5)]
    }else if(f=="R5"){
        start <- s[-c(1,3,4,5)]
        design <- d[,-c(1,3,4,5)]
    }else if(f=="R6"){
        start <- s[-c(2,3,4,5)]
        design <- cbind(d[,1]+d[2]+d[,3],d[,-c(1,2,3,4,5)])
    }else if(f=="R7"){
        start <- s[-c(1,2,3,4,5)]
        design <- d[,-c(1,2,3,4,5)]
    }else{stop("(Sub)model must be in R1-R7")}
    if(!is.matrix(design)){design <- matrix(design)}
    return(list(start=start,design=design))
}

#' Get start parmeters and design matrix for chosen (sub)model
#' #'
#' @param f chosen (sub)model in M1-M5 or R1-R7
#' @param s start parameters
#' @param d augmented design matrix
#' @param m Model - either "add" or "rec"
#' @return List of two, start is the start parameters for the sub model and design is the design matrix for the submodel.
get_submodel <- function(f, s, d, m) {
  if (m == "add") {
    get_submodel_add(f, s, d)
  } else if (m == "rec") {
    get_submodel_rec(f, s, d)
  } else {
    stop("Model must be either 'add' or 'rec', please.")
  }
}

#' Get character vector with names of all submodels for the additive or for the recessive model
#'
#' @param mod Model ("add" or "rec")
get_models <- function(mod) {
  if (mod == "add") {
    return(c("M0","M1","M2","M3","M4","M5"))
  } else if (mod == "rec") {
    return(c("R0","R1","R2","R3","R4","R5","R6","R7"))
  } else {
    stop("Model must be either 'add' or 'rec', please.")
  }
}

#' Potentially to be used to check it the tests the user require are well defined
get_tests <- function(mod) {
  if (mod == "add") {
    small <- c("M1","M5","M2","M3","M4",rep("M5",3))
    large <- c("M0",rep("M1",4),"M2","M3","M4")
  } else if (mod == "rec") {
    small <- c("R1","R2","R3","R4","R5","R6","R7","R4","R6","R7","R5","R6","R7","R7","R7","R7")
    large <- c("R0",rep("R1",6),rep("R2",3),rep("R3",3),"R4","R5","R6")
  } else {
    stop("Model must be either 'add' or 'rec', please.")
  }
  return(paste(large,small,sep="v"))
}


#' Calculating score for getting observed information quantitative traits
#'
calc_score_vector_lm <- function(esti, phen, dsgn, prior) {
  ## Joint distribution of phenotypes and states
  ##print(esti)
  joint_pheno_state <- make_joint_lm(params = esti, phenos = phen, design = dsgn, prior = prior)
  ## Marginal post probs of y (state summed out) - used for normalizing below
  prob_pheno <- colSums(joint_pheno_state)
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  etas <- c(dsgn %*% (esti[-length(esti)]))
  std_resid_wgt <- ((phen-etas) / (esti[length(esti)])^2) * wgts
  matrix(colSums(matrix(c(dsgn * std_resid_wgt), nrow=4)), nrow=length(prob_pheno))
}


#' Calculating score for getting observed information dichotomuous traits
#'
calc_score_vector_bin <- function(esti, phen, dsgn, prior) {
  ## Joint distribution of phenotypes and states
  joint_pheno_state <- make_joint_bin(params = esti, phenos = phen, design = dsgn, prior = prior)
  ## Marginal post probs of y (state summed out) - used for normalizing below
  prob_pheno <- colSums(joint_pheno_state)
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  etas <- c(dsgn %*% esti)
  std_resid_wgt <- (phen-etas) * wgts
  matrix(colSums(matrix(c(dsgn * std_resid_wgt), nrow=4)), nrow=length(prob_pheno))
}

#' Calculating score for getting observed information matrix
calc_score_vector <- function(esti, phen, dsgn, prior, quan) {
  if (quan) {
    return(calc_score_vector_lm(esti, phen, dsgn, prior))
  } else {
    return(calc_score_vector_bin(esti, phen, dsgn, prior))
  }
}

#' Calculating population specific frequencies for a single site using known admixture proportions and genotypes using an EM algorithm
#' 
#' @param start starting guess for population specific frequencies for site
#' @param geno genotypes of site
#' @param a admixture proportions of pop1
#' @return return vector of population specific frequencies (f(pop1), f(pop2))
emFreqAdmix <- function( start, geno, a ){
  q <- rbind(a, 1-a)
  tol <- 1e-4
  aNorm <- t(q) %*% start
  bNorm <- t(q) %*% (1-start)
  updateF <- function(k){
    ag <- geno*q[k,]*start[k]/aNorm
    bg <- (2-geno)*q[k,]*(1-start[k])/bNorm
    fnew <- sum(ag)/(sum(ag) + sum(bg))
    return(fnew)
  }
  f <- sapply(1:2, updateF)
  f[f<tol] <- tol
  f[1-f > 1-tol] <- 1-tol
  return(f)
}


#' Fit ancestry specific models and test different hypothesis (done after removing obs with missing data!)
#'
#' @param gs Genotype vector length N
#' @param ys Phenotype vector of length N
#' @param qvec Admixture proportions from population 1 vector length N
#' @param mafs Frequency of the counted alleles in pop1 and pop2 vector length N
#' @param ancprobs Probabilities of ancestry for each haplotype, matrix (3 cols x N rows) OR (4 cols x N rows) 
#' @param covs Covariate matrix. If missing then intercept is modelled else if intercept is desired then its design matrix must be included in covs matrix.
asamap <- function(gs, ys, # mandatory arguments 
                   qvec, mafs, ancprobs, covs, start, # can be missing
                   model='add', qu=TRUE, maxIter=40, tol=10e-5,
                   se=FALSE, model_list, test_list) {
  if(missing(model_list))
    model_list <- get_models(mod = model) # list of models that should be fit
  if(missing(test_list))
    test_list <- get_tests(mod = model) # List of model comparisons that must be done
  if(missing(covs))
    covs<-matrix(1,nrow=length(ys),ncol=1)
  if(missing(mafs)){
    ## will estimate population specific frequencies from genos and admix proportions
    keep<-!is.na(gs) & !is.na(qvec)
    mafs <- squarem2(c(0.1, 0.1), emFreqAdmix, geno = gs[keep], a = qvec[keep] )$par
  }
  if(!missing(qvec) & missing(ancprobs)){
    prior <- make_mixture_props(geno = gs, admix = qvec, freq = mafs)  
  } else if(missing(qvec) & !missing(ancprobs)){
    prior <- make_mixture_propsV2(geno = gs, ancprobs =  ancprobs)  
  } else {
    print("You should either specify admix proportions (qvec),")
    print("or ancestry probabilities (ancprobs)")
    q()
  }
  
  ## Augmented dataset 4 times size of original
  pheno <- rep(ys,each=4)
  design <- make_design(gs = gs, covs = covs, model = model)
  
  if(missing(start)){
    start <- runif(ncol(design),-1,1)
    if(qu){
      start <- c(start,sd(ys))
    }
  }
  models <- list()
  tests <- list()
  for (mdl in model_list) {

      ## emil
      print("###########################")
      print("###########################")      
      print(mdl)
      print("###########################")
      print("###########################")
      
      submodel <- get_submodel(f = mdl, s = start, d = design, m = model)      
      models[[mdl]] <- control_em(initial = submodel$start, phe = pheno, desi = submodel$design, pri = prior,
                                  qua = qu, tole = tol, maxI = maxIter)

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
  return_list <- list(models = models, tests = tests)
  if(se){ #assumes that first model is the full model
    score_vectors <- calc_score_vector(esti = models[[model_list[1]]]$par, phen = pheno, dsgn = design, prior = prior, quan = qu)
    std_errors <- sqrt(diag(solve(t(score_vectors) %*% score_vectors)))
    return_list$se <- std_errors
  }
  return(return_list)
}

#' Calculate minus log of likelihood used for dosage tests for quantitative traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes length N
#' @param design Expected design matrix, N rows
#' @return Minus log of likelihood
dosage_mll_lm <- function(params, phenos, design) {
  probs<-dnorm(x = phenos,
               mean = design%*%(params[-length(params)]),
               sd = params[length(params)])
  return(-sum(log(probs)))
}

#' Calculate minus log of likelihood used for dosage tests for dichotomuous traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects)
#' @param phenos Individual phenotypes length N
#' @param design Expected design matrix N rows
#' @return Minus log of likelihood
dosage_mll_bin <- function(params, phenos, design) {
  probs <- dbinom(x=phenos,
                  size=1,
                  prob=plogis(design%*%params))
  return(-sum(log(probs)))
}

#' Estimate parameters and calculate minus log of likelihood used for dosage tests
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects (and last sd for quantitative traits)
#' @param phenos Individual phenotypes length N
#' @param design Expected design matrix N rows
#' @return List of Parameters (par) and Minus log of likelihood (value)
dosage_fit <- function(phenos, design, qu){
  if(qu){
    fit <- fit_lm(design = design, pheno = phenos)
    sigm <- sqrt(sum(fit$residuals^2) / fit$df.residual)
    para <- c(fit$coefficients, sigm)
    mll <- dosage_mll_lm(params = para, phenos = phenos, design = design)
  }else{
    fit <- fit_bin(design = design, pheno = phenos)
    para <- fit$coefficients
    mll <- dosage_mll_bin(params = para, phenos = phenos, design = design)
  }
  return(list(par = para, value = mll))
}

asamap_dosage <-function(gs, ys, # mandatory arguments 
                         qvec, mafs, ancprobs, covs, start, # can be missing
                         model='add', qu=TRUE, model_list, test_list) {
  if(missing(model_list))
    model_list <- get_models(mod = model)
  if(missing(test_list))
    test_list <- get_tests(mod = model)
  if(missing(covs)){
    covs<-matrix(1,nrow=length(ys),ncol=1)
  }
  if(missing(mafs)){
    ## will estimate population specific frequencies from genos and admix proportions
    keep<-!is.na(gs) & !is.na(qvec)
    mafs <- squarem2(c(0.1, 0.1), emFreqAdmix, geno = gs[keep], a = qvec[keep] )$par
  }
  if(!missing(qvec) & missing(ancprobs)){
    prior <- make_mixture_props(geno = gs, admix = qvec, freq = mafs)  
  } else if(missing(qvec) & !missing(ancprobs)){
    prior <- make_mixture_propsV2(geno = gs, ancprobs =  ancprobs)  
  } else {
    print("You should either specify admix proportions (qvec),")
    print("or ancestry probabilities (ancprobs)")
    q()
  }
  
  # Augmented design matrix 4 rows per observations
  design <- make_design(gs = gs, covs = covs, model = model)
  
  if(missing(start)){
    start <- runif(ncol(design),-1,1)
    if(qu){
      start <- c(start,sd(ys))
    }
  }
  models <- list()
  tests <- list()
  for(mdl in model_list){
    submodel <- get_submodel(f = mdl, s = start, d = design, m = model)
    expected_design <- apply(X = submodel$design, MARGIN = 2, FUN = function(x){colSums(matrix(x*prior, nrow=4))})
    models[[mdl]] <- dosage_fit(phenos = ys, design = expected_design, qu = qu)
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

