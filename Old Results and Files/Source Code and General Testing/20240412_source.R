library(gtools) # for dirichlet distribution functions
library(matrixStats) # for logSumExp function (helps with rounding errors)
library(parallel) # for parallelization of computations in R (e.g., replace "lapply" with "mclapply")
library(tidyverse) # for plotting and data management functionalities
library(reshape2) # for the melt function
library(abind) # for binding together arrays (increasing data storage as necessary)
library(invgamma) # for sampling inverse-gamma random variables
library(coda) # for mcmc convergence tools
library(scales) # for nice ggplot labeling of log-scaled axes
library(gridExtra) # for side-by-side ggplots
library(ggridges) # create ridgeplots for leniency/harshness
library(label.switching)
library(rankrate)

dbtl <- function(rankings,worth,log=FALSE){
  
  ## Inputs:
  ## (1) rankings, an IxJ matrix of object orderings, where the (i,j) entry is judge i's object in rank j. May also be a vector if I=1.
  ##     NA denotes no object assigned a specific rank. An "assignments" attribute may be used to signify if judges have access to
  ##     only specific objects. "assignments" should also be an IxJ matrix of booleans, where the (i,j) entry denotes if judge i has
  ##     access to object j.
  ## (2) worth, a vector of J positive values, denoting the relative "worths" of each object.
  ## (3) log, a boolean denoting if the loglikelihood should be returned.
  
  ## Data Input Checks
  if(is.vector(rankings)){rankings <- matrix(rankings,nrow=1)}
  if(!is.matrix(rankings)){stop("rankings must be a matrix of rankings")}
  if(is.null(attr(rankings,"assignments"))){attr(rankings,"assignments") <- matrix(TRUE,nrow=nrow(rankings),ncol=ncol(rankings))}
  if(!all(dim(rankings)==dim(attr(rankings,"assignments")))){stop("assignments and rankings must have equal dimensions")}
  if(ncol(rankings)!=length(worth)){stop("rankings must have precisely as many columns as the length of worth")}
  
  ## Calculate Convenient Quantities Now
  I <- nrow(rankings); J <- ncol(rankings)
  log_worth <- log(worth)
  R <- apply(rankings,1,function(ranking){sum(!is.na(ranking))})
  judge_sum_worths <- apply(attr(rankings,"assignments"),1,function(assign){sum(worth[assign])})
  
  ## First Term in log decomposition of joint loglikelihood
  sum_log_numerators <- sum(log_worth[rankings],na.rm=T)
  
  ## Second Term in log decomposition of joint loglikelihood
  denominators <- matrix(c(judge_sum_worths,rep(NA,I*(J-1))),nrow=I,ncol=J)
  for(i in 1:I){
    if(R[i]>1){for(r in 1:(R[i]-1)){
      denominators[i,r+1] <- denominators[i,r]-worth[rankings[i,r]]
    }}
  }
  sum_log_denominators <- sum(log(denominators),na.rm=T)
  
  ## Final density calculation
  log_density <- sum_log_numerators-sum_log_denominators
  if(log){
    return(log_density)
  }else{
    return(exp(log_density))
  }
  
}
dbtlb <- function(rankings,ratings,p,theta,M,log=FALSE){
  
  ## Inputs:
  ## (1) rankings, an IxJ matrix of object orderings, where the (i,j) entry is judge i's object in rank j. May also be a vector if I=1.
  ##     NA denotes no object assigned a specific rank. An "assignments" attribute may be used to signify if judges have access to
  ##     only specific objects. "assignments" should also be an IxJ matrix of booleans, where the (i,j) entry denotes if judge i has
  ##     access to object j.
  ## (2) ratings, an IxJ matrix of ratings, where the (i,j) entry is judge i's rating for rank j, between 0 (best) and M (worst).
  ## (3) p, a vector of J values in the unit interval, denoting the relative quality of each object.
  ## (4) theta, a positive consensus scale parameter.
  ## (5) M, the maximum (worst) integer score
  ## (6) log, a boolean denoting if the loglikelihood should be returned.
  
  # Data Input Checks (rankings are checked by the dbtl function called herein)
  if(is.vector(ratings)){ratings <- matrix(ratings,nrow=1)}
  if(!is.matrix(ratings)){stop("ratings must be a matrix")}
  if(length(p)!=ncol(ratings)){stop("length(p) must equal ncol(ratings)")}
  
  worth <- exp(-theta*p)
  log_density <- dbtl(rankings=rankings,worth=worth,log=T)+
    sum(apply(ratings,1,function(x){dbinom(x=x,size=M,prob=p,log=T)}),na.rm=T)
  
  if(log){
    return(log_density)
  }else{
    return(exp(log_density))
  }
}

rbtlb <- function(I,p,theta,M,R=length(p),assignments=NULL){
  
  J <- length(p)
  if(length(R)==1){R <- rep(R,I)}
  if(is.null(assignments)){assignments <- matrix(TRUE,nrow=I,ncol=J)}
  
  # Draw ratings
  ratings <- as.matrix(t(replicate(I,rbinom(n=J,size=M,prob=p))))
  ratings[!assignments] <- NA
  
  # Draw rankings
  rankings <- matrix(NA,nrow=I,ncol=J)
  for(i in 1:I){
    which_i <- which(assignments[i,])
    rankings[i,1:length(which_i)] <- sample(which_i,prob=exp(-theta*p[which_i]))
    if(R[i]<J){rankings[i,(R[i]+1):J] <- NA}
  }
  attr(rankings,"assignments") <- assignments
  
  return(list(ratings=ratings,rankings=rankings,M=M))
}
expit <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}

btlb_map <- function(rankings,ratings,M,K,tol,verbose=FALSE,
                     delta_gamma,beta_mu,beta_nu,mu_0j,delta_theta){
  
  # Initialize
  I <- nrow(ratings); J <- ncol(ratings)
  gamma_curr <- rgamma(1,delta_gamma[1],delta_gamma[2])
  sigma2_mu_curr <- rinvgamma(1,beta_mu[1],beta_mu[2])
  sigma2_nu_curr <- rinvgamma(1,beta_nu[1],beta_nu[2])
  alpha_curr <- c(rdirichlet(1,rep(gamma_curr,K)))
  mu_curr <- matrix(rnorm(n=J*K,mean=mu_0j,sd=sqrt(sigma2_mu_curr)),nrow=J,ncol=K)
  nu_curr <- rnorm(I,mean=0,sd=sqrt(sigma2_nu_curr))
  theta_curr <- rgamma(K,delta_theta[1],delta_theta[2])
  old_objective <- -Inf
  continue <- TRUE
  
  # EM Algorithm
  while(continue){
    
    # E-STEP STARTS HERE
    
    log_density <- matrix(NA,nrow=I,ncol=K)
    for(i in 1:I){
      rankings_i <- rankings[i,,drop=FALSE]
      attr(rankings_i,"assignments") <- attr(rankings,"assignments")[i,,drop=FALSE]
      for(k in 1:K){
        log_density[i,k] <- dbtlb(rankings=rankings_i,ratings=ratings[i,,drop=FALSE],p=expit(mu_curr[,k]+nu_curr[i]),theta=theta_curr[k],M=M,log=T)
      }
    }
    log_num <- t(log(alpha_curr)+t(log_density))
    z_hat <- exp(log_num-apply(log_num,1,logSumExp)); rm(log_num)
    
    # M-STEP STARTS HERE
    
    ## updating alpha
    alpha_curr <- (gamma_curr-1+apply(z_hat,2,sum))/(K*gamma_curr-K+I)
    alpha_curr[alpha_curr<0]<-0.01
    alpha_curr <- alpha_curr/sum(alpha_curr)
    
    ## updating gamma
    gamma_curr <- optimize(f=function(gamma_new){lgamma(gamma_new*K) - K*lgamma(gamma_new) + (delta_gamma[1]-1)*log(gamma_new) + (sum(log(alpha_curr))-delta_gamma[2])*gamma_new},
                           interval=qgamma(c(0.0001,.9999),delta_gamma[1],delta_gamma[2]),maximum=T)$maximum
    
    ## updating sigma2's
    sigma2_mu_curr <- (beta_mu[2]+sum(c(mu_curr-mu_0j)^2)/2)/(beta_mu[1]+(J*K)/2 + 1) ## posterior mode, based on conjugacy
    sigma2_nu_curr <- (beta_nu[2]+sum(nu_curr^2)/2)/(beta_nu[1]+I/2 + 1) ## posterior mode, based on conjugacy
    
    ## updating mu and theta
    mutheta_new <- lapply(1:K,function(k){
      optim(par=c(mu_curr[,k],theta_curr[k]),
            fn=function(param){
              log_i <- unlist(lapply(1:I,function(i){
                rankings_i <- rankings[i,,drop=FALSE]
                attr(rankings_i,"assignments") <- attr(rankings,"assignments")[i,,drop=FALSE]
                dbtlb(rankings=rankings_i,ratings=ratings[i,,drop=FALSE],p=expit(param[1:J]+nu_curr[i]),theta=param[J+1],M=M,log=TRUE)}))
              
              sum(z_hat[,k]*log_i)+sum(dnorm(param[1:J]-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+dgamma(param[J+1],delta_theta[1],delta_theta[2],log=T)
            },
            method="L-BFGS-B",lower=c(rep(-Inf,J),1e-8),upper=c(rep(Inf,J),qgamma(.9999,delta_theta[1],delta_theta[2])),
            control=list(fnscale=-1))$par
    })
    mu_curr <- matrix(unlist(mutheta_new),ncol=K)[1:J,,drop=FALSE]
    theta_curr <- matrix(unlist(mutheta_new),ncol=K)[J+1,]
    
    ## updating nu's
    nu_new <- matrix(unlist(lapply(1:I,function(i){
      rankings_i <- rankings[i,,drop=FALSE]
      attr(rankings_i,"assignments") <- attr(rankings,"assignments")[i,,drop=FALSE]
      optimize(f=function(nu_new){dnorm(nu_new,mean=0,sd=sqrt(sigma2_nu_curr),log=T)+
          sum(z_hat[i,]*unlist(lapply(1:K,function(k){dbtlb(rankings=rankings_i,ratings=ratings[i,,drop=FALSE],p=expit(mu_curr[,k]+nu_new),theta=theta_curr[k],M=M,log=TRUE)})))},
          interval=qnorm(c(.0001,.9999),mean=0,sd=sqrt(sigma2_nu_curr)),maximum=TRUE)
    })),nrow=2)
    nu_curr <- nu_new[1,]
    
    # CALCULATION OF LOGLIKELIHOOD STARTS HERE - CHECK TO SEE IF EM SHOULD CONTINUE
    
    new_objective <- sum(nu_new[2,])+
      sum(log(alpha_curr)*apply(z_hat,2,sum))+
      dgamma(gamma_curr,delta_gamma[1],delta_gamma[2],log=T)+
      dinvgamma(sigma2_mu_curr,beta_mu[1],beta_mu[2],log=T)+
      dinvgamma(sigma2_nu_curr,beta_nu[1],beta_nu[2],log=T)+
      log(ddirichlet(alpha_curr,rep(gamma_curr,K)))+
      sum(dnorm(c(mu_curr-mu_0j),0,sqrt(sigma2_mu_curr),log=T))+
      sum(dgamma(theta_curr,delta_theta[1],delta_theta[2],log=T))
    
    if(verbose){print(new_objective-old_objective)}
    if(0 < (new_objective-old_objective) & (new_objective-old_objective) < tol){
      continue <- FALSE
    }else{
      old_objective <- new_objective
    }
  }
  
  return(list(gamma=gamma_curr,sigma2_mu=sigma2_mu_curr,sigma2_nu=sigma2_nu_curr,
              alpha=alpha_curr,z_hat=z_hat,mu=mu_curr,nu=nu_curr,theta=theta_curr,
              K=K,loglikelihood=new_objective))
}
btlb_map_noratereffects <- function(rankings,ratings,M,K,tol,verbose=FALSE,
                         delta_gamma,beta_mu,mu_0j,delta_theta){
  dbtlb_fast_log <- function(rankings,ratings,log_worth,worth,p,M,J,R){
    
    denominators <- c(sum(worth),rep(NA,J-1))
    if(R>1){for(r in 1:(R-1)){denominators[r+1] <- denominators[r]-worth[rankings[1,r]]}}
    
    return(sum(log_worth[rankings],na.rm=T)-sum(log(denominators),na.rm=T)+
             sum(apply(ratings,1,function(x){dbinom(x=x,size=M,prob=p,log=T)}),na.rm=T))
  }
  # Initialize
  continue <- TRUE
  I <- nrow(ratings); J <- ncol(ratings)
  R <- apply(rankings,1,function(pi){J-sum(is.na(pi))})
  gamma_curr <- rgamma(1,delta_gamma[1],delta_gamma[2])
  sigma2_mu_curr <- rinvgamma(1,beta_mu[1],beta_mu[2])
  alpha_curr <- rep(1/K,K)
  z_hat <- rdirichlet(I,alpha_curr)
  mu_curr <- matrix(0,nrow=J,ncol=K); theta_curr <- rep(delta_theta[1]/delta_theta[2],K)
  for(k in 1:K){
    mutheta_new <- optim(par=c(mu_curr[,k],theta_curr[k]),
                         fn=function(param){
                           mu_tmp <- param[1:J]
                           theta_tmp <- param[J+1]
                           p_tmp <- expit(mu_tmp)
                           log_worth_tmp <- -theta_tmp*p_tmp
                           worth_tmp <- exp(log_worth_tmp)
                           
                           log_dens <- unlist(lapply(1:I,function(i){dbtlb_fast_log(rankings[i,,drop=FALSE],ratings[i,,drop=FALSE],log_worth_tmp,worth_tmp,p_tmp,M,J,R[i])}))
                           sum(z_hat[,k]*log_dens)+sum(dnorm(mu_tmp,mu_0j,sd=sqrt(sigma2_mu_curr),log=T))+dgamma(theta_tmp,delta_theta[1],delta_theta[2],log=T)
                         },
                         method="L-BFGS-B",lower=c(rep(-10,J),1e-8),upper=c(rep(10,J),100),
                         control=list(fnscale=-1))$par
    mu_curr[,k] <- mutheta_new[1:J]
    theta_curr[k] <- mutheta_new[J+1]
  }
  
  log_density <- matrix(NA,nrow=I,ncol=K)
  for(k in 1:K){
    p_tmp <- expit(mu_curr[,k])
    log_worth_tmp <- -theta_curr[k]*p_tmp
    worth_tmp <- exp(log_worth_tmp)
    for(i in 1:I){
      log_density[i,k] <- dbtlb_fast_log(rankings[i,,drop=FALSE],ratings[i,,drop=FALSE],log_worth_tmp,worth_tmp,p_tmp,M,J,R[i])
    }}
  old_objective <- sum(apply(z_hat,2,sum)*log(alpha_curr))+
    sum(log_density*z_hat)+
    dgamma(gamma_curr,delta_gamma[1],delta_gamma[2],log=T)+
    dinvgamma(sigma2_mu_curr,beta_mu[1],beta_mu[2],log=T)+
    log(ddirichlet(alpha_curr,rep(gamma_curr,K)))+
    sum(dnorm(c(mu_curr-mu_0j),0,sqrt(sigma2_mu_curr),log=T))+
    sum(dgamma(theta_curr,delta_theta[1],delta_theta[2],log=T))
  
  # EM Algorithm
  while(continue){
    
    # E-STEP STARTS HERE
    log_num <- t(log(alpha_curr)+t(log_density))
    z_hat <- exp(log_num-apply(log_num,1,logSumExp)); rm(log_num)
    
    # M-STEP STARTS HERE
    
    ## updating alpha
    alpha_curr <- (gamma_curr-1+apply(z_hat,2,sum))/(K*gamma_curr-K+I)
    alpha_curr[alpha_curr<0]<-0.01
    alpha_curr <- alpha_curr/sum(alpha_curr)
    
    ## updating gamma
    gamma_curr <- optimize(f=function(gamma_new){lgamma(gamma_new*K) - K*lgamma(gamma_new) + (delta_gamma[1]-1)*log(gamma_new) + (sum(log(alpha_curr))-delta_gamma[2])*gamma_new},
                           interval=qgamma(c(0.0001,.9999),delta_gamma[1],delta_gamma[2]),maximum=T)$maximum
    
    ## updating mu and theta
    for(k in 1:K){
      mutheta_new <- optim(par=c(mu_curr[,k],theta_curr[k]),
                           fn=function(param){
                             mu_tmp <- param[1:J]
                             theta_tmp <- param[J+1]
                             p_tmp <- expit(mu_tmp)
                             log_worth_tmp <- -theta_tmp*p_tmp
                             worth_tmp <- exp(log_worth_tmp)
                             
                             log_dens <- unlist(lapply(1:I,function(i){dbtlb_fast_log(rankings[i,,drop=FALSE],ratings[i,,drop=FALSE],log_worth_tmp,worth_tmp,p_tmp,M,J,R[i])}))
                             sum(z_hat[,k]*log_dens)+sum(dnorm(mu_tmp,mu_0j,sd=sqrt(sigma2_mu_curr),log=T))+dgamma(theta_tmp,delta_theta[1],delta_theta[2],log=T)
                           },
                           method="L-BFGS-B",lower=c(rep(-10,J),1e-8),upper=c(rep(10,J),100),
                           control=list(fnscale=-1))$par
      mu_curr[,k] <- mutheta_new[1:J]
      theta_curr[k] <- mutheta_new[J+1]
    }
    
    ## updating sigma2 - posterior mode, based on conjugacy
    sigma2_mu_curr <- (beta_mu[2]+sum(c(mu_curr-mu_0j)^2)/2)/(beta_mu[1]+(J*K)/2 + 1)
    
    
    # CALCULATION OF LOGLIKELIHOOD STARTS HERE - CHECK TO SEE IF EM SHOULD CONTINUE
    for(k in 1:K){
      p_tmp <- expit(mu_curr[,k])
      log_worth_tmp <- -theta_curr[k]*p_tmp
      worth_tmp <- exp(log_worth_tmp)
      for(i in 1:I){
        log_density[i,k] <- dbtlb_fast_log(rankings[i,,drop=FALSE],ratings[i,,drop=FALSE],log_worth_tmp,worth_tmp,p_tmp,M,J,R[i])
      }}
    new_objective <- sum(apply(z_hat,2,sum)*log(alpha_curr))+
      sum(log_density*z_hat)+
      dgamma(gamma_curr,delta_gamma[1],delta_gamma[2],log=T)+
      dinvgamma(sigma2_mu_curr,beta_mu[1],beta_mu[2],log=T)+
      log(ddirichlet(alpha_curr,rep(gamma_curr,K)))+
      sum(dnorm(c(mu_curr-mu_0j),0,sqrt(sigma2_mu_curr),log=T))+
      sum(dgamma(theta_curr,delta_theta[1],delta_theta[2],log=T))
    
    
    if(verbose){print(new_objective-old_objective)}
    if(0 < (new_objective-old_objective) & (new_objective-old_objective) < tol){
      continue <- FALSE
    }else{
      old_objective <- new_objective
    }
  }
  
  return(list(gamma=gamma_curr,sigma2_mu=sigma2_mu_curr,
              alpha=alpha_curr,z_hat=z_hat,mu=mu_curr,theta=theta_curr,
              K=K,loglikelihood=new_objective))
}
btlb_mfmm <- function(rankings,ratings,M,K_start,init,delta_theta,delta_gamma,lambda,beta_mu,beta_nu,mu_0j,
                      T1,T2,sigma_mu,sigma_nu,sigma_theta,sigma_gamma,rater_effects=TRUE){
  # Internally-Used Functions
  log_dbtlb_nodatachecks <- function(ranking,rating,assignment,p,worth,logworth,J,R,M){
    # rankings density - numerator
    log_numerator <- sum(logworth[ranking],na.rm=T)
    # rankings density - denominator
    summed_worths_assignments <- exp(logSumExp(logworth[assignment]))
    summed_worths_uptocurrent <- c(0,exp(unlist(lapply(1:(J-1),function(j){logSumExp(logworth[ranking[1:j]])}))))
    if(R<J){summed_worths_uptocurrent[(R+1):J] <- NA}
    log_denominator <- sum(log(summed_worths_assignments-summed_worths_uptocurrent),na.rm=T)
    # ratings density
    log_binomial <- sum(dbinom(rating,M,p,log=T),na.rm=T)
    # overall density 
    return(log_numerator - log_denominator + log_binomial)
  }
  
  # Data Checks
  if(!is.matrix(rankings) | !is.matrix(ratings)){stop("rankings and ratings must be matrices")}
  I <- nrow(rankings); J <- ncol(rankings)
  R <- apply(rankings,1,function(ranking){sum(!is.na(ranking))})
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of equal dimension")}
  if(M %% 1 != 0 | M<=0){stop("M must be a positive integer")}
  
  if(is.null(attr(rankings,"assignments"))){
    assignments <- matrix(TRUE,nrow=nrow(rankings),ncol=ncol(rankings))
  }else{assignments <- attr(rankings,"assignments")}
  if(any(!is.logical(assignments))){stop("assignments must be boolean-valued")}
  if(any(dim(rankings) != dim(assignments))){stop("rankings and assignments must have equal dimensions")}
  if(any(!is.na(ratings[assignments==FALSE]))){stop("No ratings may be assignment to an object a judge is not assigned to")}
  if(any(unlist(lapply(1:I,function(i){ifelse(length(setdiff(na.exclude(rankings[i,]),which(assignments[i,]==TRUE)))>0,TRUE,FALSE)})))){
    stop("No ranks may be given to an object a judge is not assigned to")}
  
  if(length(delta_theta)!=2){stop("delta_theta must be a vector of length 2")}
  if(length(delta_gamma)!=2){stop("delta_gamma must be a vector of length 2")}
  if(length(lambda)!=1 | lambda <=0){stop("lambda must be a positive number.")}
  if(length(beta_mu)!=2){stop("beta_mu must be a vector of length 2")}
  if(length(beta_nu)!=2){stop("beta_nu must be a vector of length 2")}
  if(length(mu_0j)!=J){stop("mu_0j must be a vector of length ncol(ratings)")}
  
  # Step 1: Initialize and Data Storage
  print("Initializing")
  K_curr <- K_start
  if(is.null(init)){
    gamma_curr <- rgamma(1,delta_gamma[1],delta_gamma[2])
    sigma2_mu_curr <- rinvgamma(1,beta_mu[1],beta_mu[2])
    sigma2_nu_curr <- rinvgamma(1,beta_nu[1],beta_nu[2])
    alpha_curr <- c(rdirichlet(1,rep(gamma_curr,K_start)))
    mu_curr <- matrix(rnorm(n=J*K_start,mean=mu_0j,sd=sqrt(sigma2_mu_curr)),nrow=J,ncol=K_start)
    nu_curr <- rnorm(I,mean=0,sd=sqrt(sigma2_nu_curr))
    theta_curr <- rgamma(K_start,delta_theta[1],delta_theta[2])
  }else{
    gamma_curr <- init$gamma
    sigma2_mu_curr <- init$sigma2_mu
    sigma2_nu_curr <- init$sigma2_nu
    alpha_curr <- init$alpha
    mu_curr <- init$mu
    nu_curr <- init$nu
    theta_curr <- init$theta
  }
  Z_curr <- rep(NA,I)
  
  Z_iters <- matrix(NA,nrow=T1,ncol=I)
  mu_iters <- array(NA,dim=c(J,K_start,T1))
  nu_iters <- matrix(NA,nrow=T1,ncol=I)
  theta_iters <- matrix(NA,nrow=T1,ncol=K_start)
  gamma_iters <- matrix(NA,nrow=T1,ncol=1)
  alpha_iters <- matrix(NA,nrow=T1,ncol=K_start)
  K_iters <- matrix(NA,nrow=T1,ncol=2,dimnames=list(1:T1,c("Kplus","K")))
  sigma2_iters <- matrix(NA,nrow=T1,ncol=2,dimnames = list(1:T1,c("mu","nu")))
  accept <- matrix(0,nrow=T1,ncol=4)
  colnames(accept) <- c("mu","theta","nu","gamma")
  timing <- matrix(NA,nrow=T1,ncol=4)
  colnames(timing) <- c("2a","2b","2c","2d")
  checkin_iters <- round(seq(0,T1,length=11))[-1]
  
  if(rater_effects==FALSE){nu_curr <- rep(0,I)}
  
  loglik_curr <- matrix(NA,nrow=I,ncol=K_start)
  for(i in 1:I){for(k in 1:K_start){
    p_i <- expit(mu_curr[,k]+nu_curr[i])
    log_worth <- -theta_curr[k]*p_i
    loglik_curr[i,k] <- log_dbtlb_nodatachecks(ranking=rankings[i,],rating=ratings[i,],assignment=assignments[i,],
                                               p=p_i,worth=exp(log_worth),logworth=log_worth,J=J,R=R[i],M=M)
    rm(p_i,log_worth)
  }}
  loglik_iters <- matrix(NA,nrow=T1,ncol=1)
  class_iters <- array(0,dim=c(T1,I,K_start))
  
  print("Starting Gibbs Iterations")
  # Step 2:
  for(t1 in 1:T1){
    
    ## Step 2(a)
    time1 <- Sys.time()
    ### 2(a)(i) Draw new class labels
    # draw new Z_curr
    loglik_curr[is.infinite(loglik_curr)] <- -Inf
    unnorm_log_classprobs <- t(loglik_curr)+log(alpha_curr)
    classprobs <- exp(t(unnorm_log_classprobs)-apply(unnorm_log_classprobs,2,logSumExp))
    if(t1>1){class_iters[t1-1,,1:K_curr] <- classprobs}
    Z_curr <- apply(classprobs,1,function(prob){sample(1:K_curr,1,prob=prob)})
    
    ### 2(a)(ii) Update Nk_curr and K_curr, plus some re-labeling
    Nk_curr <- unlist(lapply(1:K_curr,function(k){sum(Z_curr==k)}))
    which_keep <- which(Nk_curr!=0)
    mu_curr <- mu_curr[,which_keep,drop=FALSE]
    theta_curr <- theta_curr[which_keep]
    loglik_curr <- loglik_curr[,which_keep,drop=FALSE]
    Z_curr <- unlist(lapply(1:I,function(i){which(which_keep == Z_curr[i])}))
    Nk_curr <- Nk_curr[which_keep]
    Kplus_curr <- length(Nk_curr)
    
    
    ## Step 2(b)
    time2 <- Sys.time()
    ### Step 2(b)(i) By class k, update (mu_k,theta_k)
    accept_mu <- accept_theta <- c()
    for(k in 1:Kplus_curr){
      ik <- which(Z_curr==k)
      curr_dens <- sum(loglik_curr[ik,k]) + 
        sum(dnorm(mu_curr[,k]-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
        dgamma(theta_curr[k],delta_theta[1],delta_theta[2],log=T)
      for(t2 in 1:T2){
        ### Step 2(b)(i) Update mu_jk
        for(j in 1:J){
          mu_new <- mu_curr[,k]
          mu_new[j] <- rnorm(1,mean=mu_curr[j,k],sd=sigma_mu)
          loglik_new <- unlist(lapply(ik,function(i){
            p_i <- expit(mu_new+nu_curr[i])
            log_worth <- -theta_curr[k]*p_i
            log_dbtlb_nodatachecks(ranking=rankings[i,],rating=ratings[i,],assignment=assignments[i,],
                                   p=p_i,worth=exp(log_worth),logworth=log_worth,J=J,R=R[i],M=M)
          }))
          new_dens <- sum(loglik_new) + 
            sum(dnorm(mu_new-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
            dgamma(theta_curr[k],delta_theta[1],delta_theta[2],log=T)
          
          diff <- new_dens-curr_dens
          if(!is.nan(diff) & !is.infinite(diff) & log(runif(1)) < diff){
            accept_mu <- c(accept_mu,1)
            mu_curr[j,k] <- mu_new[j]
            curr_dens <- new_dens
            loglik_curr[ik,k] <- loglik_new
          }else{accept_mu <- c(accept_mu,0)}
        }
        ### Step 2(b)(ii) Update each theta_k
        theta_new <- rnorm(1,mean=theta_curr[k],sd=sigma_theta)
        if(theta_new >0){
          loglik_new <- unlist(lapply(ik,function(i){
            p_i <- expit(mu_curr[,k]+nu_curr[i])
            log_worth <- -theta_new*p_i
            log_dbtlb_nodatachecks(ranking=rankings[i,],rating=ratings[i,],assignment=assignments[i,],
                                   p=p_i,worth=exp(log_worth),logworth=log_worth,J=J,R=R[i],M=M)
          }))
          new_dens <- sum(loglik_new) + 
            sum(dnorm(mu_curr[,k]-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
            dgamma(theta_new,delta_theta[1],delta_theta[2],log=T)
          
          diff <- new_dens-curr_dens
          if(!is.nan(diff) & !is.infinite(diff) & log(runif(1)) < diff){
            accept_theta <- c(accept_theta,1)
            theta_curr[k] <- theta_new
            curr_dens <- new_dens
            loglik_curr[ik,k] <- loglik_new
          }else{accept_theta <- c(accept_theta,0)}
        }else{accept_theta <- c(accept_theta,0)}
      }
      loglik_curr[-ik,k] <- unlist(lapply(setdiff(1:I,ik),function(i){
        p_i <- expit(mu_curr[,k]+nu_curr[i])
        log_worth <- -theta_curr[k]*p_i
        log_dbtlb_nodatachecks(ranking=rankings[i,],rating=ratings[i,],assignment=assignments[i,],
                               p=p_i,worth=exp(log_worth),logworth=log_worth,J=J,R=R[i],M=M)
      }))
    }
    accept[t1,1:2] <- c(mean(accept_mu),mean(accept_theta))
    
    ### Step 2(b)(ii) By judge i, update nu_i
    accept_nu <- c()
    if(rater_effects){
      for(i in 1:I){
        Z_i <- Z_curr[i]
        nu_new <- rnorm(1,mean=nu_curr[i],sd=sigma_nu)
        p_i <- expit(mu_curr[,Z_i]+nu_new)
        log_worth <- -theta_curr[Z_i]*p_i
        loglik_new <- log_dbtlb_nodatachecks(ranking=rankings[i,],rating=ratings[i,],assignment=assignments[i,],
                                             p=p_i,worth=exp(log_worth),logworth=log_worth,J=J,R=R[i],M=M)
        
        diff <- loglik_new-loglik_curr[i,Z_i] + dnorm(nu_new,mean=0,sd=sqrt(sigma2_nu_curr),log=T)-dnorm(nu_curr[i],mean=0,sd=sqrt(sigma2_nu_curr),log=T)
        if(!is.nan(diff) & !is.infinite(diff) & log(runif(1)) < diff){
          accept_nu <- c(accept_nu,1)
          nu_curr[i] <- nu_new
          loglik_curr[i,Z_i] <- loglik_new
          loglik_curr[i,-Z_i] <- NA
        }else{accept_nu <- c(accept_nu,0)}
      }
      accept[t1,3] <- mean(accept_nu)
    }
    
    ### Step 2(b)(iii) Update sigma2_mu and sigma2_nu
    sigma2_mu_curr <- rinvgamma(1,beta_mu[1]+J*K_curr/2,beta_mu[2]+sum((mu_curr-mu_0j)^2)/2)
    sigma2_nu_curr <- rinvgamma(1,beta_nu[1]+I/2,beta_nu[2]+sum(nu_curr^2)/2)
    
    
    ## Step 2(c)
    time3 <- Sys.time()
    ### Step 2(c)(i) Update K_curr
    probs <- unlist(lapply(Kplus_curr:(Kplus_curr+100),function(k){
      dpois(k-1,lambda,log=TRUE)+lfactorial(k)-lfactorial(k-Kplus_curr)+lgamma(gamma_curr*k)-lgamma(I+gamma_curr*k)
    }))
    K_curr <- sample(Kplus_curr:(Kplus_curr+100),1,prob = exp(probs-logSumExp(probs)))
    if(K_curr > ncol(alpha_iters)){
      mu_iters <- abind(mu_iters, array(NA, replace(dim(mu_iters), 2, K_curr-ncol(alpha_iters))), along = 2)
      theta_iters <- cbind(theta_iters,matrix(NA,nrow=T1,ncol=K_curr-ncol(alpha_iters)))
      class_iters <- abind(class_iters, array(0, replace(dim(class_iters), 3, K_curr-ncol(alpha_iters))), along = 3)
      alpha_iters <- cbind(alpha_iters,matrix(NA,nrow=T1,ncol=K_curr-ncol(alpha_iters)))
    }
    
    ### Step 2(c)(ii) Update gamma_curr
    gamma_new <- abs(rnorm(1,gamma_curr,sd=sigma_gamma))
    dens_curr <- dgamma(gamma_curr,delta_gamma[1],delta_gamma[2],log=T)+lgamma(gamma_curr*K_curr)-
      lgamma(I+gamma_curr*K_curr)+sum(lgamma(Nk_curr+gamma_curr)-lgamma(gamma_curr))
    dens_new <-  dgamma(gamma_new, delta_gamma[1],delta_gamma[2],log=T)+lgamma(gamma_new*K_curr)-
      lgamma(I+gamma_new*K_curr )+sum(lgamma(Nk_curr+gamma_new )-lgamma(gamma_new))
    if(log(runif(1))<(dens_new-dens_curr)){
      accept[t1,4] <- 1
      gamma_curr <- gamma_new
    }
    
    ## Step 2(d)
    time4 <- Sys.time()
    ### Step 2(d)(i) Create new (mu_k,theta_k) if K > K+
    if(K_curr > Kplus_curr){
      for(k in (Kplus_curr+1):K_curr){
        mu_curr <- cbind(mu_curr,rnorm(J,mean=mu_0j,sd=sqrt(sigma2_mu_curr)))
        theta_curr <- c(theta_curr,rgamma(1,delta_theta[1],delta_theta[2]))
        Nk_curr <- c(Nk_curr,0)
      }
      loglik_curr <- cbind(loglik_curr,matrix(NA,nrow=I,ncol=K_curr-Kplus_curr))
    }
    # fill in loglik_curr matrix
    for(i in 1:I){
      which_k <- which(is.na(loglik_curr[i,]))
      for(k in which_k){
        p_i <- expit(mu_curr[,k]+nu_curr[i])
        log_worth <- -theta_curr[k]*p_i
        loglik_curr[i,k] <- log_dbtlb_nodatachecks(ranking=rankings[i,],rating=ratings[i,],assignment=assignments[i,],
                                                   p=p_i,worth=exp(log_worth),logworth=log_worth,J=J,R=R[i],M=M)
      }
    }
    
    ### Step 2(d)(ii) Update alpha
    alpha_curr <- c(rdirichlet(1,gamma_curr+Nk_curr))
    time5 <- Sys.time()
    
    ## Store Current Values!
    Z_iters[t1,] <- Z_curr
    mu_iters[,1:K_curr,t1] <- mu_curr
    nu_iters[t1,] <- nu_curr
    theta_iters[t1,1:K_curr] <- theta_curr
    gamma_iters[t1,] <- gamma_curr
    alpha_iters[t1,1:K_curr] <- alpha_curr
    K_iters[t1,] <- c(Kplus_curr,K_curr)
    sigma2_iters[t1,] <- c(sigma2_mu_curr,sigma2_nu_curr)
    loglik_iters[t1,] <- sum(unlist(lapply(1:I,function(i){log(alpha_curr[Z_curr[i]])+loglik_curr[i,Z_curr[i]]})))+
      log(ddirichlet(alpha_curr,rep(gamma_curr,K_curr)))+
      dgamma(gamma_curr,delta_gamma[1],delta_gamma[2],log=T)+
      dpois(K_curr-1,lambda,log=T)+
      sum(dgamma(theta_curr,delta_theta[1],delta_theta[2],log=T))+
      sum(dnorm(mu_curr-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
      sum(dnorm(nu_curr,mean=0,sd=sqrt(sigma2_nu_curr),log=T))+
      dinvgamma(sigma2_mu_curr,beta_mu[1],beta_mu[2],log=T)+
      dinvgamma(sigma2_nu_curr,beta_nu[1],beta_nu[2],log=T)
    timing[t1,] <- difftime(c(time2,time3,time4,time5),c(time1,time2,time3,time4),units="secs")
    if(t1 %in% checkin_iters){
      text <- paste0("Status: ",which(t1 == checkin_iters)*10,"% complete. Current class sizes: ",paste0(Nk_curr,collapse=", "))
      print(text)
    }
  }
  # updating final class probabilities
  loglik_curr[is.infinite(loglik_curr)] <- -Inf
  unnorm_log_classprobs <- t(loglik_curr)+log(alpha_curr)
  classprobs <- exp(t(unnorm_log_classprobs)-apply(unnorm_log_classprobs,2,logSumExp))
  class_iters[t1,,1:K_curr] <- classprobs
  
  print("Done!")
  
  # Return Results!
  list(gamma=gamma_iters,sigma2=sigma2_iters,alpha=alpha_iters,
       mu=mu_iters,nu=nu_iters,theta=theta_iters,Z=Z_iters,K=K_iters,loglik=loglik_iters,class=class_iters,
       accept=accept,timing=timing)
}
btlb_mfmm_noratereffects <- function(rankings,ratings,M,K_start,init,delta_theta,delta_gamma,lambda,beta_mu,mu_0j,
                                     T1,T2,sigma_mu,sigma_theta,sigma_gamma){
  # Internally-Used Functions
  dbtlb_fast_log <- function(rankings,ratings,log_worth,worth,p,M,J,R){
    
    denominators <- c(sum(worth),rep(NA,J-1))
    if(R>1){for(r in 1:(R-1)){denominators[r+1] <- denominators[r]-worth[rankings[1,r]]}}
    
    return(sum(log_worth[rankings],na.rm=T)-sum(log(denominators),na.rm=T)+
             sum(apply(ratings,1,function(x){dbinom(x=x,size=M,prob=p,log=T)}),na.rm=T))
  }
  
  # Data Checks
  if(!is.matrix(rankings) | !is.matrix(ratings)){stop("rankings and ratings must be matrices")}
  I <- nrow(rankings); J <- ncol(rankings)
  R <- apply(rankings,1,function(ranking){sum(!is.na(ranking))})
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of equal dimension")}
  if(M %% 1 != 0 | M<=0){stop("M must be a positive integer")}
  
  if(is.null(attr(rankings,"assignments"))){
    assignments <- matrix(TRUE,nrow=nrow(rankings),ncol=ncol(rankings))
  }else{assignments <- attr(rankings,"assignments")}
  if(any(!is.logical(assignments))){stop("assignments must be boolean-valued")}
  if(any(dim(rankings) != dim(assignments))){stop("rankings and assignments must have equal dimensions")}
  if(any(!is.na(ratings[assignments==FALSE]))){stop("No ratings may be assignment to an object a judge is not assigned to")}
  if(any(unlist(lapply(1:I,function(i){ifelse(length(setdiff(na.exclude(rankings[i,]),which(assignments[i,]==TRUE)))>0,TRUE,FALSE)})))){
    stop("No ranks may be given to an object a judge is not assigned to")}
  
  if(length(delta_theta)!=2){stop("delta_theta must be a vector of length 2")}
  if(length(delta_gamma)!=2){stop("delta_gamma must be a vector of length 2")}
  if(length(lambda)!=1 | lambda <=0){stop("lambda must be a positive number.")}
  if(length(beta_mu)!=2){stop("beta_mu must be a vector of length 2")}
  if(length(mu_0j)!=J){stop("mu_0j must be a vector of length ncol(ratings)")}
  
  # Step 1: Initialize and Data Storage
  print("Initializing")
  K_curr <- K_start
  if(is.null(init)){
    gamma_curr <- rgamma(1,delta_gamma[1],delta_gamma[2])
    sigma2_mu_curr <- rinvgamma(1,beta_mu[1],beta_mu[2])
    alpha_curr <- c(rdirichlet(1,rep(gamma_curr,K_start)))
    mu_curr <- matrix(rnorm(n=J*K_start,mean=mu_0j,sd=sqrt(sigma2_mu_curr)),nrow=J,ncol=K_start)
    theta_curr <- rgamma(K_start,delta_theta[1],delta_theta[2])
  }else{
    gamma_curr <- init$gamma
    sigma2_mu_curr <- init$sigma2_mu
    alpha_curr <- init$alpha
    mu_curr <- init$mu
    theta_curr <- init$theta
  }
  Z_curr <- rep(NA,I)
  
  Z_iters <- matrix(NA,nrow=T1,ncol=I)
  mu_iters <- array(NA,dim=c(J,K_start,T1))
  theta_iters <- matrix(NA,nrow=T1,ncol=K_start)
  gamma_iters <- matrix(NA,nrow=T1,ncol=1)
  alpha_iters <- matrix(NA,nrow=T1,ncol=K_start)
  K_iters <- matrix(NA,nrow=T1,ncol=2,dimnames=list(1:T1,c("Kplus","K")))
  sigma2_iters <- matrix(NA,nrow=T1,ncol=1,dimnames = list(1:T1,c("mu")))
  accept <- matrix(0,nrow=T1,ncol=3)
  colnames(accept) <- c("mu","theta","gamma")
  timing <- matrix(NA,nrow=T1,ncol=4)
  colnames(timing) <- c("2a","2b","2c","2d")
  checkin_iters <- round(seq(0,T1,length=11))[-1]
  
  loglik_curr <- matrix(NA,nrow=I,ncol=K_start)
  for(k in 1:K_start){
    p_k <- expit(mu_curr[,k])
    log_worth_k <- -theta_curr[k]*p_k
    worth_k <- exp(log_worth_k)
    for(i in 1:I){
      loglik_curr[i,k] <- dbtlb_fast_log(rankings=rankings[i,,drop=FALSE],ratings=ratings[i,,drop=FALSE],
                                         log_worth=log_worth_k,worth=worth_k,p=p_k,M=M,J=J,R=R[i])
    }}
  loglik_iters <- matrix(NA,nrow=T1,ncol=1)
  class_iters <- array(0,dim=c(T1,I,K_start))
  
  print("Starting Gibbs Iterations")
  # Step 2:
  for(t1 in 1:T1){
    
    ## Step 2(a)
    time1 <- Sys.time()
    ### 2(a)(i) Draw new class labels
    # draw new Z_curr
    loglik_curr[is.infinite(loglik_curr)] <- -Inf
    unnorm_log_classprobs <- t(loglik_curr)+log(alpha_curr)
    classprobs <- exp(t(unnorm_log_classprobs)-apply(unnorm_log_classprobs,2,logSumExp))
    if(t1>1){class_iters[t1-1,,1:K_curr] <- classprobs}
    Z_curr <- apply(classprobs,1,function(prob){sample(1:K_curr,1,prob=prob)})
    
    ### 2(a)(ii) Update Nk_curr and K_curr, plus some re-labeling
    Nk_curr <- unlist(lapply(1:K_curr,function(k){sum(Z_curr==k)}))
    which_keep <- which(Nk_curr!=0)
    mu_curr <- mu_curr[,which_keep,drop=FALSE]
    theta_curr <- theta_curr[which_keep]
    loglik_curr <- loglik_curr[,which_keep,drop=FALSE]
    Z_curr <- unlist(lapply(1:I,function(i){which(which_keep == Z_curr[i])}))
    Nk_curr <- Nk_curr[which_keep]
    Kplus_curr <- length(Nk_curr)
    
    
    ## Step 2(b)
    time2 <- Sys.time()
    ### Step 2(b)(i) By class k, update (mu_k,theta_k)
    accept_mu <- accept_theta <- c()
    for(k in 1:Kplus_curr){
      ik <- which(Z_curr==k)
      curr_dens <- sum(loglik_curr[ik,k]) + 
        sum(dnorm(mu_curr[,k]-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
        dgamma(theta_curr[k],delta_theta[1],delta_theta[2],log=T)
      p_k <- expit(mu_curr[,k])
      log_worth_k <- -theta_curr[k]*p_k
      worth_k <- exp(log_worth_k)
      
      for(t2 in 1:T2){
        ### Step 2(b)(i) Update mu_jk
        for(j in 1:J){
          mu_new <- mu_curr[,k];mu_new[j] <- rnorm(1,mean=mu_curr[j,k],sd=sigma_mu)
          p_k_new <- p_k; p_k_new[j] <- expit(mu_new[j])
          log_worth_k_new <- log_worth_k; log_worth_k_new[j] <- -theta_curr[k]*p_k_new[j]
          worth_k_new <- worth_k; worth_k_new[j] <- exp(log_worth_k_new[j])
          
          loglik_new <- unlist(lapply(ik,function(i){
            dbtlb_fast_log(rankings=rankings[i,,drop=FALSE],ratings=ratings[i,,drop=FALSE],
                           log_worth=log_worth_k_new,worth=worth_k_new,p=p_k_new,M=M,J=J,R=R[i])
          }))
          new_dens <- sum(loglik_new) + 
            sum(dnorm(mu_new-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
            dgamma(theta_curr[k],delta_theta[1],delta_theta[2],log=T)
          
          diff <- new_dens-curr_dens
          if(!is.nan(diff) & !is.infinite(diff) & log(runif(1)) < diff){
            accept_mu <- c(accept_mu,1)
            mu_curr[j,k] <- mu_new[j]
            p_k <- p_k_new
            log_worth_k <- log_worth_k_new
            worth_k <- worth_k_new
            curr_dens <- new_dens
            loglik_curr[ik,k] <- loglik_new
          }else{accept_mu <- c(accept_mu,0)}
        }
        ### Step 2(b)(ii) Update each theta_k
        theta_new <- rnorm(1,mean=theta_curr[k],sd=sigma_theta)
        if(theta_new >0){
          
          log_worth_k_new <- -theta_new*p_k
          worth_k_new <- exp(log_worth_k_new)
          
          loglik_new <- unlist(lapply(ik,function(i){
            dbtlb_fast_log(rankings=rankings[i,,drop=FALSE],ratings=ratings[i,,drop=FALSE],
                           log_worth=log_worth_k_new,worth=worth_k_new,p=p_k,M=M,J=J,R=R[i])
          }))
          
          new_dens <- sum(loglik_new) + 
            sum(dnorm(mu_curr[,k]-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
            dgamma(theta_new,delta_theta[1],delta_theta[2],log=T)
          
          diff <- new_dens-curr_dens
          if(!is.nan(diff) & !is.infinite(diff) & log(runif(1)) < diff){
            accept_theta <- c(accept_theta,1)
            theta_curr[k] <- theta_new
            log_worth_k <- log_worth_k_new
            worth_k <- worth_k_new
            curr_dens <- new_dens
            loglik_curr[ik,k] <- loglik_new
          }else{accept_theta <- c(accept_theta,0)}
        }else{accept_theta <- c(accept_theta,0)}
      }
      loglik_curr[-ik,k] <- unlist(lapply(setdiff(1:I,ik),function(i){
        dbtlb_fast_log(rankings=rankings[i,,drop=FALSE],ratings=ratings[i,,drop=FALSE],
                       log_worth=log_worth_k,worth=worth_k,p=p_k,M=M,J=J,R=R[i])
      }))
    }
    accept[t1,1:2] <- c(mean(accept_mu),mean(accept_theta))
    
    ### Step 2(b)(iii) Update sigma2_mu
    sigma2_mu_curr <- rinvgamma(1,beta_mu[1]+J*K_curr/2,beta_mu[2]+sum((mu_curr-mu_0j)^2)/2)
    
    
    ## Step 2(c)
    time3 <- Sys.time()
    ### Step 2(c)(i) Update K_curr
    probs <- unlist(lapply(Kplus_curr:(Kplus_curr+100),function(k){
      dpois(k-1,lambda,log=TRUE)+lfactorial(k)-lfactorial(k-Kplus_curr)+lgamma(gamma_curr*k)-lgamma(I+gamma_curr*k)
    }))
    K_curr <- sample(Kplus_curr:(Kplus_curr+100),1,prob = exp(probs-logSumExp(probs)))
    if(K_curr > ncol(alpha_iters)){
      mu_iters <- abind(mu_iters, array(NA, replace(dim(mu_iters), 2, K_curr-ncol(alpha_iters))), along = 2)
      theta_iters <- cbind(theta_iters,matrix(NA,nrow=T1,ncol=K_curr-ncol(alpha_iters)))
      class_iters <- abind(class_iters, array(0, replace(dim(class_iters), 3, K_curr-ncol(alpha_iters))), along = 3)
      alpha_iters <- cbind(alpha_iters,matrix(NA,nrow=T1,ncol=K_curr-ncol(alpha_iters)))
    }
    
    ### Step 2(c)(ii) Update gamma_curr
    gamma_new <- rnorm(1,gamma_curr,sd=sigma_gamma)
    dens_curr <- dgamma(gamma_curr,delta_gamma[1],delta_gamma[2],log=T)+lgamma(gamma_curr*K_curr)-
      lgamma(I+gamma_curr*K_curr)+sum(lgamma(Nk_curr+gamma_curr)-lgamma(gamma_curr))
    if(gamma_new>0){
      dens_new <-  dgamma(gamma_new, delta_gamma[1],delta_gamma[2],log=T)+lgamma(gamma_new*K_curr)-
        lgamma(I+gamma_new*K_curr )+sum(lgamma(Nk_curr+gamma_new )-lgamma(gamma_new))
      if(log(runif(1))<(dens_new-dens_curr)){
        accept[t1,3] <- 1
        gamma_curr <- gamma_new
      }
    }
    
    ## Step 2(d)
    time4 <- Sys.time()
    ### Step 2(d)(i) Create new (mu_k,theta_k) if K > K+
    if(K_curr > Kplus_curr){
      for(k in (Kplus_curr+1):K_curr){
        mu_curr <- cbind(mu_curr,rnorm(J,mean=mu_0j,sd=sqrt(sigma2_mu_curr)))
        theta_curr <- c(theta_curr,rgamma(1,delta_theta[1],delta_theta[2]))
        Nk_curr <- c(Nk_curr,0)
      }
      loglik_curr <- cbind(loglik_curr,matrix(NA,nrow=I,ncol=K_curr-Kplus_curr))
    }
    # fill in loglik_curr matrix
    which_k <- which(apply(loglik_curr,2,function(loglik){any(is.na(loglik))}))
    for(k in which_k){
      p_k <- expit(mu_curr[,k])
      log_worth_k <- -theta_curr[k]*p_k
      worth_k <- exp(log_worth_k)
      
      loglik_curr[,k] <- unlist(lapply(1:I,function(i){
        dbtlb_fast_log(rankings=rankings[i,,drop=FALSE],ratings=ratings[i,,drop=FALSE],log_worth=log_worth_k,worth=worth_k,p=p_k,M=M,J=J,R=R[i])
      }))
    }
    
    ### Step 2(d)(ii) Update alpha
    alpha_curr <- c(rdirichlet(1,gamma_curr+Nk_curr))
    time5 <- Sys.time()
    
    ## Store Current Values!
    Z_iters[t1,] <- Z_curr
    mu_iters[,1:K_curr,t1] <- mu_curr
    theta_iters[t1,1:K_curr] <- theta_curr
    gamma_iters[t1,] <- gamma_curr
    alpha_iters[t1,1:K_curr] <- alpha_curr
    K_iters[t1,] <- c(Kplus_curr,K_curr)
    sigma2_iters[t1,] <- c(sigma2_mu_curr)
    
    loglik_iters[t1,] <- sum(unlist(lapply(1:I,function(i){log(alpha_curr[Z_curr[i]])+loglik_curr[i,Z_curr[i]]})))+
      log(ddirichlet(alpha_curr,rep(gamma_curr,K_curr)))+
      dgamma(gamma_curr,delta_gamma[1],delta_gamma[2],log=T)+
      dpois(K_curr-1,lambda,log=T)+
      sum(dgamma(theta_curr,delta_theta[1],delta_theta[2],log=T))+
      sum(dnorm(mu_curr-mu_0j,mean=0,sd=sqrt(sigma2_mu_curr),log=T))+
      dinvgamma(sigma2_mu_curr,beta_mu[1],beta_mu[2],log=T)
    timing[t1,] <- difftime(c(time2,time3,time4,time5),c(time1,time2,time3,time4),units="secs")
    if(t1 %in% checkin_iters){
      text <- paste0("Status: ",which(t1 == checkin_iters)*10,"% complete. Current class sizes: ",paste0(Nk_curr,collapse=", "))
      print(text)
      print(paste0("Iteration: ",t1))
    }
  }
  # updating final class probabilities
  loglik_curr[is.infinite(loglik_curr)] <- -Inf
  unnorm_log_classprobs <- t(loglik_curr)+log(alpha_curr)
  classprobs <- exp(t(unnorm_log_classprobs)-apply(unnorm_log_classprobs,2,logSumExp))
  class_iters[t1,,1:K_curr] <- classprobs
  
  print("Done!")
  
  # Return Results!
  list(gamma=gamma_iters,sigma2=sigma2_iters,alpha=alpha_iters,
       mu=mu_iters,theta=theta_iters,Z=Z_iters,K=K_iters,loglik=loglik_iters,class=class_iters,
       accept=accept,timing=timing)
}
stephens_labelswitching <- function(btlb_mfmm_output,burn_prop=0,thin=1){
  total_iters <- dim(btlb_mfmm_output$mu)[3]
  keep_iters <- seq(ceiling(burn_prop*total_iters)+1,total_iters,by=thin)
  res_label <- stephens(btlb_mfmm_output$class[keep_iters,,])
  I <- dim(btlb_mfmm_output$Z)[2]
  for(iter in seq_along(keep_iters)){
    btlb_mfmm_output$mu[,,keep_iters[iter]] <- btlb_mfmm_output$mu[,res_label$permutations[iter,],keep_iters[iter]]
    btlb_mfmm_output$alpha[keep_iters[iter],] <- btlb_mfmm_output$alpha[keep_iters[iter],res_label$permutations[iter,]]
    btlb_mfmm_output$theta[keep_iters[iter],] <- btlb_mfmm_output$theta[keep_iters[iter],res_label$permutations[iter,]]
    
    new_Z <- btlb_mfmm_output$Z[iter,]
    for(i in 1:I){new_Z[i] <- res_label$permutations[iter,new_Z[i]]}
    btlb_mfmm_output$Z[iter,] <- new_Z
  }
  
  return(list(gamma=btlb_mfmm_output$gamma[keep_iters,,drop=FALSE],
              sigma2=btlb_mfmm_output$sigma2[keep_iters,],alpha=btlb_mfmm_output$alpha[keep_iters,],
              mu=btlb_mfmm_output$mu[,,keep_iters],nu=btlb_mfmm_output$nu[keep_iters,],
              theta=btlb_mfmm_output$theta[keep_iters,],Z=btlb_mfmm_output$Z[keep_iters,],
              K=btlb_mfmm_output$K[keep_iters,],loglik=btlb_mfmm_output$loglik[keep_iters,,drop=FALSE],
              class=btlb_mfmm_output$class,accept=btlb_mfmm_output$accept,timing=btlb_mfmm_output$timing))
}