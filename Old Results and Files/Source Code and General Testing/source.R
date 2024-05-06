#### Packages ####
library(gtools) # for dirichlet distribution functions
library(matrixStats) #for logSumExp function (helps with rounding errors)
library(fipp) #for calculating prior on Kplus
library(salso) #for cluster assignment plotting
library(reshape2)
library(ggplot2)
library(gridExtra)
library(dplyr)

#### BTL-B Functions ####

# Random Data Generation
rbtlb <- function(I,p,theta,M,R=length(p),seed=NULL){
  if(!is.null(seed)){set.seed(seed)}
  
  J <- length(p)
  X <- as.matrix(t(replicate(I,rbinom(n=J,size=M,prob=p))))
  Pi <- as.matrix(t(replicate(I,{sample(1:J,R,prob=exp(-theta*p))})))
  
  return(list(X=X,Pi=Pi))
}

# Density Functions
dbtlb <- function(Pi,X,p,theta,M,log=FALSE,Pi_full=NULL){
  
  #check data format (Pi is checked by dbtl function called herein)
  if(is.vector(X)){X <- matrix(X,nrow=1)}
  if(!is.matrix(X)){stop("X must be a matrix of ratings")}
  if(length(p)!=ncol(X)){stop("length(p) must equal ncol(X)")}
  
  worth <- exp(-theta*p)
  logd <- dbtl(Pi=Pi,worth=worth,log=T,Pi_full=Pi_full)+
    sum(apply(X,1,function(x){dbinom(x=x,size=M,prob=p,log=T)}),na.rm=T)
  
  if(log){
    return(logd)
  }else{
    return(exp(logd))
  }
}
dbtl <- function(Pi,worth,log=FALSE,Pi_full=NULL){
  if(is.vector(Pi)){Pi <- matrix(Pi,nrow=1)}
  if(is.vector(Pi_full)){Pi_full <- matrix(Pi_full,nrow=1)}
  if(!is.matrix(Pi)){stop("Pi must be a matrix of (partial) rankings")}
  if(!is.null(Pi_full)){if(nrow(Pi)!=nrow(Pi_full)){stop("nrow(Pi) must equal nrow(Pi_full).")}}
  if(length(worth)<max(c(0,Pi,Pi_full),na.rm=T)){stop("Ensure one worth parameter for every object in Pi and/or Pi_full")}
  
  if(sum(apply(Pi,1,function(pi){!all(is.na(pi))}))==0){ #if all rankings empty, return likelihood of 1 (convenience)
    if(log){return(0)}else{return(1)}
  }
  
  if(!is.null(Pi_full)){
    logd <- sum(unlist(lapply(1:nrow(Pi),function(i){
      pi <- na.exclude(Pi[i,])
      if(length(pi)==0){return(0)
      }else{
        num <- worth[pi]
        denom_terms <- na.exclude(worth[Pi_full[i,]])
        return(sum(log(num))-sum(unlist(lapply(1:length(num),function(j){log(sum(denom_terms[j:length(denom_terms)]))}))))
      }
    })))
  }else{
    logd <- sum(apply(Pi,1,function(pi){
      if(length(na.exclude(pi))==0){return(0)
      }else{
        pi <- na.exclude(pi)
        num <- worth[pi]
        denom_terms <- c(worth[pi],worth[setdiff(1:length(worth),pi)])
        return(sum(log(num)-log(rev(cumsum(rev(denom_terms)))[1:length(pi)])))
      }
    }))
  }
  
  if(log){return(logd)
  }else{return(exp(logd))}
}

# MAP Estimation
btlb_map <- function(Pi,X,M,Pi_full=NULL,K,xi1,xi2,a,b,gamma1,gamma2,
                     tol=1,maxit=50,verbose=TRUE,seed=NULL){
  
  
  # function to calculate objective function
  get_obj <- function(ptheta,pi,gamma){
    lik_terms <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      pk <- ptheta[1:J,k]
      thetak <- ptheta[J+1,k]
      worthk <- exp(-thetak*pk)
      lik_terms[,k] <- log(pi[k])+unlist(lapply(1:I,function(i){
        dbtl(Pi=Pi[i,],worth=worthk,Pi_full=Pi_full[i,],log=T)+sum(dbinom(x=X[i,],size=M,prob=pk,log=T),na.rm=T)
      }))
    }
    lik_term <- sum(apply(lik_terms,1,logSumExp))
    prior_term <- dgamma(gamma,xi1,xi2,log=TRUE)+log(ddirichlet(pi,rep(gamma,K)))+
      sum(dbeta(ptheta[1:J,],a,b,log=TRUE))+sum(dgamma(ptheta[J+1,],gamma1,gamma2,log=TRUE))
    currobj <- lik_term+prior_term
    return(currobj)
  }
  
  # model constants
  I <- nrow(X)
  J <- ncol(X)
  
  # initializers
  if(!is.null(seed)){set.seed(seed)}
  gamma <- rgamma(1,xi1,xi2)
  pi <- rep(1/K,K)
  ptheta <- matrix(c(rbeta(J*K,a,b),rgamma(K,gamma1,gamma2)),byrow=TRUE,nrow=J+1,ncol=K)
  Z <- matrix(NA,nrow=I,ncol=K)
  
  #get current loglikelihood
  curr_obj <- get_obj(ptheta,pi,gamma)
  
  currdiff <- Inf
  iter <- 0
  while(currdiff>tol & iter<=maxit){
    iter <- iter + 1
    
    # E-Step: Update Z
    t1 <- Sys.time()
    for(i in 1:I){
      lprobs <- unlist(lapply(1:K,function(k){
        log(pi[k])+dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=ptheta[1:J,k],
                         theta=ptheta[J+1,k],M=M,log=T,Pi_full=Pi_full[i,])
      }))
      Z[i,] <- exp(lprobs-logSumExp(lprobs))
    }

    
    # M-Step:
    
    ## Update pi
    t2 <- Sys.time()
    pi <- (gamma-1+apply(Z,2,sum))/(I - K + gamma*K)
    if(round(sum(pi),5)!=1){print("Error! Something with pi is wrong")}
    
    ## Update gamma
    t3 <- Sys.time()
    gamma <- optimize(f=function(gamma){lgamma(gamma*K)-K*lgamma(gamma)+(xi1-1)*log(gamma)+gamma*(sum(log(pi))-xi2)},
                      interval=c(1e-8,qgamma(1-1e-8,xi1,xi2)),maximum=TRUE)$maximum
    
    ## Update (p,theta)
    t4 <- Sys.time()
    for(k in 1:K){
      #print(k)
      zhat_k <- Z[,k]
      # which_vals <- which(zhat_k>.001)
      # if(length(which_vals)==0){which_vals <- 1:I}
      obj <- function(par){
        worthk <- exp(-par[1:J]*par[J+1])
        lik_term <- sum(unlist(lapply(1:I,function(i){
          zhat_k[i]*(dbtl(Pi=Pi[i,],worth=worthk,log=T,Pi_full=Pi_full[i,])+sum(X[i,]*log(par[1:J])+(M-X[i,])*log(1-par[1:J]),na.rm=T))
        })))
        prior_term <- sum(dbeta(par[1:J],a,b,log=TRUE))+dgamma(par[J+1],gamma1,gamma2,log=TRUE)
        return(-lik_term-prior_term)
      }
      ptheta[,k] <- optim(ptheta[,k],obj,method="L-BFGS-B",lower=rep(0.00001,J+1),upper=c(rep(1-.00001,J),qgamma(.9999,gamma1,gamma2)))$par
    }
    
    

   
    # Iteration Updates
    t5 <- Sys.time()
    new_obj <- get_obj(ptheta,pi,gamma)
    currdiff <- new_obj-curr_obj
    curr_obj <- new_obj
    t6 <- Sys.time()
    
    if(verbose){
      print(paste0("End of Iteration ",iter))
      print("Time to update Z, pi, gamma,ptheta,objective: ")
      print(round(c(t2-t1,t3-t2,t4-t3,t5-t4,t6-t5),2))
      print(paste0("Difference: ",round(currdiff,3)))
    }
  }
  
  return(list(p=ptheta[1:J,],theta=ptheta[J+1,],pihat=pi,gammahat=gamma,Z=Z,obj=curr_obj))
  
}

# MFM Estimation
btlb_mfm <- function(Pi,X,M,Pi_full=NULL,lambda,a,b,gamma1,gamma2,xi1,xi2,
                     startK = 1, mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                     max_iters = 100, mh_iters = 5, burn = 0.5,thin = 2,seed = NULL){
  
  print("Initializing Chain")
  if(!is.null(seed)){set.seed(seed)}
  
  ## Determine Constants
  I <- nrow(X)
  J <- ncol(X)
  R <- ncol(Pi)
  
  ## Set Up Data Storage
  which_keep <- seq(from=round(burn*max_iters*mh_iters),to=max_iters * mh_iters,by=thin)
  K_all <- c()
  Kplus_all <- c()
  gamma_all <- c()
  pi_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  p_all <- array(NA,dim=c(J,I,length(which_keep)))
  theta_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  Z_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  accept_p <- c()
  accept_theta <- c()
  accept_gamma <- c()
  
  ## Initialize at Random
  K <- startK
  gamma <- rgamma(1,xi1,xi2)
  pi <- c(rdirichlet(1,rep(gamma,K)))
  ptheta <- matrix(NA,nrow=J+1,ncol=K)
  for(k in 1:K){ptheta[,k] <- c(rbeta(J,a,b),rgamma(1,gamma1,gamma2))}
  Z <- rep(NA,I)
  
  print("Starting Chain")
  iter <- 1
  progress_iters <- round(seq(0,max(which_keep),length=11))[-1]
  while(iter <= max(which_keep)){
    
    ## Step 1: Update Labels
    probs <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      pk <- ptheta[1:J,k]
      thetak <- ptheta[J+1,k]
      worthk <- exp(-thetak*pk)
      probs[,k] <- exp(log(pi[k])+unlist(lapply(1:I,function(i){
        dbtl(Pi=Pi[i,],worth=worthk,Pi_full=Pi_full[i,],log=T)+sum(dbinom(x=X[i,],size=M,prob=pk,log=T),na.rm=T)
      })))
    }
    Z <- apply(probs,1,function(prob){sample.int(K,1,prob=prob)})
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)})) #update Nk, Kplus, and reorder
    Kplus <- sum(Nk>0)
    ptheta <- as.matrix(ptheta[,which(Nk>0)])
    Ztmp <- Z
    for(k in 1:K){Ztmp[which(Z==k)] <- which(c(which(Nk>0),setdiff(1:K,which(Nk>0)))==k)}
    Z <- Ztmp
    rm(Ztmp)
    
    
    ## Step 2: Update non-empty component parameters
    for(mh_iter in 0:(mh_iters-1)){
      for(k in 1:Kplus){
        whichk <- which(Z == k)
        constant1 <- a+apply(matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        constant2 <- b+apply(M-matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        Pi_mat <- matrix(Pi[whichk,],nrow=length(whichk))
        if(is.null(Pi_full)){Pi_full_mat <- NULL}else{Pi_full_mat <- matrix(Pi_full[whichk,],nrow=length(whichk))}
        
        for(j in 1:J){ #update each p_jk
          prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)
          while(prop_pjk<=0 | prop_pjk>=1){prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)}
          prop_p <- ptheta[1:J,k]
          prop_p[j] <- prop_pjk
          curr_worth <- exp(-ptheta[J+1,k]*ptheta[1:J,k])
          prop_worth <- exp(-ptheta[J+1,k]*prop_p)
          
          logprob_prop <- dbtl(Pi=Pi_mat,worth=prop_worth,log=T,Pi_full=Pi_full_mat)+
            (constant1[j])*log(prop_pjk)+(constant2[j])*log(1-prop_pjk)
          logprob_curr <- dbtl(Pi=Pi_mat,worth=curr_worth,log=T,Pi_full=Pi_full_mat)+
            (constant1[j])*log(ptheta[j,k])+(constant2[j])*log(1-ptheta[j,k])
          
          u <- runif(1)
          if(log(u) <  logprob_prop-logprob_curr){
            accept_p <- c(accept_p,1)
            ptheta[j,k] <- prop_pjk
          }else{accept_p <- c(accept_p,0)}
        }
        #update theta_k
        prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)
        while(prop_thetak<=0){prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)}
        prop_worth <- exp(-prop_thetak*ptheta[1:J,k])
        curr_worth <- exp(-ptheta[J+1,k]*ptheta[1:J,k])
        
        logprob_prop <- dbtl(Pi_mat,worth=prop_worth,log=T,Pi_full=Pi_full_mat)+
          (gamma1-1)*log(prop_thetak)-gamma2*prop_thetak
        logprob_curr <- dbtl(Pi_mat,worth=curr_worth,log=T,Pi_full=Pi_full_mat)+
          (gamma1-1)*log(ptheta[J+1,k])-gamma2*ptheta[J+1,k]
        
        u <- runif(1)
        if(log(u) <  logprob_prop-logprob_curr){
          accept_theta <- c(accept_theta,1)
          ptheta[J+1,k] <- prop_thetak
        }else{accept_theta <- c(accept_theta,0)}
      }
      if(iter+mh_iter %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        p_all[,1:Kplus,it] <- ptheta[1:J,] #update only the first Kplus parameters,
        theta_all[it,1:Kplus] <- ptheta[J+1,] #(will update the rest later!)
      }
    }
    
    
    ## Step 3: Update K and gamma
    
    #update K
    probs <- unlist(lapply(Kplus:(Kplus+100),function(k){
      dpois(k-1,lambda,log=T)+lfactorial(k)-lfactorial(k-Kplus)+lgamma(gamma*k)-lgamma(I+gamma*k)
    }))
    K <- sample(Kplus:(Kplus+100),1,prob = exp(probs-logSumExp(probs)))
    
    #update gamma
    prop_gamma <- rnorm(1,gamma,mh_gamma)
    while(prop_gamma<=0){prop_gamma <- rnorm(1,gamma,mh_gamma)}
    logprob_prop <- dgamma(prop_gamma,xi1,xi2,log=TRUE)+lgamma(prop_gamma*K)-lgamma(I+prop_gamma*K)+
      sum(lgamma(Nk[1:Kplus]+prop_gamma)-lgamma(prop_gamma))
    logprob_curr <- dgamma(gamma,3,2,log=TRUE)+lgamma(gamma*K)-lgamma(I+gamma*K)+
      sum(lgamma(Nk[1:Kplus]+gamma)-lgamma(gamma))
    u <- runif(1)
    if(log(u) <  logprob_prop-logprob_curr){
      accept_gamma <- c(accept_gamma,1)
      gamma <- prop_gamma
    }else{accept_gamma <- c(accept_gamma,0)}
    
    ## Step 4: Update empty components
    if(K>Kplus){
      for(k in (Kplus+1):K){ptheta <- cbind(ptheta,c(rbeta(J,a,b),rgamma(1,gamma1,gamma2)))}
    }
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)}))
    pi <- c(rdirichlet(1,gamma+Nk))
    
    if(length(pi)!=K | length(Nk)!=K | ncol(ptheta)!=K ){stop("something wrong!")}
    if(Kplus > K){stop("something wrong with Kplus!")}
    
    ## Save Values and Update Counter
    for(mh_iter in 0:(mh_iters-1)){
      if((iter+mh_iter) %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        K_all <- c(K_all,K)
        Kplus_all <- c(Kplus_all,Kplus)
        gamma_all <- c(gamma_all,gamma)
        pi_all[it,1:K] <- pi
        Z_all[it,] <- Z
        if(K>Kplus){
          p_all[,(Kplus+1):K,it] <- ptheta[1:J,(Kplus+1):K] #update only the first Kplus parameters,
          theta_all[it,(Kplus+1):K] <- ptheta[J+1,(Kplus+1):K] #(will update the rest later!)
        }
      }}
    
    iter <- iter + mh_iters
    if(length(progress_iters)>0 & iter >= progress_iters[1]){
      print(paste0((11-length(progress_iters))*10,"% Complete: Iteration ",iter-1," out of ",max(which_keep)))
      print(paste0("Current K = ",K,"; Kplus = ",Kplus))
      progress_iters <- progress_iters[-1]
    }
  }
  print(paste0("Done! Saving ",length(which_keep)," estimate iterations after burning/thinning"))
  return(list(Z=Z_all,p=p_all[,1:max(K_all),],theta=theta_all[,1:max(K_all)],
              pi=pi_all[,1:max(K_all)],K=K_all,Kplus=Kplus_all,gamma=gamma_all,
              accept_p=(cumsum(accept_p)/1:length(accept_p))[seq(round(length(accept_p)*burn),length(accept_p),by=thin)],
              accept_theta=(cumsum(accept_theta)/1:length(accept_theta))[seq(round(length(accept_theta)*burn),length(accept_theta),by=thin)],
              accept_gamma=(cumsum(accept_gamma)/1:length(accept_gamma))[seq(round(length(accept_gamma)*burn),length(accept_gamma),by=thin)],
              max_iters=max_iters,mh_iters=mh_iters,burn=burn,thin=thin,seed=seed))
}

# Finite Mixture (Fixed K) Estimation
btlb_fm <- function(Pi,X,M,Pi_full=NULL,K,a,b,gamma1,gamma2,xi1,xi2,
                    mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                    max_iters = 100, mh_iters = 5, burn = 0.5, thin = 2, seed = NULL){
  
  print("Initializing Chain")
  if(!is.null(seed)){set.seed(seed)}
  
  ## Determine Constants
  I <- nrow(X)
  J <- ncol(X)
  R <- ncol(Pi)
  
  ## Set Up Data Storage
  which_keep <- seq(from=round(burn*max_iters*mh_iters),to=max_iters * mh_iters,by=thin)
  gamma_all <- c()
  pi_all <- matrix(NA,nrow=length(which_keep),ncol=K)
  p_all <- array(NA,dim=c(J,K,length(which_keep)))
  theta_all <- matrix(NA,nrow=length(which_keep),ncol=K)
  Z_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  accept_p <- c()
  accept_theta <- c()
  accept_gamma <- c()
  
  ## Initialize at Random
  gamma <- rgamma(1,xi1,xi2)
  pi <- c(rdirichlet(1,rep(gamma,K)))
  ptheta <- matrix(NA,nrow=J+1,ncol=K)
  for(k in 1:K){ptheta[,k] <- c(rbeta(J,a,b),rgamma(1,gamma1,gamma2))}
  Z <- rep(NA,I)
  
  print("Starting Chain")
  iter <- 1
  progress_iters <- round(seq(0,max(which_keep),length=11))[-1]
  while(iter <= max(which_keep)){
    
    ## Step 1: Update Labels
    probs <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      pk <- ptheta[1:J,k]
      thetak <- ptheta[J+1,k]
      worthk <- exp(-thetak*pk)
      probs[,k] <- exp(log(pi[k])+unlist(lapply(1:I,function(i){
        dbtl(Pi=Pi[i,],worth=worthk,Pi_full=Pi_full[i,],log=T)+sum(dbinom(x=X[i,],size=M,prob=pk,log=T),na.rm=T)
      })))
    }
    Z <- apply(probs,1,function(prob){sample.int(K,1,prob=prob)})
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)}))
    
    
    ## Step 2: Update component parameters
    for(mh_iter in 0:(mh_iters-1)){
      for(k in 1:K){
        whichk <- which(Z == k)
        if(length(whichk)==0 ){whichk <- 1:I}
        constant1 <- a+apply(matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        constant2 <- b+apply(M-matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        Pi_mat <- matrix(Pi[whichk,],nrow=length(whichk))
        if(is.null(Pi_full)){Pi_full_mat <- NULL}else{Pi_full_mat <- matrix(Pi_full[whichk,],nrow=length(whichk))}
        
        for(j in 1:J){ #update each p_jk
          prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)
          while(prop_pjk<=0 | prop_pjk>=1){prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)}
          prop_p <- ptheta[1:J,k]
          prop_p[j] <- prop_pjk
          curr_worth <- exp(-ptheta[J+1,k]*ptheta[1:J,k])
          prop_worth <- exp(-ptheta[J+1,k]*prop_p)
          
          logprob_prop <- dbtl(Pi=Pi_mat,worth=prop_worth,log=T,Pi_full=Pi_full_mat)+
            (constant1[j])*log(prop_pjk)+(constant2[j])*log(1-prop_pjk)
          logprob_curr <- dbtl(Pi=Pi_mat,worth=curr_worth,log=T,Pi_full=Pi_full_mat)+
            (constant1[j])*log(ptheta[j,k])+(constant2[j])*log(1-ptheta[j,k])
          
          u <- runif(1)
          if(log(u) <  logprob_prop-logprob_curr){
            accept_p <- c(accept_p,1)
            ptheta[j,k] <- prop_pjk
          }else{accept_p <- c(accept_p,0)}
        }
        #update theta_k
        prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)
        while(prop_thetak<=0){prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)}
        prop_worth <- exp(-prop_thetak*ptheta[1:J,k])
        curr_worth <- exp(-ptheta[J+1,k]*ptheta[1:J,k])
        
        logprob_prop <- dbtl(Pi_mat,worth=prop_worth,log=T,Pi_full=Pi_full_mat)+
          (gamma1-1)*log(prop_thetak)-gamma2*prop_thetak
        logprob_curr <- dbtl(Pi_mat,worth=curr_worth,log=T,Pi_full=Pi_full_mat)+
          (gamma1-1)*log(ptheta[J+1,k])-gamma2*ptheta[J+1,k]
        
        u <- runif(1)
        if(log(u) <  logprob_prop-logprob_curr){
          accept_theta <- c(accept_theta,1)
          ptheta[J+1,k] <- prop_thetak
        }else{accept_theta <- c(accept_theta,0)}
      }
      
      if(iter+mh_iter %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        p_all[,,it] <- ptheta[1:J,] 
        theta_all[it,] <- ptheta[J+1,]
      }
    }
    
    ## Step 3: Update gamma
    prop_gamma <- rnorm(1,gamma,mh_gamma)
    while(prop_gamma<=0){prop_gamma <- rnorm(1,gamma,mh_gamma)}
    logprob_prop <- dgamma(prop_gamma,xi1,xi2,log=TRUE)+lgamma(prop_gamma*K)-lgamma(I+prop_gamma*K)+
      sum(lgamma(Nk[1:K]+prop_gamma)-lgamma(prop_gamma))
    logprob_curr <- dgamma(gamma,3,2,log=TRUE)+lgamma(gamma*K)-lgamma(I+gamma*K)+
      sum(lgamma(Nk[1:K]+gamma)-lgamma(gamma))
    u <- runif(1)
    if(log(u) <  logprob_prop-logprob_curr){
      accept_gamma <- c(accept_gamma,1)
      gamma <- prop_gamma
    }else{accept_gamma <- c(accept_gamma,0)}
    
    ## Step 4: Update pi
    pi <- c(rdirichlet(1,gamma+Nk))
    if(length(pi)!=K | length(Nk)!=K | ncol(ptheta)!=K ){stop("something wrong!")}
    
    ## Save Values and Update Counter
    for(mh_iter in 0:(mh_iters-1)){
      if((iter+mh_iter) %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        gamma_all <- c(gamma_all,gamma)
        pi_all[it,1:K] <- pi
        Z_all[it,] <- Z
      }}
    
    iter <- iter + mh_iters
    if(length(progress_iters)>0 & iter >= progress_iters[1]){
      print(paste0((11-length(progress_iters))*10,"% Complete: Iteration ",iter-1," out of ",max(which_keep)))
      progress_iters <- progress_iters[-1]
    }
  }
  print(paste0("Done! Saving ",length(which_keep)," estimate iterations after burning/thinning"))
  return(list(Z=Z_all,p=p_all,theta=theta_all,pi=pi_all,gamma=gamma_all,
              accept_p=(cumsum(accept_p)/1:length(accept_p))[seq(round(length(accept_p)*burn),length(accept_p),by=thin)],
              accept_theta=(cumsum(accept_theta)/1:length(accept_theta))[seq(round(length(accept_theta)*burn),length(accept_theta),by=thin)],
              accept_gamma=(cumsum(accept_gamma)/1:length(accept_gamma))[seq(round(length(accept_gamma)*burn),length(accept_gamma),by=thin)],
              max_iters=max_iters,mh_iters=mh_iters,burn=burn,thin=thin,seed=seed))
  
  
}

# Goodness-of-Fit Functions
sample_postpred <- function(posterior,post_samples=10,reps=5){
  J <- nrow(posterior$p)
  X <- matrix(NA,nrow=0,ncol=J)
  if(is.vector(posterior$pi)){posterior$pi <- matrix(posterior$pi)}
  if(is.matrix(posterior$p)){posterior$p <- array(posterior$p,dim=c(nrow(posterior$p),1,ncol(posterior$p)))}
  if(is.vector(posterior$theta)){posterior$theta <- matrix(posterior$theta)}
  which_k <- 1:nrow(posterior$pi)
  for(iteration in sample(which_k,post_samples,replace=TRUE)){
    K <- sum(!is.na(posterior$pi[iteration,]))
    X <- rbind(X,t(replicate(reps,{rbinom(J,M,posterior$p[,sample(1:K,1,prob=posterior$pi[iteration,1:K],replace=T),iteration])})))
  }
  Pi <- matrix(NA,nrow=0,ncol=J)
  for(iteration in sample(which_k,post_samples,replace=TRUE)){
    K <- sum(!is.na(posterior$pi[iteration,]))
    Pi <- rbind(Pi,t(replicate(reps,{
      z <- sample(1:K,1,prob=posterior$pi[iteration,1:K],replace=T)
      sample(1:J,J,prob=exp(-posterior$theta[iteration,z]*posterior$p[,z,iteration]))
    })))
  }
  return(list(X=X,Pi=Pi))
}
get_gof <- function(posterior,post_samples=1000,reps=nrow(X),X,Pi,Pi_full=NULL,seed=NULL){
  
  if(!is.null(seed)){set.seed(seed)}
  
  #sample from posterior predictive
  postpred <- t(replicate(post_samples,{
    data <- sample_postpred(posterior,1,reps)
    c(apply(data$X,2,mean),apply(data$X,2,var))
  }))
  J <- dim(posterior$p)[1]
  
  #calculate summary statistics for rating mean and variance
  postpred_mean <- melt(postpred[,1:J])
  postpred_var <- melt(postpred[,(J+1):(2*J)])
  mean_data <- apply(X,2,function(x){mean(x,na.rm=T)})
  var_data <- apply(X,2,function(x){var(x,na.rm=T)})
  
  #calculate summary statistics for pairwise comparisons
  pairs_data <- matrix(NA,nrow=0,ncol=2)
  for(i in 1:(J-1)){for(j in (i+1):J){pairs_data <- rbind(pairs_data,c(i,j))}}
  obs <- t(apply(pairs_data,1,function(ij){
    i <- ij[1]
    j <- ij[2]
    
    if(!is.null(Pi_full)){
      which <- which(apply(Pi_full,1,function(pi){i %in% pi & j %in% pi}))
    }else{which <- 1:nrow(Pi)}
    if(length(which) == 0){return(c(NA,0))}
    
    rank_i <- apply(Pi[which,],1,function(pi){if(i %in% pi){return(which(pi==i))}else{return(NA)}})
    rank_j <- apply(Pi[which,],1,function(pi){if(j %in% pi){return(which(pi==j))}else{return(NA)}})
    
    pairs <- unlist(lapply(1:length(rank_i),function(index){
      if(!is.na(rank_i[index]) & !is.na(rank_j[index])){return(ifelse(rank_i[index] < rank_j[index],1,0))}
      if(!is.na(rank_i[index]) & is.na(rank_j[index])){return(1)}
      if(is.na(rank_i[index]) & !is.na(rank_j[index])){return(0)}
      if(is.na(rank_i[index]) & is.na(rank_j[index])){return(NA)}
    }))
    
    return(c(mean(pairs,na.rm=T),sum(!is.na(pairs))))
  }))
  postpred <- sample_postpred(posterior,post_samples,reps)
  post <- apply(pairs_data,1,function(ij){
    i <- ij[1]
    j <- ij[2]
    rank_i <- apply(postpred$Pi,1,function(pi){which(pi==i)})
    rank_j <- apply(postpred$Pi,1,function(pi){which(pi==j)})
    mean(rank_i < rank_j)
  })
  pairs_data <- as.data.frame(cbind(pairs_data,obs,post))
  names(pairs_data) <- c("i","j","obs","num_comparisons","post")
  
  #create plots
  p1<-ggplot(postpred_mean,aes(factor(Var2),value))+geom_violin()+
    geom_point(data=data.frame(Var2=1:J,value=mean_data),aes(factor(Var2),value),color="red")+
    ylim(c(0,max(postpred_mean[,"value"],mean_data,40)))+
    xlab("Proposal")+ylab("Posterior Predictive Mean")+ggtitle("Ratings: Mean")
  p2<-ggplot(postpred_var,aes(factor(Var2),value))+geom_violin()+
    ylim(c(0,max(postpred_var[,"value"],var_data)))+
    geom_point(data=data.frame(Var2=1:J,value=var_data),aes(factor(Var2),value),color="red")+
    xlab("Proposal")+ylab("Posterior Predictive Variance")+ggtitle("Ratings: Variance")
  p3<-ggplot(pairs_data,aes(x=obs,y=post))+geom_point()+
    geom_abline(slope=1,intercept=0,linetype=2)+
    xlim(c(0,1))+ylim(c(0,1))+
    xlab("Observed Probabilities")+
    ylab("Posterior Predictive Expected Probabilities")+
    ggtitle("Rankings: Pairwise Probabilities")
  
  list(mean_data=mean_data,var_data=var_data,pairs_data=pairs_data,p1=p1,p2=p2,p3=p3)
}


#### BTL-B Examples ####

if(FALSE){ #change to "TRUE" to run
  
  # data generation
  set.seed(1)
  M <- 4
  dat1 <- rbtlb(I=20,p=runif(10),theta=10,M=M)
  dat2 <- rbtlb(I=10,p=runif(10),theta=20,M=M)
  X <- rbind(dat1$X,dat2$X)
  Pi <- rbind(dat1$Pi,dat2$Pi)
  rm(dat1,dat2)
  
  # density functions
  dbtlb(Pi=Pi,X=X,p=runif(ncol(X)),theta=10,M=M,log=TRUE)
  dbtl(Pi=Pi,worth=exp(-10*runif(ncol(X))),log=TRUE)
  
  # map
  map <- btlb_map(Pi=Pi,X=X,M=M,Pi_full=NULL,K=2,
                  xi1=2,xi2=3,a=2,b=2,gamma1=10,gamma2=.5,
                  tol=.01,maxit=50,verbose=TRUE,seed=1)
  
  # calculation of prior on K+
  pmfstatic2 <- nClusters(Kplus=1:5,N=nrow(X),type="static",gamma=1,maxK=20)
  dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = 1))
  ggplot(data.frame(K=1:length(dens),Mass=dens),aes(K,Mass))+geom_line()+geom_point()+
    ylab("Prior Mass on K+")+xlab("K+")
  
  # mfm
  mfm <- btlb_mfm(Pi=Pi,X=X,M=M,Pi_full=NULL,lambda=1,a=2,b=2,gamma1=10,gamma2=.5,xi1=2,xi2=2,
                  startK = 1, mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                  max_iters = 200, mh_iters = 5, burn = 0.5,thin = 2,seed = 1)

  gof <- get_gof(posterior=mfm,post_samples=500,reps=5,X=X,Pi=Pi,Pi_full=NULL)
  #grid.arrange(gof$p1,gof$p2,gof$p3,nrow=1)
  
  # fixed K
  fm <- btlb_fm(Pi=Pi,X=X,M=M,Pi_full=NULL,K=2,
                a=2,b=2,gamma1=10,gamma2=.5,xi1=2,xi2=2,
                mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                max_iters = 200, mh_iters = 5, burn = 0.5, thin = 2, seed = 1)
}

#### MB Functions ####

# Random Data Generation
rmallow <- function(I,pi0,theta,R=length(pi0)){
  tmp<-function(pi0,theta){
    J<-length(pi0)
    Vj <- pi1 <- c()
    for(j in 1:(J-1)){
      probs <- exp(-theta*(0:(J-j)))/sum(exp(-theta*(0:(J-j))))
      Vj <- c(Vj,sample.int(J-j+1,size=1,prob=probs)-1)
    }
    Vj <- c(Vj, 0)
    for(j in 1:J){
      pi1 <- c(pi1, pi0[Vj[j]+1])
      pi0 <- setdiff(pi0, pi0[Vj[j]+1])
    }
    return(pi1)
  }
  Pi <- t(replicate(I, tmp(pi0, theta)))[,1:R]
  if(I==1){Pi <- matrix(Pi,nrow=1,ncol=R)}
  return(Pi)
}
rmb <- function(I,p,pi0=NULL,theta,M,R=length(p)){
  J <- length(p)
  X <- as.matrix(t(replicate(I, rbinom(n = J, size = M, prob = p))))
  
  ties <- length(setdiff(1:J,rank(p,ties.method="min")))>0
  if(ties & is.null(pi0)){stop("Ties present in p. Set pi0 to break ties")}
  if(ties & !is.null(pi0)){if(any(diff(p[pi0])<0)){stop("pi0 is incompatible with p")}}
  if(!ties & is.null(pi0)){pi0 <- order(p)}
  if(!ties & !is.null(pi0)){if(any(order(p)!=pi0)){stop("pi0 is incompatible with p")}}
  
  Pi <- rmallow(I=I,pi0=pi0,theta=theta,R=R)
  return(list(X=X,Pi=Pi))
}

# Density Functions
dmb <- function (Pi, X, p, pi0 = NULL, theta, M, log = FALSE){
  
  #check data format (Pi is checked by dmallow function called herein)
  if(is.vector(X)){X <- matrix(X,nrow=1)}
  if(!is.matrix(X)){stop("X must be a matrix of ratings")}
  if(length(p)!=ncol(X)){stop("length(p) must equal ncol(X)")}
  
  if(length(p)!=ncol(X)){stop("p must equal ncol(X)")}
  ties <- length(setdiff(1:ncol(X),rank(p,ties.method="min")))>0
  if(ties & is.null(pi0)){stop("Ties present in p. Set pi0 to break ties")}
  if(ties & !is.null(pi0)){if(any(diff(p[pi0])<0)){stop("pi0 is incompatible with p")}}
  if(!ties & is.null(pi0)){pi0 <- order(p)}
  if(!ties & !is.null(pi0)){if(any(order(p)!=pi0)){stop("pi0 is incompatible with p")}}
  
  logd <- dmallow(Pi=Pi,pi0=pi0,theta=theta,log=T) + 
    sum(apply(X, 1, function(x){dbinom(x, M, p, log = T)}), na.rm=T)
  
  if(log){return(logd)}else{return(exp(logd))}
}
dmallow <- function(Pi, pi0, theta, log = FALSE){
  
  if(is.vector(Pi)){Pi <- matrix(Pi,nrow=1)}
  if(!is.matrix(Pi)){stop("Pi must be a matrix of (partial) rankings")}
  
  if(sum(apply(Pi,1,function(pi){!all(is.na(pi))}))==0){ #if all rankings empty, return likelihood of 1 (convenience)
    if(log){return(0)}else{return(1)}
  }
  
  J <- length(pi0)
  R <- apply(Pi, 1, function(pi){length(na.exclude(pi))})
  R[which(R == J - 1)] <- J
  logd <- -theta * sum(apply(Pi, 1, function(pi){kendall(pi, pi0)})) - 
    sum(unlist(lapply(R, function(r){psi(theta, J, r, log = T)})))
  
  if(log){return(logd)}else{return(exp(logd))}
}
kendall <- function(pi, pi0) {
  pi <- na.exclude(as.vector(pi))
  pi0 <- as.vector(pi0)
  R <- length(pi)
  J <- length(pi0)
  if (length(setdiff(1:J, pi0)) != 0) {
    stop("pi0 must be a complete ranking")
  }
  if (any(pi > J, na.rm = T)) {
    stop("pi cannot contain items not in pi0")
  }
  if (R > J) {
    stop("R must be <= J")
  }
  dist <- 0
  for (r in 1:R) {
    dist <- dist + (which(pi0 == pi[r]) - 1)
    pi0 <- setdiff(pi0, pi[r])
  }
  return(dist)
}
psi <- function(theta, J, R, log = FALSE){
  if (R > J) {
    stop("R must be <= J")
  }
  if (R <= 0 | J <= 0) {
    stop("R and J must be positive integers")
  }
  if (theta <= 0) {
    stop("theta must be >=0")
  }
  log_psi <- sum(log((1 - exp(-theta * (J - 1:R + 1)))) - log((1 - 
                                                                 exp(-theta))))
  if (log) {
    return(log_psi)
  }
  else {
    return(exp(log_psi))
  }
}

# MAP Estimation
mb_map <- function(Pi,X,M,K,gamma,a,b,gamma1,gamma2,
                     tol=1,maxit=50,verbose=TRUE,seed=NULL){
  
  
  # function to calculate objective function
  get_obj <- function(ptheta,pi){
    lik_terms <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      pk <- ptheta[1:J,k]
      thetak <- ptheta[J+1,k]
      lik_terms[,k] <- log(pi[k])+unlist(lapply(1:I,function(i){
        dmb(Pi=Pi[i,],X=X[i,],p=pk,theta=thetak,M=M,log=T)
      }))
    }
    lik_term <- sum(apply(lik_terms,1,logSumExp))
    prior_term <- log(ddirichlet(pi,rep(gamma,K)))+sum(dbeta(ptheta[1:J,],a,b,log=TRUE))+
      sum(dgamma(ptheta[J+1,],gamma1,gamma2,log=TRUE))
    currobj <- lik_term+prior_term
    return(currobj)
  }
  
  # model constants
  I <- nrow(X)
  J <- ncol(X)
  
  # initializers
  if(!is.null(seed)){set.seed(seed)}
  pi <- rep(1/K,K)
  ptheta <- matrix(c(rbeta(J*K,a,b),rgamma(K,gamma1,gamma2)),byrow=TRUE,nrow=J+1,ncol=K)
  Z <- matrix(NA,nrow=I,ncol=K)
  
  #get current loglikelihood
  curr_obj <- get_obj(ptheta,pi)
  
  currdiff <- Inf
  iter <- 0
  while(currdiff>tol & iter<=maxit){
    iter <- iter + 1
    
    # E-Step: Update Z
    t1 <- Sys.time()
    for(i in 1:I){
      lprobs <- unlist(lapply(1:K,function(k){
        log(pi[k])+dmb(Pi=Pi[i,],X=X[i,],p=ptheta[1:J,k],theta=ptheta[J+1,k],M=M,log=T)
      }))
      Z[i,] <- exp(lprobs-logSumExp(lprobs))
    }
    
    
    # M-Step:
    ## Update (p,theta)
    t2 <- Sys.time()
    for(k in 1:K){
      #print(k)
      zhat_k <- Z[,k]
      # which_vals <- which(zhat_k>.001)
      # if(length(which_vals)==0){which_vals <- 1:I}
      obj <- function(par){
        lik_term <- sum(unlist(lapply(1:I,function(i){
          zhat_k[i]*dmb(Pi=Pi[i,],X=X[i,],p=par[1:J],theta=par[J+1],M=M,log=T)
        })))
        prior_term <- sum(dbeta(par[1:J],a,b,log=TRUE))+dgamma(par[J+1],gamma1,gamma2,log=TRUE)
        return(-lik_term-prior_term)
      }
      ptheta[,k] <- optim(ptheta[,k],obj,method="L-BFGS-B",lower=seq(1e-8,1e-7,length=J+1),
                          upper=c(1-seq(1e-8,1e-7,length=J),qgamma(.9999,gamma1,gamma2)))$par
    }
    
    ## Update pi
    t3 <- Sys.time()
    pi <- (gamma-1+apply(Z,2,sum))/(I - K + gamma*K)
    if(round(sum(pi),5)!=1){print("Error! Something with pi is wrong")}
    
    
    # Iteration Updates
    t4 <- Sys.time()
    new_obj <- get_obj(ptheta,pi)
    currdiff <- new_obj-curr_obj
    curr_obj <- new_obj
    t5 <- Sys.time()
    
    if(verbose){
      print(paste0("End of Iteration ",iter))
      print("Time to update Z, ptheta,pi,objective: ")
      print(round(c(t2-t1,t3-t2,t4-t3,t5-t4),2))
      print(paste0("Difference: ",round(currdiff,3)))
    }
  }
  
  return(list(p=ptheta[1:J,],theta=ptheta[J+1,],pihat=pi,Z=Z,obj=curr_obj))
  
}

# MFM Estimation
mb_mfm <- function(Pi,X,M,lambda,a,b,gamma1,gamma2,xi1,xi2,
                     startK = 1, mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                     max_iters = 100, mh_iters = 5, burn = 0.5,thin = 2,seed = NULL){
  
  print("Initializing Chain")
  if(!is.null(seed)){set.seed(seed)}
  
  ## Determine Constants
  I <- nrow(X)
  J <- ncol(X)
  R <- ncol(Pi)
  
  ## Set Up Data Storage
  which_keep <- seq(from=round(burn*max_iters*mh_iters),to=max_iters * mh_iters,by=thin)
  K_all <- c()
  Kplus_all <- c()
  gamma_all <- c()
  pi_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  p_all <- array(NA,dim=c(J,I,length(which_keep)))
  theta_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  Z_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  accept_p <- c()
  accept_theta <- c()
  accept_gamma <- c()
  
  ## Initialize at Random
  K <- startK
  gamma <- rgamma(1,xi1,xi2)
  pi <- c(rdirichlet(1,rep(gamma,K)))
  ptheta <- matrix(NA,nrow=J+1,ncol=K)
  for(k in 1:K){ptheta[,k] <- c(rbeta(J,a,b),rgamma(1,gamma1,gamma2))}
  Z <- rep(NA,I)
  
  print("Starting Chain")
  iter <- 1
  progress_iters <- round(seq(0,max(which_keep),length=11))[-1]
  while(iter <= max(which_keep)){
    
    ## Step 1: Update Labels
    probs <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      probs[,k] <- exp(log(pi[k])+unlist(lapply(1:I,function(i){
        dmb(Pi=Pi[i,],X=X[i,],p=ptheta[1:J,k],theta=ptheta[J+1,k],M=M,log=T)
      })))
    }
    Z <- apply(probs,1,function(prob){sample.int(K,1,prob=prob)})
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)})) #update Nk, Kplus, and reorder
    Kplus <- sum(Nk>0)
    ptheta <- as.matrix(ptheta[,which(Nk>0)])
    Ztmp <- Z
    for(k in 1:K){Ztmp[which(Z==k)] <- which(c(which(Nk>0),setdiff(1:K,which(Nk>0)))==k)}
    Z <- Ztmp
    rm(Ztmp)
    
    
    ## Step 2: Update non-empty component parameters
    for(mh_iter in 0:(mh_iters-1)){
      for(k in 1:Kplus){
        whichk <- which(Z == k)
        constant1 <- a+apply(matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        constant2 <- b+apply(M-matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        Pi_mat <- matrix(Pi[whichk,],nrow=length(whichk))
        
        for(j in 1:J){ #update each p_jk
          prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)
          while(prop_pjk<=0 | prop_pjk>=1){prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)}
          prop_p <- ptheta[1:J,k]
          prop_p[j] <- prop_pjk
          prop_pi0 <- order(prop_p)
          curr_pi0 <- order(ptheta[1:J,k])
          if(any(prop_pi0!=curr_pi0)){
            logprob_prop <- dmallow(Pi=Pi_mat,pi0=prop_pi0,theta=ptheta[J+1,k],log=T)+
              (constant1[j])*log(prop_pjk)+(constant2[j])*log(1-prop_pjk)
            logprob_curr <- dmallow(Pi=Pi_mat,pi0=curr_pi0,theta=ptheta[J+1,k],log=T)+
              (constant1[j])*log(ptheta[j,k])+(constant2[j])*log(1-ptheta[j,k])
          }else{
            logprob_prop <- (constant1[j])*log(prop_pjk)+(constant2[j])*log(1-prop_pjk)
            logprob_curr <- (constant1[j])*log(ptheta[j,k])+(constant2[j])*log(1-ptheta[j,k])
          }
          
          u <- runif(1)
          if(log(u) <  logprob_prop-logprob_curr){
            accept_p <- c(accept_p,1)
            ptheta[j,k] <- prop_pjk
          }else{accept_p <- c(accept_p,0)}
        }
        #update theta_k
        prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)
        while(prop_thetak<=0){prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)}
        curr_pi0 <- order(ptheta[1:J,k])
        
        
        logprob_prop <- dmallow(Pi=Pi_mat,pi0=curr_pi0,theta=prop_thetak,log=TRUE)+
          (gamma1-1)*log(prop_thetak)-gamma2*prop_thetak
        logprob_curr <- dmallow(Pi=Pi_mat,pi0=curr_pi0,theta=ptheta[J+1,k],log=TRUE)+
          (gamma1-1)*log(ptheta[J+1,k])-gamma2*ptheta[J+1,k]
          
        u <- runif(1)
        if(log(u) <  logprob_prop-logprob_curr){
          accept_theta <- c(accept_theta,1)
          ptheta[J+1,k] <- prop_thetak
        }else{accept_theta <- c(accept_theta,0)}
      }
      if(iter+mh_iter %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        p_all[,1:Kplus,it] <- ptheta[1:J,] #update only the first Kplus parameters,
        theta_all[it,1:Kplus] <- ptheta[J+1,] #(will update the rest later!)
      }
    }
    
    
    ## Step 3: Update K and gamma
    
    #update K
    probs <- unlist(lapply(Kplus:(Kplus+100),function(k){
      dpois(k-1,lambda,log=T)+lfactorial(k)-lfactorial(k-Kplus)+lgamma(gamma*k)-lgamma(I+gamma*k)
    }))
    K <- sample(Kplus:(Kplus+100),1,prob = exp(probs-logSumExp(probs)))
    
    #update gamma
    prop_gamma <- rnorm(1,gamma,mh_gamma)
    while(prop_gamma<=0){prop_gamma <- rnorm(1,gamma,mh_gamma)}
    logprob_prop <- dgamma(prop_gamma,xi1,xi2,log=TRUE)+lgamma(prop_gamma*K)-lgamma(I+prop_gamma*K)+
      sum(lgamma(Nk[1:Kplus]+prop_gamma)-lgamma(prop_gamma))
    logprob_curr <- dgamma(gamma,3,2,log=TRUE)+lgamma(gamma*K)-lgamma(I+gamma*K)+
      sum(lgamma(Nk[1:Kplus]+gamma)-lgamma(gamma))
    u <- runif(1)
    if(log(u) <  logprob_prop-logprob_curr){
      accept_gamma <- c(accept_gamma,1)
      gamma <- prop_gamma
    }else{accept_gamma <- c(accept_gamma,0)}
    
    ## Step 4: Update empty components
    if(K>Kplus){
      for(k in (Kplus+1):K){ptheta <- cbind(ptheta,c(rbeta(J,a,b),rgamma(1,gamma1,gamma2)))}
    }
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)}))
    pi <- c(rdirichlet(1,gamma+Nk))
    
    if(length(pi)!=K | length(Nk)!=K | ncol(ptheta)!=K ){stop("something wrong!")}
    if(Kplus > K){stop("something wrong with Kplus!")}
    
    ## Save Values and Update Counter
    for(mh_iter in 0:(mh_iters-1)){
      if((iter+mh_iter) %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        K_all <- c(K_all,K)
        Kplus_all <- c(Kplus_all,Kplus)
        gamma_all <- c(gamma_all,gamma)
        pi_all[it,1:K] <- pi
        Z_all[it,] <- Z
        if(K>Kplus){
          p_all[,(Kplus+1):K,it] <- ptheta[1:J,(Kplus+1):K] #update only the first Kplus parameters,
          theta_all[it,(Kplus+1):K] <- ptheta[J+1,(Kplus+1):K] #(will update the rest later!)
        }
      }}
    
    iter <- iter + mh_iters
    if(length(progress_iters)>0 & iter >= progress_iters[1]){
      print(paste0((11-length(progress_iters))*10,"% Complete: Iteration ",iter-1," out of ",max(which_keep)))
      print(paste0("Current K = ",K,"; Kplus = ",Kplus))
      progress_iters <- progress_iters[-1]
    }
  }
  print(paste0("Done! Saving ",length(which_keep)," estimate iterations after burning/thinning"))
  return(list(Z=Z_all,p=p_all[,1:max(K_all),],theta=theta_all[,1:max(K_all)],
              pi=pi_all[,1:max(K_all)],K=K_all,Kplus=Kplus_all,gamma=gamma_all,
              accept_p=(cumsum(accept_p)/1:length(accept_p))[seq(round(length(accept_p)*burn),length(accept_p),by=thin)],
              accept_theta=(cumsum(accept_theta)/1:length(accept_theta))[seq(round(length(accept_theta)*burn),length(accept_theta),by=thin)],
              accept_gamma=(cumsum(accept_gamma)/1:length(accept_gamma))[seq(round(length(accept_gamma)*burn),length(accept_gamma),by=thin)],
              max_iters=max_iters,mh_iters=mh_iters,burn=burn,thin=thin,seed=seed))
}

# Finite Mixture (Fixed K) Estimation
mb_fm <- function(Pi,X,M,K,a,b,gamma1,gamma2,xi1,xi2,
                   mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                   max_iters = 100, mh_iters = 5, burn = 0.5,thin = 2,seed = NULL){
  
  print("Initializing Chain")
  if(!is.null(seed)){set.seed(seed)}
  
  ## Determine Constants
  I <- nrow(X)
  J <- ncol(X)
  R <- ncol(Pi)
  
  ## Set Up Data Storage
  which_keep <- seq(from=round(burn*max_iters*mh_iters),to=max_iters * mh_iters,by=thin)
  gamma_all <- c()
  pi_all <- matrix(NA,nrow=length(which_keep),ncol=K)
  p_all <- array(NA,dim=c(J,K,length(which_keep)))
  theta_all <- matrix(NA,nrow=length(which_keep),ncol=K)
  Z_all <- matrix(NA,nrow=length(which_keep),ncol=I)
  accept_p <- c()
  accept_theta <- c()
  accept_gamma <- c()
  
  ## Initialize at Random
  gamma <- rgamma(1,xi1,xi2)
  pi <- c(rdirichlet(1,rep(gamma,K)))
  ptheta <- matrix(NA,nrow=J+1,ncol=K)
  for(k in 1:K){ptheta[,k] <- c(rbeta(J,a,b),rgamma(1,gamma1,gamma2))}
  Z <- rep(NA,I)
  
  print("Starting Chain")
  iter <- 1
  progress_iters <- round(seq(0,max(which_keep),length=11))[-1]
  while(iter <= max(which_keep)){
    
    ## Step 1: Update Labels
    probs <- matrix(NA,nrow=I,ncol=K)
    for(k in 1:K){
      probs[,k] <- exp(log(pi[k])+unlist(lapply(1:I,function(i){
        dmb(Pi=Pi[i,],X=X[i,],p=ptheta[1:J,k],theta=ptheta[J+1,k],M=M,log=T)
      })))
    }
    Z <- apply(probs,1,function(prob){sample.int(K,1,prob=prob)})
    Nk <- unlist(lapply(1:K,function(k){sum(Z==k)})) 
    
    
    ## Step 2: Update component parameters
    for(mh_iter in 0:(mh_iters-1)){
      for(k in 1:K){
        whichk <- which(Z == k)
        if(length(whichk)==0 ){whichk <- 1:I}
        constant1 <- a+apply(matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        constant2 <- b+apply(M-matrix(X[whichk,],nrow=length(whichk)),2,sum,na.rm=T)-1
        Pi_mat <- matrix(Pi[whichk,],nrow=length(whichk))
        
        for(j in 1:J){ #update each p_jk
          prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)
          while(prop_pjk<=0 | prop_pjk>=1){prop_pjk <- rnorm(1,ptheta[j,k],mh_pjk)}
          prop_p <- ptheta[1:J,k]
          prop_p[j] <- prop_pjk
          prop_pi0 <- order(prop_p)
          curr_pi0 <- order(ptheta[1:J,k])
          if(any(prop_pi0!=curr_pi0)){
            logprob_prop <- dmallow(Pi=Pi_mat,pi0=prop_pi0,theta=ptheta[J+1,k],log=T)+
              (constant1[j])*log(prop_pjk)+(constant2[j])*log(1-prop_pjk)
            logprob_curr <- dmallow(Pi=Pi_mat,pi0=curr_pi0,theta=ptheta[J+1,k],log=T)+
              (constant1[j])*log(ptheta[j,k])+(constant2[j])*log(1-ptheta[j,k])
          }else{
            logprob_prop <- (constant1[j])*log(prop_pjk)+(constant2[j])*log(1-prop_pjk)
            logprob_curr <- (constant1[j])*log(ptheta[j,k])+(constant2[j])*log(1-ptheta[j,k])
          }
          
          u <- runif(1)
          if(log(u) <  logprob_prop-logprob_curr){
            accept_p <- c(accept_p,1)
            ptheta[j,k] <- prop_pjk
          }else{accept_p <- c(accept_p,0)}
        }
        #update theta_k
        prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)
        while(prop_thetak<=0){prop_thetak <- rnorm(1,ptheta[J+1,k],mh_thetak)}
        curr_pi0 <- order(ptheta[1:J,k])
        
        
        logprob_prop <- dmallow(Pi=Pi_mat,pi0=curr_pi0,theta=prop_thetak,log=TRUE)+
          (gamma1-1)*log(prop_thetak)-gamma2*prop_thetak
        logprob_curr <- dmallow(Pi=Pi_mat,pi0=curr_pi0,theta=ptheta[J+1,k],log=TRUE)+
          (gamma1-1)*log(ptheta[J+1,k])-gamma2*ptheta[J+1,k]
        
        u <- runif(1)
        if(log(u) <  logprob_prop-logprob_curr){
          accept_theta <- c(accept_theta,1)
          ptheta[J+1,k] <- prop_thetak
        }else{accept_theta <- c(accept_theta,0)}
      }
      if(iter+mh_iter %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        p_all[,,it] <- ptheta[1:J,] 
        theta_all[it,] <- ptheta[J+1,] 
      }
    }
    
    
    ## Step 3: Update gamma
    prop_gamma <- rnorm(1,gamma,mh_gamma)
    while(prop_gamma<=0){prop_gamma <- rnorm(1,gamma,mh_gamma)}
    logprob_prop <- dgamma(prop_gamma,xi1,xi2,log=TRUE)+lgamma(prop_gamma*K)-lgamma(I+prop_gamma*K)+
      sum(lgamma(Nk[1:K]+prop_gamma)-lgamma(prop_gamma))
    logprob_curr <- dgamma(gamma,3,2,log=TRUE)+lgamma(gamma*K)-lgamma(I+gamma*K)+
      sum(lgamma(Nk[1:K]+gamma)-lgamma(gamma))
    u <- runif(1)
    if(log(u) <  logprob_prop-logprob_curr){
      accept_gamma <- c(accept_gamma,1)
      gamma <- prop_gamma
    }else{accept_gamma <- c(accept_gamma,0)}
    
    ## Step 4: Update pi
    pi <- c(rdirichlet(1,gamma+Nk))
    if(length(pi)!=K | length(Nk)!=K | ncol(ptheta)!=K ){stop("something wrong!")}
    
    ## Save Values and Update Counter
    for(mh_iter in 0:(mh_iters-1)){
      if((iter+mh_iter) %in% which_keep){
        it <- which(which_keep==(iter+mh_iter))
        gamma_all <- c(gamma_all,gamma)
        pi_all[it,1:K] <- pi
        Z_all[it,] <- Z
      }}
    
    iter <- iter + mh_iters
    if(length(progress_iters)>0 & iter >= progress_iters[1]){
      print(paste0((11-length(progress_iters))*10,"% Complete: Iteration ",iter-1," out of ",max(which_keep)))
      progress_iters <- progress_iters[-1]
    }
  }
  print(paste0("Done! Saving ",length(which_keep)," estimate iterations after burning/thinning"))
  return(list(Z=Z_all,p=p_all,theta=theta_all,pi=pi_all,gamma=gamma_all,
              accept_p=(cumsum(accept_p)/1:length(accept_p))[seq(round(length(accept_p)*burn),length(accept_p),by=thin)],
              accept_theta=(cumsum(accept_theta)/1:length(accept_theta))[seq(round(length(accept_theta)*burn),length(accept_theta),by=thin)],
              accept_gamma=(cumsum(accept_gamma)/1:length(accept_gamma))[seq(round(length(accept_gamma)*burn),length(accept_gamma),by=thin)],
              max_iters=max_iters,mh_iters=mh_iters,burn=burn,thin=thin,seed=seed))
}


#### MB Examples ####

if(FALSE){ #change to "TRUE" to run
  
  # data generation
  set.seed(1)
  M <- 4
  dat1 <- rmb(I=20,p=runif(5),theta=2,M=M)
  dat2 <- rmb(I=10,p=runif(5),theta=1,M=M)
  X <- rbind(dat1$X,dat2$X)
  Pi <- rbind(dat1$Pi,dat2$Pi)
  rm(dat1,dat2)
  
  # density functions
  dmb(Pi=Pi,X=X,p=runif(ncol(X)),theta=1,M=M,log=TRUE)
  dmallow(Pi=Pi,pi0=sample(1:ncol(X)),theta=1,log=TRUE)
  
  # map
  map <- mb_map(Pi=Pi,X=X,M=M,K=2,gamma=1,a=2,b=2,gamma1=3,gamma2=1,
                tol=.01,maxit=50,verbose=TRUE,seed=1)
 
  # calculation of prior on K+
  pmfstatic2 <- nClusters(Kplus=1:5,N=nrow(X),type="static",gamma=1,maxK=20)
  dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = 1))
  ggplot(data.frame(K=1:length(dens),Mass=dens),aes(K,Mass))+geom_line()+geom_point()+
    ylab("Prior Mass on K+")+xlab("K+")
  
  # mfm
  mfm <- mb_mfm(Pi=Pi,X=X,M=M,lambda=1,a=2,b=2,gamma1=3,gamma2=1,xi1=2,xi2=2,
                startK = 1, mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
                max_iters = 100, mh_iters = 10, burn = 0.5,thin = 2,seed = 1)
  
  # fm
  fm <- mb_fm(Pi=Pi,X=X,M=M,K=2,a=2,b=2,gamma1=3,gamma2=1,xi1=2,xi2=2,
              mh_pjk = 0.01, mh_thetak = 1, mh_gamma = 0.1,
              max_iters = 200, mh_iters = 5, burn = 0.5, thin = 2, seed = 1)
}

