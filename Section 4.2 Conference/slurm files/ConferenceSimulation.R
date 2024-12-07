### Initial Setup ####
source("source.R")
iter <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

# creating matrix to store results for each simulation
results <- matrix(NA,nrow=0,ncol=255)

# setting constants/parameters for all combinations of M, R, theta
I <- J <- 50
delta_gamma <- c(2,2)
beta_mu <- c(5,10)
beta_nu <- c(10,1)
mu_0j <- rep(0,J)
delta_theta <- c(5,0.25)
K <- 1
tol <- 1
verbose <- FALSE

### Iterating over M, R, theta ####
for(M in c(4,9)){
  for(R in c(4,8,12,24)){
    for(theta in c(1,10,20,40)){
      
      #initialize
      set.seed(M*R*theta*iter)
      print(paste0("M=",M,", R=",R,", theta=",theta,", iter=",iter))
      
      #make assignments at random
      continue <- FALSE
      while(continue==FALSE){
        assignments <- matrix(FALSE,nrow=50,ncol=50)
        tryCatch({
          for(i in 1:I){
            prob <- R-apply(assignments,2,sum)[which(apply(assignments,2,sum)<R)]
            assignments[i,sample(which(apply(assignments,2,sum)<R),R,prob = prob)] <- TRUE
          }
        },error=function(e){})
        if(all(apply(assignments,2,sum)==R)){continue <- TRUE}
      }
      
      # draw parameters
      sigma2_mu <- 1.5^2
      sigma2_nu <- 0.15^2
      mu <- rnorm(J,mean=0,sd=sqrt(sigma2_mu))
      nu <- rnorm(n=I,mean=0,sd=sqrt(sigma2_nu))
      
      # draw data
      ratings <- rankings <- matrix(NA,nrow=I,ncol=J)
      for(i in 1:I){
        data_i <- rbtlb(I=1,p=expit(mu+nu[i]),theta=theta,M=M,R=4,assignments=assignments[i,,drop=FALSE])
        ratings[i,] <- data_i$ratings
        rankings[i,] <- data_i$rankings
      }
      attr(rankings,"assignments") <- assignments
      
      # fit model
      tryCatch({
        map_est <- btlb_map(rankings,ratings,M,K,tol,verbose,delta_gamma,beta_mu,beta_nu,mu_0j,delta_theta)
        results <- rbind(results,c(M,R,iter,mu,c(map_est$mu),
                                   theta,map_est$theta,
                                   apply(ratings,2,function(x){mean(x,na.rm=TRUE)}),
                                   nu,c(map_est$nu)))
      },error=function(e){
        tryCatch({
          map_est <- btlb_map(rankings,ratings,M,K,tol,verbose,delta_gamma,beta_mu,beta_nu,mu_0j,delta_theta)
          results <- rbind(results,c(M,R,iter,mu,c(map_est$mu),
                                     theta,map_est$theta,
                                     apply(ratings,2,function(x){mean(x,na.rm=TRUE)}),
                                     nu,c(map_est$nu)))
        },error=function(e){
          print("Fail MAP Search for M,R,theta,iter:")
          print(c(M,R,theta,iter))
        })
      })
    }
  }
}

### Saving Results
names(results) <- c("M","R","iter",paste0("mu",1:50),paste0("mu_map",1:50),
                    "theta","theta_map",paste0("Xbar",1:50),paste0("nu",1:50),paste0("nu_map",1:50))
save(results,file=paste0("res_",iter,".Rdata"))
