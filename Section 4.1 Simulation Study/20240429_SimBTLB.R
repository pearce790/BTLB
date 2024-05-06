### Initialize: ####
source("20240429_source.R")
parameters <- expand.grid(I = c(25,50,100),K = c(1,2,4),seed = 1:20)
iteration <- as.numeric(commandArgs(trailingOnly=TRUE)[1]) 
I <- parameters[iteration,1]
K <- parameters[iteration,2]
seed <- parameters[iteration,3]
iteration_character <- paste0(paste0(rep(0,3-nchar(iteration)),collapse=""),iteration,sep = "")

#### Draw Parameters ####
set.seed(seed)
J <- 8
R <- J
M <- 9
beta_mu <- c(5,10)
beta_nu <- c(5,1)
mu_0j <- rep(0,J)
sigma2_mu <- rinvgamma(1,beta_mu[1],beta_mu[2])
sigma2_nu <- rinvgamma(1,beta_nu[1],beta_nu[2])
mu <- matrix(rnorm(J*K,mean=mu_0j,sd=sqrt(sigma2_mu)),nrow=J)
nu <- rnorm(n=I,mean=0,sd=sqrt(sigma2_nu))
theta <- rep(20,K)
alpha <- rep(1/K,K)

#### Case 1: True BTL-Binomial ####
set.seed(seed)
Z <- rep(1:K,ceiling(I/K))[1:I]
ratings <- rankings <- matrix(NA,nrow=I,ncol=J)
for(i in 1:I){
  data_i <- rbtlb(I=1,p=expit(mu[,Z[i]]+nu[i]),theta=theta[Z[i]],M=M,R=R)
  ratings[i,] <- data_i$ratings
  rankings[i,] <- data_i$rankings
}
fitted_model <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=8,init=NULL,
                          delta_theta=c(5,0.25),delta_gamma=c(2,4),lambda=1,beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=1000,T2=2,sigma_mu=0.2,sigma_nu=0.2,sigma_theta=3,sigma_gamma=0.5)
fitted_model <- stephens_labelswitching(fitted_model,burn_prop=0.8,thin=1)
fitted_model$class <- NULL
res <- list(model="BTL-Binomial",rater_effects=TRUE,ratings=ratings,rankings=rankings,Z=Z,I=I,K=K,seed=seed,
            sigma2_mu=sigma2_mu,sigma2_nu=sigma2_nu,mu=mu,nu=nu,theta=theta,alpha=alpha,
            results=fitted_model)
save(res,file=paste0("results/res_SimBTLB_20240429_Case1_",iteration_character,".Rdata"))

#### Case 2: Misspecified Ratings ####
set.seed(seed)
Z <- rep(1:K,ceiling(I/K))[1:I]
ratings <- rankings <- matrix(NA,nrow=I,ncol=J)
for(i in 1:I){
  data_i <- rbtlb(I=1,p=expit(mu[,Z[i]]),theta=theta[Z[i]],M=M,R=R)
  ratings[i,] <- data_i$ratings
  rankings[i,] <- data_i$rankings
}
fitted_model <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=8,init=NULL,
                          delta_theta=c(5,0.25),delta_gamma=c(2,4),lambda=1,beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=1000,T2=2,sigma_mu=0.2,sigma_nu=0.2,sigma_theta=3,sigma_gamma=0.5)
fitted_model <- stephens_labelswitching(fitted_model,burn_prop=0.8,thin=1)
fitted_model$class <- NULL
res <- list(model="Misspecified\nRatings",rater_effects=FALSE,ratings=ratings,rankings=rankings,Z=Z,I=I,K=K,seed=seed,
            sigma2_mu=sigma2_mu,sigma2_nu=sigma2_nu,mu=mu,nu=rep(0,I),theta=theta,alpha=alpha,
            results=fitted_model)
save(res,file=paste0("results/res_SimBTLB_20240429_Case2_",iteration_character,".Rdata"))

#### Case 3: Misspecified Rankings ####
set.seed(seed)
Z <- rep(1:K,ceiling(I/K))[1:I]
ratings <- rankings <- matrix(NA,nrow=I,ncol=J)
for(i in 1:I){
  data_i <- rbtlb(I=1,p=expit(mu[,Z[i]]+nu[i]),theta=theta[Z[i]],M=M,R=R)
  ratings[i,] <- data_i$ratings
  rankings[i,] <- rmall(I=1,pi0=order(expit(mu[,Z[i]]+nu[i])),theta=1.5,R=R)
}
fitted_model <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=8,init=NULL,
                          delta_theta=c(5,0.25),delta_gamma=c(2,4),lambda=1,beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=1000,T2=2,sigma_mu=0.2,sigma_nu=0.2,sigma_theta=3,sigma_gamma=0.5)
fitted_model <- stephens_labelswitching(fitted_model,burn_prop=0.8,thin=1)
fitted_model$class <- NULL
res <- list(model="Misspecified\nRankings",rater_effects=TRUE,ratings=ratings,rankings=rankings,Z=Z,I=I,K=K,seed=seed,
            sigma2_mu=sigma2_mu,sigma2_nu=sigma2_nu,mu=mu,nu=nu,theta=theta,alpha=alpha,
            results=fitted_model)
save(res,file=paste0("results/res_SimBTLB_20240429_Case3_",iteration_character,".Rdata"))

#### Case 4: BTL-Binomial with rater effects, but assume none ####
set.seed(seed)
Z <- rep(1:K,ceiling(I/K))[1:I]
ratings <- rankings <- matrix(NA,nrow=I,ncol=J)
for(i in 1:I){
  data_i <- rbtlb(I=1,p=expit(mu[,Z[i]]+nu[i]),theta=theta[Z[i]],M=M,R=R)
  ratings[i,] <- data_i$ratings
  rankings[i,] <- data_i$rankings
}
fitted_model <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=8,init=NULL,
                          delta_theta=c(5,0.25),delta_gamma=c(2,4),lambda=1,beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=1000,T2=2,sigma_mu=0.2,sigma_nu=0.2,sigma_theta=3,sigma_gamma=0.5, rater_effects = FALSE)
fitted_model <- stephens_labelswitching(fitted_model,burn_prop=0.8,thin=1)
fitted_model$class <- NULL
res <- list(model="Case 4",rater_effects=FALSE,ratings=ratings,rankings=rankings,Z=Z,I=I,K=K,seed=seed,
            sigma2_mu=sigma2_mu,sigma2_nu=sigma2_nu,mu=mu,nu=nu,theta=theta,alpha=alpha,
            results=fitted_model)
save(res,file=paste0("results/res_SimBTLB_20240429_Case4_",iteration_character,".Rdata"))

