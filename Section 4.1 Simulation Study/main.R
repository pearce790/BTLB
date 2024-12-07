#### Instructions ####
# This .R script may be run in its entirety to replicate the results presented in 
# Section 4.1 of the manuscript. Note that running this file involves estimating 720 
# unique models to data, which may be time-consuming and slow. Specifically,
# if the script is run without parallelization, it would take approximately 
# 1.5--2 days to run on a standard but high-powered 2023 MacBook Pro. Parallelization 
# along the outer loop in the code that follows can reduce total computation time to 
# approximately 30 minutes. Associated slurm files are included in the "slurm files" folder. 
# 
# Additionally, note that the output files from our 720 simulation scenarios are
# included in the provided "results_PrimarySimulation" folder, which allows the
# user to skip model fitting entirely and replicate the creation of plots by
# skipping directly to the "Figures" section of the .R script, if desired.

#### Load Source Files and Setup Simulation Scenarios ####
source("source.R")
library(RColorBrewer)
iterations <- unlist(lapply(paste0("00",as.character(1:180)),function(val){str_sub(val,-3,-1)}))
parameters <- expand.grid(I = c(25,50,100),K = c(1,2,4),seed = 1:20)

#### Code to Run Simulations ####
for(iteration in 1:nrow(parameters)){
  ### Initialize ###
  I <- parameters[iteration,1]
  K <- parameters[iteration,2]
  seed <- parameters[iteration,3]
  iteration_character <- paste0(paste0(rep(0,3-nchar(iteration)),collapse=""),iteration,sep = "")
  
  ### Draw Parameters ###
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
  
  ### Case 1: True BTL-Binomial ###
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
  save(res,file=paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase1_",iteration_character,".Rdata"))
  
  ### Case 2: Misspecified Ratings ###
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
  save(res,file=paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase2_",iteration_character,".Rdata"))
  
  ### Case 3: Misspecified Rankings ###
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
  save(res,file=paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase3_",iteration_character,".Rdata"))
  
  ### Case 4: BTL-Binomial with rater effects, but assume none ###
  set.seed(seed)
  Z <- rep(1:K,ceiling(I/K))[1:I]
  ratings <- rankings <- matrix(NA,nrow=I,ncol=J)
  for(i in 1:I){
    data_i <- rbtlb(I=1,p=expit(mu[,Z[i]]+nu[i]),theta=theta[Z[i]],M=M,R=R)
    ratings[i,] <- data_i$ratings
    rankings[i,] <- data_i$rankings
  }
  fitted_model <- btlb_mfmm_noratereffects(rankings=rankings,ratings=ratings,M=M,K_start=8,init=NULL,
                            delta_theta=c(5,0.25),delta_gamma=c(2,4),lambda=1,beta_mu=beta_mu,mu_0j=mu_0j,
                            T1=1000,T2=2,sigma_mu=0.2,sigma_theta=3,sigma_gamma=0.5)
  fitted_model <- stephens_labelswitching(fitted_model,burn_prop=0.8,thin=1)
  fitted_model$class <- NULL
  res <- list(model="Case 4",rater_effects=FALSE,ratings=ratings,rankings=rankings,Z=Z,I=I,K=K,seed=seed,
              sigma2_mu=sigma2_mu,sigma2_nu=sigma2_nu,mu=mu,nu=nu,theta=theta,alpha=alpha,
              results=fitted_model)
  save(res,file=paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase4_",iteration_character,".Rdata"))
}

#### Figures ####
## Figure 1 ####
results <- matrix(nrow=0,ncol=9)
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase1_",iter,".Rdata"))
    
    MAE_K <- mean(abs(res$results$K[,1]-res$K))
    MAE_nu <- mean(abs(apply(res$results$nu,2,mean)-res$nu))
    
    reorderings <- permutations(res$K,res$K)
    postmean_mu <- apply(res$results$mu[,1:res$K,,drop=FALSE],c(1,2),function(x){mean(x,na.rm=T)})
    order <- reorderings[which.min(apply(reorderings,1,function(order){median(abs(postmean_mu[,order]-res$mu))})),]
    MAE_mu <- mean(abs(postmean_mu[,order]-res$mu))
    
    postmean_theta <- apply(res$results$theta[,1:res$K,drop=FALSE],2,function(x){mean(x,na.rm=T)})
    MAE_theta <- mean(abs(postmean_theta[order]-res$theta))
    
    results <- rbind(results,c(res$model,res$rater_effects,res$I,res$K,res$seed,MAE_K,MAE_nu,MAE_mu,MAE_theta))
  }, error = function(e) {print(iter)})
}
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase2_",iter,".Rdata"))
    
    MAE_K <- mean(abs(res$results$K[,1]-res$K))
    MAE_nu <- mean(abs(apply(res$results$nu,2,mean)-res$nu))
  
    postmean_mu <- apply(res$results$mu[,1:res$K,,drop=FALSE],c(1,2),function(x){mean(x,na.rm=T)})
    postmean_mu <- postmean_mu[,apply(postmean_mu,2,function(x){!any(is.nan(x))}),drop=FALSE]
    reorderings <- permutations(ncol(postmean_mu),ncol(postmean_mu))
    order <- reorderings[which.min(apply(reorderings,1,function(order){median(abs(postmean_mu[,order]-res$mu))})),]
    MAE_mu <- mean(abs(postmean_mu[,order]-res$mu))
    
    postmean_theta <- apply(res$results$theta[,1:res$K,drop=FALSE],2,function(x){mean(x,na.rm=T)})
    MAE_theta <- mean(abs(postmean_theta[order]-res$theta))
    
    results <- rbind(results,c("Misspecified\nRatings",res$rater_effects,res$I,res$K,res$seed,MAE_K,MAE_nu,MAE_mu,MAE_theta))
  }, error = function(e){print(iter)} )
}
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase3_",iter,".Rdata"))
    
    MAE_K <- mean(abs(res$results$K[,1]-res$K))
    MAE_nu <- mean(abs(apply(res$results$nu,2,mean)-res$nu))
    
    postmean_mu <- apply(res$results$mu[,1:res$K,,drop=FALSE],c(1,2),function(x){mean(x,na.rm=T)})
    postmean_mu <- postmean_mu[,apply(postmean_mu,2,function(x){!any(is.nan(x))}),drop=FALSE]
    reorderings <- permutations(ncol(postmean_mu),ncol(postmean_mu))
    order <- reorderings[which.min(apply(reorderings,1,function(order){median(abs(postmean_mu[,order]-res$mu))})),]
    MAE_mu <- mean(abs(postmean_mu[,order]-res$mu))

    postmean_theta <- apply(res$results$theta[,1:res$K,drop=FALSE],2,function(x){mean(x,na.rm=T)})
    MAE_theta <- mean(abs(postmean_theta[order]-res$theta))
    
    results <- rbind(results,c("Misspecified\nRankings",res$rater_effects,res$I,res$K,res$seed,MAE_K,MAE_nu,MAE_mu,MAE_theta))
  }, error = function(e) {print(iter)})
}
rm(res)
results <- as.data.frame(results)
names(results) <- c("model","raters","I","K","seed","MAE_K","MAE_nu","MAE_mu","MAE_theta")
results$model <- factor(results$model,levels=c("BTL-Binomial","Misspecified\nRatings","Misspecified\nRankings"))
results$raters <- as.logical(results$raters)
results$I <- factor(as.numeric(results$I))
results$K <- factor(as.numeric(results$K),labels=c(expression(K[0]~"= 1"),expression(K[0]~"= 2"),expression(K[0]~"= 4")))
results$seed <- as.numeric(results$seed)
results$MAE_K <- as.numeric(results$MAE_K)
results$MAE_nu <- as.numeric(results$MAE_nu)
results$MAE_mu <- as.numeric(results$MAE_mu)
results$MAE_theta <- as.numeric(results$MAE_theta)
results[results$model=="Misspecified\nRankings","MAE_theta"] <- NA
plot_results <- melt(results,id.vars=c(1,2,3,4,5))
plot_results$variable<- factor(plot_results$variable,levels=c("MAE_K","MAE_nu","MAE_mu","MAE_theta"),
                               labels=c(expression("MAE("~K^~"+"~")"),expression("MAE("~nu~")"),
                                        expression("MAE("~mu~")"),expression("MAE("~theta~")")))
ggplot(plot_results,aes(color=I,y=value,x=model))+geom_boxplot(outlier.size = .5,outlier.alpha=0.5)+
  facet_grid(variable~K,scales="free_y",labeller = label_parsed)+
  scale_color_manual(values=c("skyblue","#3182BD","#08519C"))+
  theme_bw()+theme(legend.position = "bottom",panel.grid.major.x = element_blank(),
                   panel.grid.minor = element_blank())+
  labs(color="Number of Judges, I",y=element_blank(),x=element_blank())
ggsave("ResultsPlots/Fig1.png",plot=last_plot(),width=10,height=4.5,units="in") 

## Figure 2 ####
results <- matrix(nrow=0,ncol=4)
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/results_PrimarySimulation/resCase4_",iter,".Rdata"))
    MAE_K <- mean(abs(res$results$K[,1]-res$K))
    results <- rbind(results,c(res$I,res$K,res$seed,MAE_K))
  }, error = function(e) print(iter))
}
results <- as.data.frame(results)
names(results) <- c("I","K","seed","MAE_K")
results$I <- factor(as.numeric(results$I))
results$K <- factor(as.numeric(results$K),labels=c(expression(K[0]~"= 1"),expression(K[0]~"= 2"),expression(K[0]~"= 4")))
results$seed <- as.numeric(results$seed)
results$MAE_K <- as.numeric(results$MAE_K)
ggplot(results,aes(x=I,y=MAE_K))+geom_boxplot(outlier.size = .5,outlier.alpha=0.5)+
  facet_grid(~K,labeller = label_parsed)+
  theme_bw()+theme(legend.position = "bottom")+scale_y_continuous(breaks=seq(0,100,by=1))+
  labs(x="Number of Judges, I",y=expression("MAE("~K^~"+"~")"))
ggsave("ResultsPlots/Fig2.png",plot=last_plot(),width=9,height=2,units="in") 

