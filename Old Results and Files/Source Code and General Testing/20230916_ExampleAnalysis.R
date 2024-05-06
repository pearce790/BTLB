source("~/BTLBinomial/20230826_source.R")

## Set Constants and Hyperparameters ####

I <- 100
J <- 8
M <- 25
R <- rep(4,I)
assignments <- matrix(FALSE,nrow=I,ncol=J)
set.seed(0)
for(i in 1:I){assignments[i,sample(1:8,6)] <- TRUE}

delta_gamma <- c(2,2)
beta_mu <- c(5,10)
beta_nu <- c(10,1)
lambda <- 2
mu_0j <- seq(-2,2,length=J)
delta_theta <- c(10,1)

## Sample Data for Testing New Formulation ####

set.seed(1)
gamma <- rgamma(1,delta_gamma[1],delta_gamma[2])
sigma2_mu <- rinvgamma(1,beta_mu[1],beta_mu[2])
sigma2_nu <- rinvgamma(1,beta_nu[1],beta_nu[2])
K <- rpois(1,lambda)+1
alpha <- c(rdirichlet(1,rep(gamma,K)))
Z <- replicate(I,sample(1:K,1,prob=alpha))
mu <- matrix(rnorm(n=J*K,mean=mu_0j,sd=sqrt(sigma2_mu)),nrow=J,ncol=K)
nu <- rnorm(I,mean=0,sd=sqrt(sigma2_nu))
theta <- rgamma(K,delta_theta[1],delta_theta[2])

rankings <- ratings <- matrix(NA,nrow=I,ncol=J)
for(i in 1:I){
  tmp <- rbtlb(I=1,p=expit(mu[,Z[i]] + nu[i]),theta=theta[Z[i]],M=M,R=R[i],assignments=matrix(assignments[i,],nrow=1))
  rankings[i,] <- tmp$rankings
  ratings[i,] <- tmp$ratings
}
attr(rankings,"assignments") <- assignments
true_values <- list(gamma=gamma,sigma2_mu=sigma2_mu,sigma2_nu=sigma2_nu,K=K,alpha=alpha,Z=Z,mu=mu,nu=nu,theta=theta)
rm(tmp,assignments,i,R,gamma,sigma2_mu,sigma2_nu,K,alpha,Z,mu,nu,theta,I,J)


## MAP Estimation (Fixed K) ####
set.seed(1)
map_hat <- btlb_map(rankings,ratings,M,K=3,tol=1,verbose=TRUE,
                    delta_gamma = delta_gamma,beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,delta_theta = delta_theta)
table(apply(map_hat$z_hat,1,which.max),true_values$Z)

## Bayesian Estimation under Telescoping Sampler (Variable K) ####

## Run short mcmc chains to ensure reasonable Metropolis-Hastings parameters, i.e., tune the acceptance probabilities
set.seed(2)
mfmm_test <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=map_hat$K,init=map_hat,
                       delta_theta = delta_theta, delta_gamma=delta_gamma,lambda=lambda,
                       T1=100,T2=2,sigma_mu=0.5,sigma_nu=0.5,sigma_theta=2.5,sigma_gamma=2.5)
ggplot(melt(apply(mfmm_test$accept,2,cummean)),aes(Var1,value,group=Var2,color=Var2))+
  theme_bw()+scale_color_brewer(type="seq",palette="Dark2")+
  scale_y_continuous(limits=c(0,1))+geom_line()+
  labs(color="Variable",x="Iteration",y="Cumulative Acceptance Probability",title="Metropolis-Hastings Acceptance Probabilities")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())

## Run longer chain to develop proper posterior distribution
set.seed(2)
mfmm_mcmc <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=map_hat$K,init=map_hat,
                       delta_theta = delta_theta, delta_gamma=delta_gamma,lambda=lambda,
                       T1=2000,T2=2,sigma_mu=0.5,sigma_nu=0.5,sigma_theta=2.5,sigma_gamma=2.5)

## Visualize Results ####
# diagnostics
g1<-ggplot(melt(apply(mfmm_mcmc$accept,2,cummean)),aes(Var1,value,group=Var2,color=Var2))+
  theme_bw()+scale_color_brewer(type="seq",palette="Dark2")+
  scale_y_continuous(limits=c(0,1))+geom_line()+
  labs(color="Variable",x="Primary Iteration",y="Cumulative Acceptance Probability",title="Metropolis-Hastings Acceptance Probabilities")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())
g2<-ggplot(melt(mfmm_mcmc$timing),aes(Var1,value,group=Var2,color=Var2))+
  theme_bw()+scale_color_brewer(type="seq",palette="Dark2")+geom_line()+
  scale_y_log10(breaks = trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x)))+
  labs(color="Algorithm Step",x="Primary Iteration",y="Time (seconds)",title="Algorithm Speed by Iteration")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())
grid.arrange(g1,g2,nrow=1)

# examine important trace plots to assess convergence
ggplot(melt(mfmm_mcmc$K),aes(Var1,value,color=factor(Var2,levels=c("K","Kplus"),labels=c("K","K+"))))+
  geom_line()+theme_bw()+
  scale_y_continuous(limits=c(1,max(mfmm_mcmc$K)),breaks=1:max(mfmm_mcmc$K))+
  scale_color_brewer(type="seq",palette="Dark2")+
  labs(x="Primary Iteration",y="Number of Classes",color="",title="Trace Plot: Number of Classes (K) and Clusters (K+)")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
ggplot(melt(mfmm_mcmc$alpha),aes(Var1,value,color=factor(Var2)))+
  geom_line()+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  labs(x="Primary Iteration",y="Class Weights",color="Class Labels",title="Trace Plot: Class Weights")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  guides(color = guide_legend(nrow = 1))
ggplot(melt(mfmm_mcmc$sigma2),aes(Var1,value,color=factor(Var2)))+
  geom_line()+theme_bw()+scale_color_brewer(type="seq",palette="Dark2")+
  labs(x="Primary Iteration",y="Variance Hyperparameter",color="",title="Trace Plot: Variance Hyperparameters")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
ggplot(melt(mfmm_mcmc$gamma),aes(Var1,value))+
  geom_line()+theme_bw()+
  labs(x="Primary Iteration",y=expression(gamma),color="",title="Trace Plot: Class Concentration")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# select burn-in, thinning, and number of clusters to display
K_display <- 3
burn_in_proportion <- 0.25
thin_by <- 1

# examine posterior summaries
number_iterations <- nrow(mfmm_mcmc$gamma)
keep_iters <- seq(burn_in_proportion*number_iterations+1,number_iterations,by=1)
keep_iters_K_display <- keep_iters[mfmm_mcmc$K[keep_iters,2] == K_display]


ggplot(melt(mfmm_mcmc$K[keep_iters,]),aes(x=1,fill=factor(value)))+geom_bar()+
  facet_grid(~factor(Var2,levels=c("K","Kplus"),labels=c("Classes","Clusters")))+
  theme_bw()+scale_fill_viridis_d()+
  labs(fill=element_blank(),y="Proportion",x=element_blank(),title="Posterior Distribution: Classes and Clusters")+
  scale_y_continuous(breaks=seq(0,length(keep_iters),length=5),labels=seq(0,1,length=5))+
  theme(legend.position = "bottom",panel.grid = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank())
gres_gamma <- ggplot(data.frame(gamma=mfmm_mcmc$gamma[keep_iters,,drop=FALSE]),aes(x=1,y=gamma))+
  geom_violin()+theme_bw()+
  labs(x="",y=expression(gamma),title="Class Concentration Parameter")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank())
gres_mu <- ggplot(data.frame(sigma2_mu=mfmm_mcmc$sigma2[keep_iters,1]),aes(x=1,y=sigma2_mu))+
  geom_violin()+theme_bw()+
  labs(x="",y=expression(sigma^2[mu]),title="Object Quality Variance")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank())
gres_nu <- ggplot(data.frame(sigma2_mu=mfmm_mcmc$sigma2[keep_iters,2]),aes(x=1,y=sigma2_mu))+
  geom_violin()+theme_bw()+
  labs(x="",y=expression(sigma^2[nu]),title="Judge Leniency Variance")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank())
grid.arrange(gres_gamma,gres_mu,gres_nu,nrow=1)

ggplot(melt(mfmm_mcmc$alpha[keep_iters_K_display,1:K_display]),aes(x=factor(Var2),y=value))+
  geom_violin()+theme_bw()+
  labs(x="Class",y="Density",title=expression("Posterior Distribution: "~alpha))+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())

ggplot(melt(mfmm_mcmc$Z[keep_iters_K_display,order(apply(mfmm_mcmc$Z[keep_iters_K_display,],2,mean))]),
       aes(x=factor(Var2),fill=factor(value),color=factor(value)))+
  geom_bar()+theme_minimal()+
  labs(x="Judges",y="Proportion",fill="Class",color="Class",title="Posterior Class Membership Probabilities",
       subtitle="Judges ordered by modal class membership")+
  scale_y_continuous(breaks=seq(0,length(keep_iters_K_display),length=5),labels=seq(0,1,length=5))+
  scale_color_brewer(type="seq",palette="Dark2")+
  scale_fill_brewer(type="seq",palette="Dark2")+
  theme(panel.grid = element_blank(),legend.position = "bottom",axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks.y=element_line())

nu_summary <- data.frame(t(apply(mfmm_mcmc$nu,2,function(x){quantile(x,c(0.25,0.5,0.75))})))
nu_summary <- nu_summary %>% arrange(X50.); nu_summary$Judge <- 1:100
names(nu_summary) <- c("Lower","Median","Upper","Judge")
ggplot(nu_summary,aes(x=Judge,y=Median,ymin=Lower,ymax=Upper))+
  geom_point(alpha=.5)+geom_errorbar(alpha=.5)+theme_bw()+
  labs(x="Judges",y="Leniency (low) vs. Harshness (high)",title=expression("50% Credible Intervals: "~nu),
       subtitle="Judges ordered by posterior median")+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        axis.text.x = element_blank())

ggplot(melt(mfmm_mcmc$mu[,1:K_display,keep_iters]),
       aes(x=factor(Var1),y=value,color=factor(Var2),fill=factor(Var2)))+
  geom_violin(alpha=.5)+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  scale_fill_brewer(type="seq",palette="Dark2")+
  geom_vline(xintercept=seq(1.5,nrow(mfmm_mcmc$mu)-.5),color="gray90")+
  labs(x="Object",y="Posterior",fill="Cluster",color="Cluster",
       title=expression("Posterior Distribution: "~mu))+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom")+
  geom_point(data=melt(true_values$mu[,c(3,1,2)]),aes(x=factor(Var1),y=value,color=factor(Var2)),
             size=2)



ggplot(melt(mfmm_mcmc$theta[keep_iters,1:K_display]),
       aes(x=factor(Var2),y=value,color=factor(Var2),fill=factor(Var2)))+
  geom_violin()+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  scale_fill_brewer(type="seq",palette="Dark2")+
  labs(x="Class",y="Consensus Scale",fill="Cluster",color="Cluster",
       title=expression("Posterior Distribution: "~theta))+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom")

true_values$alpha[c(2,3,1)]
true_values$theta[c(3,1,2)]

