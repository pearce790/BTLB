## Source Data/Code ####
source("20240429_source.R")
load("Section 4.3 Panel Review/Data/AIBS_Clean.RData")


## Set Constants and Hyperparameters ####
I <- nrow(rankings)
J <- ncol(rankings)

delta_gamma <- c(2,4)
beta_mu <- c(5,10)
beta_nu <- c(5,1)
lambda <- 2
mu_0j <- rep(0,J)
delta_theta <- c(10,0.5)

### Priors
library(fipp)
set.seed(1)
pmfstatic2 <- nClusters(Kplus=1:8,N=I,type="static",gamma=0.25,maxK=20)
dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = lambda))
prior_dens <- data.frame(K=1:length(dens),gamma=0.25,dens=dens)
for(gamma in c(0.5,1,2)){
  pmfstatic2 <- nClusters(Kplus=1:8,N=I,type="static",gamma=gamma,maxK=20)
  dens <- pmfstatic2(priorK = dpois, priorKparams = list(lambda = lambda))
  prior_dens <- rbind(prior_dens,data.frame(K=1:length(dens),gamma=gamma,dens=dens))
}
prior_plot <- grid.arrange(
  ggplot(data.frame(gamma=rgamma(100000,delta_gamma[1],delta_gamma[2])),aes(x=gamma))+
    geom_histogram(aes(y=..density..),binwidth = 0.05)+scale_x_continuous(limits=c(0,3))+
    theme_bw()+labs(x=expression(gamma),y="Density",title=expression("Prior on "~gamma)),
  ggplot(data.frame(K=rpois(100000,lambda)+1),aes(x=K))+
    geom_histogram(aes(y=..density..),binwidth=0.25)+scale_x_continuous(breaks=1:8,limits=c(.5,8.5))+
    theme_bw()+labs(x="K",y="Density",title="Prior on K"),
  ggplot(prior_dens,aes(x=K,y=dens,group=factor(gamma),color=factor(gamma)))+
    geom_line()+geom_point()+scale_x_continuous(breaks=c(1:9))+
    theme_bw()+labs(x="K+",y="Mass",color=expression(gamma),title=expression("Prior on K+|"~gamma))+
    theme(legend.position = c(.8,.6)),
  nrow=1)
ggsave("ResultsPlots/AIBS_prior.png",prior_plot,width=10,height=3,units="in") 


## MAP Estimation (to obtain initializers) ####
set.seed(1)
map_res <- replicate(10,btlb_map(rankings,ratings,M,K=2,tol=1,verbose=FALSE,
                                delta_gamma = delta_gamma,beta_mu=beta_mu,beta_nu=beta_nu,
                                mu_0j=mu_0j,delta_theta = delta_theta))
map_hat_K2 <- map_res[,which.max(map_res["loglikelihood",])]
map_res <- replicate(10,btlb_map(rankings,ratings,M,K=3,tol=1,verbose=FALSE,
                                 delta_gamma = delta_gamma,beta_mu=beta_mu,beta_nu=beta_nu,
                                 mu_0j=mu_0j,delta_theta = delta_theta))
map_hat_K3 <- map_res[,which.max(map_res["loglikelihood",])]


## MFMM Estimation starting at K=2, K=3, and K=17 ####
set.seed(2)
mfmm_mcmc_K2 <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=map_hat_K2$K,init=map_hat_K2,
                          delta_theta=delta_theta,delta_gamma=delta_gamma,lambda=lambda,
                          beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=100000,T2=4,sigma_mu=0.5,sigma_nu=0.2,sigma_theta=2,sigma_gamma=1)
mfmm_mcmc_K3 <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=map_hat_K3$K,init=map_hat_K3,
                          delta_theta=delta_theta,delta_gamma=delta_gamma,lambda=lambda,
                          beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=100000,T2=4,sigma_mu=0.5,sigma_nu=0.2,sigma_theta=2,sigma_gamma=1)
mfmm_mcmc_K17 <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=I,init=NULL,
                       delta_theta=delta_theta,delta_gamma=delta_gamma,lambda=lambda,
                       beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                       T1=100000,T2=4,sigma_mu=0.5,sigma_nu=0.2,sigma_theta=2,sigma_gamma=1)

save.image("Section 4.3 Panel Review/20240203_AIBSAnalysis.RData")
# load("Section 4.3 Panel Review/20240203_AIBSAnalysis.RData")

## Combine Chains ####
n_iters <- nrow(mfmm_mcmc_K2$gamma)
burn_prop <- 0.5
keep_iters <- seq(n_iters*burn_prop+1,n_iters,by=10)

mfmm_mcmc_K3$Z[keep_iters,][c(which(mfmm_mcmc_K3$Z[keep_iters,]==1),which(mfmm_mcmc_K3$Z[keep_iters,]==2),
                              which(mfmm_mcmc_K3$Z[keep_iters,]==3),which(mfmm_mcmc_K3$Z[keep_iters,]==4))] <-
  c(rep(2,sum(mfmm_mcmc_K3$Z[keep_iters,]==1)),rep(1,sum(mfmm_mcmc_K3$Z[keep_iters,]==2)),
    rep(3,sum(mfmm_mcmc_K3$Z[keep_iters,]==3)),rep(4,sum(mfmm_mcmc_K3$Z[keep_iters,]==4)))
mfmm_mcmc <- list(gamma=matrix(c(mfmm_mcmc_K2$gamma[keep_iters,],mfmm_mcmc_K3$gamma[keep_iters,],mfmm_mcmc_K17$gamma[keep_iters,])),
                  sigma2=rbind(mfmm_mcmc_K2$sigma2[keep_iters,],mfmm_mcmc_K3$sigma2[keep_iters,],mfmm_mcmc_K17$sigma2[keep_iters,]),
                  alpha=rbind(mfmm_mcmc_K2$alpha[keep_iters,c(2,1,3:11)],mfmm_mcmc_K3$alpha[keep_iters,1:11],mfmm_mcmc_K17$alpha[keep_iters,c(2,1,3:11)]),
                  mu=abind(mfmm_mcmc_K2$mu[,c(2,1,3:11),keep_iters],mfmm_mcmc_K3$mu[,1:11,keep_iters],mfmm_mcmc_K17$mu[,c(2,1,3:11),keep_iters]),
                  nu=rbind(mfmm_mcmc_K2$nu[keep_iters,],mfmm_mcmc_K3$nu[keep_iters,],mfmm_mcmc_K17$nu[keep_iters,]),
                  theta=rbind(mfmm_mcmc_K2$theta[keep_iters,c(2,1,3:11)],mfmm_mcmc_K3$theta[keep_iters,1:11],mfmm_mcmc_K17$theta[keep_iters,c(2,1,3:11)]),
                  Z=rbind(mfmm_mcmc_K2$Z[keep_iters,],mfmm_mcmc_K3$Z[keep_iters,],mfmm_mcmc_K17$Z[keep_iters,]),
                  K=rbind(mfmm_mcmc_K2$K[keep_iters,],mfmm_mcmc_K3$K[keep_iters,],mfmm_mcmc_K17$K[keep_iters,]),
                  loglik=rbind(mfmm_mcmc_K2$loglik_iters[keep_iters,],mfmm_mcmc_K3$loglik_iters[keep_iters,],mfmm_mcmc_K17$loglik_iters[keep_iters,]),
                  accept=rbind(mfmm_mcmc_K2$accept[keep_iters,],mfmm_mcmc_K3$accept[keep_iters,],mfmm_mcmc_K17$accept[keep_iters,]),
                  timing=rbind(mfmm_mcmc_K2$timing[keep_iters,],mfmm_mcmc_K3$timing[keep_iters,],mfmm_mcmc_K17$timing[keep_iters,]))



## MH Diagnostics and Timings ####
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
AIBS_alg <- grid.arrange(g1,g2,nrow=1)
ggsave("ResultsPlots/AIBS_time.png",AIBS_alg,width=10,height=3,units="in") 

## Posterior Distributions ####
plotK<-ggplot(data.frame(clusters=mfmm_mcmc$K[,1]),aes(x=clusters))+geom_bar()+
  theme_bw()+labs(x="K+",y="Probability")+
  #labs(fill=element_blank(),y="Proportion",x=element_blank(),title="Posterior Distribution: Classes and Clusters")+
  scale_y_continuous(breaks=seq(0,nrow(mfmm_mcmc$K),length=5),labels=seq(0,1,length=5))

plotZ<-ggplot(melt(mfmm_mcmc$Z),
              aes(x=factor(Var2),
                  fill=factor(value,levels=c(2,1,3,4),labels=c(1:4)),
                  color=factor(value,levels=c(2,1,3,4),labels=c(1:4))))+
  geom_bar()+theme_minimal()+
  labs(x="Judges",y="Proportion",fill="Class",color="Class")+
  scale_y_continuous(breaks=seq(0,nrow(mfmm_mcmc$Z),length=5),labels=seq(0,1,length=5))+
  scale_color_brewer(type="seq",palette="Dark2")+
  scale_fill_brewer(type="seq",palette="Dark2")+
  theme(panel.grid = element_blank(),legend.position = "bottom",
        axis.ticks.x = element_blank(),axis.ticks.y=element_line())

class <- data.frame(Var2=1:I,pred_class=round(apply(mfmm_mcmc$Z,2,median)))
nu_long <- melt(mfmm_mcmc$nu) %>% left_join(class,by="Var2")
plot_nu<-ggplot(nu_long,aes(y=factor(Var2,levels=order(apply(mfmm_mcmc$nu,2,median))),x=value,
                            fill=factor(pred_class,levels=c(2,1),labels=c(1,2))))+
  geom_density_ridges(bandwidth=.04,alpha=0.75)+theme_minimal()+
  scale_fill_brewer(type="seq",palette="Dark2")+
  scale_x_continuous(limits=c(-.75,0.75))+
  labs(y="Judge",x="Leniency (low) vs. Harshness (high)",fill="Class")+
  theme(panel.grid = element_blank(),legend.position = "right")
ggsave("ResultsPlots/resAIBS_2.png",plot_nu,width=10,height=4,units="in")

plot_mu<-ggplot(melt(mfmm_mcmc$mu[,1:2,]),
                aes(x=factor(Var1,levels=order(apply(mfmm_mcmc$mu[,1,],1,median))),
                    y=value,color=factor(Var2),fill=factor(Var2)))+
  geom_violin()+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  scale_fill_brewer(type="seq",palette="Dark2")+
  geom_vline(xintercept=seq(1.5,nrow(mfmm_mcmc$mu)-.5),color="gray90")+
  labs(x="Object",y=expression(mu[j]),fill="Cluster",color="Cluster")+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom")

plot_theta<-ggplot(melt(mfmm_mcmc$theta[,1:2]),
                   aes(x=factor(Var2),y=value,color=factor(Var2),fill=factor(Var2)))+
  geom_violin()+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  scale_fill_brewer(type="seq",palette="Dark2")+
  labs(x="Class",y=expression(theta),fill="Class",color="Class")+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom")

plot1 <- grid.arrange(
  grid.arrange(plotK,plotZ+theme(legend.position = "right"),widths=c(0.2,0.8)),
  grid.arrange(plot_mu+theme(legend.position = "none"),plot_theta+theme(legend.position = "right"),
               widths=c(0.8,0.2))
)
ggsave("ResultsPlots/resAIBS_1.png",plot1,width=12,height=5,units="in")

gres_gamma <- ggplot(data.frame(gamma=mfmm_mcmc$gamma),aes(x=1,y=gamma))+
  geom_violin()+theme_bw()+
  labs(x="",title=expression(gamma),y="Posterior")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank())
gres_mu <- ggplot(data.frame(sigma2_mu=mfmm_mcmc$sigma2[,1]),aes(x=1,y=sigma2_mu))+
  geom_violin()+theme_bw()+scale_y_continuous(limits=c(0,max(mfmm_mcmc$sigma2[,1])))+
  labs(x="",title=expression(sigma[mu]^2),y=element_blank())+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank())
gres_nu <- ggplot(data.frame(sigma2_mu=mfmm_mcmc$sigma2[,2]),aes(x=1,y=sigma2_mu))+
  geom_violin()+theme_bw()+scale_y_continuous(limits=c(0,max(mfmm_mcmc$sigma2[,1])))+
  labs(x="",title=expression(sigma[nu]^2),y=element_blank())+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank())
gres_alpha <- ggplot(melt(mfmm_mcmc$alpha[,1:4]),aes(x=factor(Var2),y=value))+
  geom_boxplot(outlier.alpha=0.01)+theme_bw()+
  labs(x="Class",y=element_blank(),title=expression(alpha))+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())
AIBS_param2<-grid.arrange(gres_alpha,gres_gamma,gres_mu,gres_nu,nrow=1,widths=c(.4, 0.2,0.2,0.2))
ggsave("ResultsPlots/resAIBS_3.png",AIBS_param2,width=10,height=3,units="in") 

## Trace Plots ####
melt_K <- melt(mfmm_mcmc$K)
melt_K$Var1 <- rep(1:(nrow(melt_K)/2),2)
traceK<-ggplot(melt_K,aes(Var1,value,color=factor(Var2,levels=c("K","Kplus"),labels=c("K","K+"))))+
  geom_line()+theme_bw()+
  scale_y_continuous(limits=c(1,max(mfmm_mcmc$K)),breaks=1:max(mfmm_mcmc$K))+
  scale_color_brewer(type="seq",palette="Dark2")+
  labs(x="Iteration",y="Value",color="",title="Trace Plot: K, K+")+
  theme(legend.position = "right",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
tracealpha<-ggplot(melt(mfmm_mcmc$alpha[,1:5]),aes(Var1,value,color=factor(Var2)))+
  geom_line()+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  labs(x="Iteration",y="Value",color="k",title=expression("Trace Plot:"~alpha))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
melt_sigma <- melt(mfmm_mcmc$sigma2)
melt_sigma$Var1 <- rep(1:(nrow(melt_sigma)/2),2)
tracesigma<-ggplot(melt_sigma,aes(Var1,value,color=factor(Var2)))+
  geom_line()+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2",labels=c(expression(sigma[mu]^2),expression(sigma[nu]^2)))+
  labs(x="Iteration",y="Value",color="",title="Trace Plot: Variance Hyperparameters")+
  theme(legend.position = "right",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
tracegamma<-ggplot(melt(mfmm_mcmc$gamma),aes(Var1,value))+
  geom_line()+theme_bw()+
  labs(x="Iteration",y="Value",color="",title=expression("Trace Plot:"~gamma))+
  theme(legend.position = "right",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
AIBS_trace1<-grid.arrange(tracegamma,tracesigma,tracealpha,traceK)
ggsave("ResultsPlots/AIBS_trace1.png",AIBS_trace1,width=10,height=5,units="in") 


AIBS_trace2<-ggplot(melt(mfmm_mcmc$mu[,1:2,seq(1,15000,by=10)]),aes(x=Var3,color=factor(Var2),y=value))+
  geom_line(alpha=0.8)+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  facet_wrap(~Var1)+labs(x="Iteration",y="Value",title=expression("Trace Plot: "~mu),color="Class")
AIBS_trace3<-ggplot(melt(mfmm_mcmc$theta[,1:2]),aes(x=Var1,color=factor(Var2),y=value))+
  geom_line(alpha=c(rep(1,15000),rep(.5,15000)))+theme_bw()+
  scale_color_brewer(type="seq",palette="Dark2")+
  labs(x="Iteration",y="Value",title=expression("Trace Plot: "~theta),color="Class")
AIBS_trace4<-ggplot(melt(mfmm_mcmc$nu[seq(1,15000,by=10),]),aes(x=Var1,y=value))+
  geom_line(alpha=0.6)+theme_bw()+
  facet_wrap(~Var2,nrow=3)+labs(x="Iteration",y="Value",title=expression("Trace Plot: "~nu))
ggsave("ResultsPlots/AIBS_trace2.png",AIBS_trace2,width=10,height=5,units="in") 
ggsave("ResultsPlots/AIBS_trace3.png",AIBS_trace3,width=10,height=5,units="in") 
ggsave("ResultsPlots/AIBS_trace4.png",AIBS_trace4,width=10,height=5,units="in") 

## Model Assessment Statistics ####
iters <- seq(1,1500,by=5)
postpred_ratings_mean <- matrix(NA,nrow=length(iters),ncol=J)
for(j in 1:J){
  ratings_j<-c()
  for(iter in iters){
    mus <- mfmm_mcmc$mu[j,3-mfmm_mcmc$Z[iter,],iter] #correcting for label switching!
    nus <- mfmm_mcmc$nu[iter,]
    ratings_j <- c(ratings_j,mean(rbinom(I,M,expit(mus+nus))))
  }
  postpred_ratings_mean[,j] <- ratings_j
}
postpred_ratings_mean <- as.data.frame(postpred_ratings_mean)
names(postpred_ratings_mean) <- 1:J
postpred_ratings_var <- matrix(NA,nrow=length(iters),ncol=J)
for(j in 1:J){
  ratings_j<-c()
  for(iter in iters){
    mus <- mfmm_mcmc$mu[j,3-mfmm_mcmc$Z[iter,],iter] #correcting for label switching!
    nus <- mfmm_mcmc$nu[iter,]
    ratings_j <- c(ratings_j,var(rbinom(I,M,expit(mus+nus))))
  }
  postpred_ratings_var[,j] <- ratings_j
}
postpred_ratings_var <- as.data.frame(postpred_ratings_var)
names(postpred_ratings_var) <- 1:J
results_rankings <- matrix(NA,nrow=0,ncol=J)
for(iter in iters){
  for(i in 1:I){
    class <- 3-mfmm_mcmc$Z[iter,i]
    mu <- mfmm_mcmc$mu[,class,iter]
    nu <- mfmm_mcmc$nu[iter,i]
    theta <- mfmm_mcmc$theta[iter,class]
    tryCatch({
      results_rankings <- rbind(results_rankings,c(rbtlb(1,expit(mu+nu),theta,M,R=J)$rankings))
    },error=function(e){})
    
  }}
results_rankings_preds <- matrix(NA,nrow=0,ncol=4)
for(j1 in 1:J){
  for(j2 in setdiff(1:J,j1)){
    ranks_j1 <- apply(results_rankings,1,function(ranking){which(ranking==j1)})
    ranks_j2 <- apply(results_rankings,1,function(ranking){which(ranking==j2)})
    true_ranks_j1<-apply(rankings,1,function(ranking){
      if(j1 %in% ranking){return(which(ranking==j1))
      }else{return(J+1)}
    })
    true_ranks_j2<-apply(rankings,1,function(ranking){
      if(j2 %in% ranking){return(which(ranking==j2))
      }else{return(J+1)}
    })
    results_rankings_preds <- rbind(
      results_rankings_preds,
      c(j1,j2,mean(ranks_j1<ranks_j2),
        mean((true_ranks_j1 < true_ranks_j2)[(true_ranks_j1 < (J+1))|(true_ranks_j2 < (J+1))]))
    )
  }
}
results_rankings_preds <- as.data.frame(results_rankings_preds)
names(results_rankings_preds) <- c("j1","j2","postpred","obs")

postpred1<- ggplot(melt(postpred_ratings_mean),aes(x=factor(variable),y=value))+
  geom_violin()+theme_bw()+labs(x="Proposal",y="Mean",title="Ratings: Mean")+
  scale_y_continuous(limits=c(0,M))+theme(panel.grid.minor = element_blank())+
  geom_point(data=data.frame(variable=factor(1:J),value=apply(ratings,2,function(x){mean(x,na.rm=T)})),
             aes(x=variable,y=value),color="red")
postpred2<- ggplot(melt(postpred_ratings_var),aes(x=factor(variable),y=value))+
  geom_violin()+theme_bw()+labs(x="Proposal",y="Variance",title="Ratings: Variance")+
  theme(panel.grid.minor = element_blank())+
  geom_point(data=data.frame(variable=factor(1:J),value=apply(ratings,2,function(x){var(x,na.rm=T)})),
             aes(x=variable,y=value),color="red")
postpred3<- ggplot(results_rankings_preds,aes(x=obs,y=postpred))+geom_point()+
  geom_abline(slope=1,intercept=0,linetype=2)+theme_bw()+
  theme(panel.grid.minor = element_blank())+
  labs(x="Observed Probabilities",y="Posterior Predictive Mean Probabilities",
       title="Rankings: Pairwise Probabilities")
AIBS_postpred<-grid.arrange(postpred1,postpred2,postpred3,nrow=1)
ggsave("ResultsPlots/AIBS_postpred.png",AIBS_postpred,width=12,height=3,units="in") 
