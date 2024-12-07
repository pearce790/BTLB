#### Instructions ####
# This .R script may be run in its entirety to replicate the results presented in 
# Section 4.3 of the manuscript and associated Appendices B.3 and B.4. Note that 
# running this file may be slow due to the large number of MCMC iterations we chose 
# to run. Still, on a relatively high-powered but standard MacBook Pro, this script 
# took approximately 45 minutes to run without any parallelization.

#### Load Source Files ####
source("source.R")
load("Section 4.3 Panel Review/Data/AIBS_Clean.RData")

#### Data Analyses ####
I <- nrow(rankings)
J <- ncol(rankings)
delta_gamma <- c(2,4)
beta_mu <- c(5,10)
beta_nu <- c(5,1)
lambda <- 2
mu_0j <- rep(0,J)
delta_theta <- c(10,0.5)

## MAP Estimation for K=2 and K=3
set.seed(1)
map_res <- replicate(10,btlb_map(rankings,ratings,M,K=2,tol=1,verbose=FALSE,
                                delta_gamma = delta_gamma,beta_mu=beta_mu,beta_nu=beta_nu,
                                mu_0j=mu_0j,delta_theta = delta_theta))
map_hat_K2 <- map_res[,which.max(map_res["loglikelihood",])]
map_res <- replicate(10,btlb_map(rankings,ratings,M,K=3,tol=1,verbose=FALSE,
                                 delta_gamma = delta_gamma,beta_mu=beta_mu,beta_nu=beta_nu,
                                 mu_0j=mu_0j,delta_theta = delta_theta))
map_hat_K3 <- map_res[,which.max(map_res["loglikelihood",])]

## MFMM Estimation with 3 chains (each starting with K=2, 3, or 17)
set.seed(2)
mfmm_mcmc_K2 <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=map_hat_K2$K,init=map_hat_K2,
                          delta_theta=delta_theta,delta_gamma=delta_gamma,lambda=lambda,
                          beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=10000,T2=4,sigma_mu=0.5,sigma_nu=0.2,sigma_theta=2,sigma_gamma=1)
mfmm_mcmc_K3 <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=map_hat_K3$K,init=map_hat_K3,
                          delta_theta=delta_theta,delta_gamma=delta_gamma,lambda=lambda,
                          beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                          T1=10000,T2=4,sigma_mu=0.5,sigma_nu=0.2,sigma_theta=2,sigma_gamma=1)
mfmm_mcmc_K17 <- btlb_mfmm(rankings=rankings,ratings=ratings,M=M,K_start=I,init=NULL,
                       delta_theta=delta_theta,delta_gamma=delta_gamma,lambda=lambda,
                       beta_mu=beta_mu,beta_nu=beta_nu,mu_0j=mu_0j,
                       T1=10000,T2=4,sigma_mu=0.5,sigma_nu=0.2,sigma_theta=2,sigma_gamma=1)

# Burn-in and label-switching
mfmm_mcmc <- mfmm_mcmc_K2
keep_iters <- c(6001:10000,16001:20000,26001:30000) #burn first 60% of each chain
mfmm_mcmc$gamma <- rbind(mfmm_mcmc$gamma,mfmm_mcmc_K3$gamma,mfmm_mcmc_K17$gamma)[keep_iters,,drop=FALSE]
mfmm_mcmc$nu <- rbind(mfmm_mcmc$nu,mfmm_mcmc_K3$nu,mfmm_mcmc_K17$nu)[keep_iters,,drop=FALSE]
mfmm_mcmc$sigma2 <- rbind(mfmm_mcmc$sigma2,mfmm_mcmc_K3$sigma2,mfmm_mcmc_K17$sigma2)[keep_iters,,drop=FALSE]
mfmm_mcmc$Z <- rbind(mfmm_mcmc$Z,mfmm_mcmc_K3$Z,mfmm_mcmc_K17$Z)[keep_iters,,drop=FALSE]
mfmm_mcmc$K <- rbind(mfmm_mcmc$K,mfmm_mcmc_K3$K,mfmm_mcmc_K17$K)[keep_iters,,drop=FALSE]
mfmm_mcmc$loglik <- rbind(mfmm_mcmc$loglik,mfmm_mcmc_K3$loglik,mfmm_mcmc_K17$loglik)[keep_iters,,drop=FALSE]
mfmm_mcmc$accept <- rbind(mfmm_mcmc$accept,mfmm_mcmc_K3$accept,mfmm_mcmc_K17$accept)[keep_iters,,drop=FALSE]
mfmm_mcmc$timing <- rbind(mfmm_mcmc$timing,mfmm_mcmc_K3$timing,mfmm_mcmc_K17$timing)[keep_iters,,drop=FALSE]
mfmm_mcmc$alpha <- rbind(cbind(mfmm_mcmc_K2$alpha,matrix(NA,nrow=10000,ncol=17-10)),
                         cbind(mfmm_mcmc_K3$alpha,matrix(NA,nrow=10000,ncol=17-9)),
                         mfmm_mcmc_K17$alpha)[keep_iters,,drop=FALSE]
mfmm_mcmc$theta <- rbind(cbind(mfmm_mcmc_K2$theta,matrix(NA,nrow=10000,ncol=17-10)),
                         cbind(mfmm_mcmc_K3$theta,matrix(NA,nrow=10000,ncol=17-9)),
                         mfmm_mcmc_K17$theta)[keep_iters,,drop=FALSE]
mfmm_mcmc$mu <- abind(abind(mfmm_mcmc_K2$mu,array(NA,dim=c(25,17-10,10000)),along=2),
                      abind(mfmm_mcmc_K3$mu,array(NA,dim=c(25,17-9,10000)),along=2),
                      mfmm_mcmc_K17$mu)[,,keep_iters,drop=FALSE]
mfmm_mcmc$class <- abind(abind(mfmm_mcmc_K2$class,array(0,dim=c(10000,17,17-10)),along=3),
                         abind(mfmm_mcmc_K3$class,array(0,dim=c(10000,17,17-9)),along=3),
                         mfmm_mcmc_K17$class,along=1)[keep_iters,,,drop=FALSE]
mfmm_mcmc <- stephens_labelswitching(mfmm_mcmc,burn_prop = 0,thin=1) # fix label switching in combined chain
mfmm_mcmc$class <- NULL # saving space!
rm(mfmm_mcmc_K2,mfmm_mcmc_K3,mfmm_mcmc_K17) # saving space!


#### Figures ####

# save color palette
colors <- rev(c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#2171b5","#08306b"))

## Main Paper Figure 4 ####
order_meanrating <- order(apply(ratings,2,function(x){mean(x,na.rm=T)}))
pEDA1<-ggplot(na.exclude(melt(ratings)),aes(x=factor(Var2,levels=order_meanrating),y=value))+
  geom_boxplot()+theme_bw()+
  scale_y_continuous(limits=c(0,M))+
  labs(x=element_blank(),y="Rating")+theme(panel.grid.minor = element_blank())
pEDA2<-ggplot(na.exclude(melt(rankings)),aes(x=factor(value,levels=order_meanrating),fill=factor(Var2,levels=6:1)))+
  geom_bar(stat="count",color="gray50")+theme_bw()+
  scale_y_continuous(limits=c(0,14),breaks=seq(0,14,by=2))+
  scale_x_discrete(drop=FALSE)+
  scale_fill_manual(values=c("#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#005A32"))+
  labs(x=element_blank(),y="Count",fill="Rank")+
  theme(legend.position="bottom",panel.grid.minor = element_blank())+
  guides(fill = guide_legend(nrow = 1,reverse=TRUE))
pEDA <- grid.arrange(pEDA1,pEDA2,nrow=2,heights=c(.45,.55))  
ggsave("ResultsPlots/Fig4.png",pEDA,width=12,height=5,units="in")

## Main Paper Figure 5 ####
plotK<-ggplot(data.frame(clusters=mfmm_mcmc$K[,1]),aes(x=clusters))+geom_bar()+
  theme_bw()+labs(x="K+",y="Probability")+
  scale_y_continuous(breaks=seq(0,nrow(mfmm_mcmc$K),length=5),labels=seq(0,1,length=5))
plotZ<-ggplot(melt(mfmm_mcmc$Z),
                aes(x=factor(Var2),
                    fill=factor(value,levels=c(1:6),labels=1:6),
                    color=factor(value,levels=c(1:6),labels=1:6)))+
  geom_bar()+theme_minimal()+
  labs(x="Judges",y="Proportion",fill="Class",color="Class")+
  scale_y_continuous(breaks=seq(0,nrow(mfmm_mcmc$Z),length=5),labels=seq(0,1,length=5))+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  theme(panel.grid = element_blank(),legend.position = "bottom",
        axis.ticks.x = element_blank(),axis.ticks.y=element_line())
plot_mu<-ggplot(melt(mfmm_mcmc$mu[,1:2,]),
                aes(x=factor(Var1,levels=order(apply(mfmm_mcmc$mu[,1,],1,median))),
                    y=value,color=factor(Var2),fill=factor(Var2)))+
  geom_violin()+theme_bw()+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  geom_vline(xintercept=seq(1.5,nrow(mfmm_mcmc$mu)-.5),color="gray90")+
  labs(x="Proposal",y=expression(mu[j]),fill="Cluster",color="Cluster")+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom")
plot_theta<-ggplot(melt(mfmm_mcmc$theta[,1:2]),
                   aes(x=factor(Var2),y=value,color=factor(Var2),fill=factor(Var2)))+
  geom_violin()+theme_bw()+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  labs(x="Class",y=expression(theta),fill="Class",color="Class")+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom")
plot1 <- grid.arrange(
  grid.arrange(plotK,plotZ+theme(legend.position = "right"),widths=c(0.2,0.8)),
  grid.arrange(plot_mu+theme(legend.position = "none"),plot_theta+theme(legend.position = "right"),
               widths=c(0.8,0.2))
)
ggsave("ResultsPlots/Fig5.png",plot1,width=12,height=5,units="in")

## Main Paper Figure 6 ####
class <- data.frame(Var2=1:I,pred_class=round(apply(mfmm_mcmc$Z,2,median)))
nu_long <- melt(mfmm_mcmc$nu) %>% left_join(class,by="Var2")
plot_nu<-ggplot(nu_long,aes(y=factor(Var2,levels=order(apply(mfmm_mcmc$nu,2,median))),x=value,
                            fill=factor(pred_class)))+
  geom_density_ridges(bandwidth=.04,alpha=0.95)+theme_minimal()+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  scale_x_continuous(limits=c(-.75,0.75))+
  labs(y="Judge",x="Leniency (low) vs. Harshness (high)",fill="Class")+
  theme(panel.grid = element_blank(),legend.position = "right")
ggsave("ResultsPlots/Fig6.png",plot_nu,width=8,height=4,units="in")

## Main Paper Figure 7 ####
sigma_mu=0.5;sigma_nu=0.2;sigma_theta=2;sigma_gamma=1
T1<-1000;T2<-5; K<-2
set.seed(20240502)
res_ALL_K2<-btlb_lcmm(rankings,ratings,M,K,delta_theta,delta_gamma,beta_mu,beta_nu,mu_0j,
                      T1,T2,sigma_mu,sigma_nu,sigma_theta,sigma_gamma)
set.seed(20240502)
res_RANKINGS_K2<-btl_lcmm(rankings,M,K,delta_gamma,beta_mu,mu_0j,
                          T1,T2,sigma_mu,sigma_gamma)
set.seed(20240502)
res_RATINGS_K2<-binom_lcmm(ratings,M,K,delta_gamma,beta_mu,beta_nu,mu_0j,
                           T1,T2,sigma_mu,sigma_nu,sigma_gamma)
big_ALL <- which.max(apply(res_ALL_K2$alpha,2,mean))
big_RATINGS <- which.max(apply(res_RATINGS_K2$alpha,2,mean))
big_RANKINGS <- which.max(apply(res_RANKINGS_K2$alpha,2,mean))
ranks_ALL <- melt(apply(res_ALL_K2$mu[,big_ALL,(T1/2+1):T1],2,rank))
order_ranks_ALL <- ranks_ALL %>% group_by(Var1) %>% summarize(meanrank=mean(value)) %>% arrange(meanrank) %>%pull(Var1)
ranks_RATINGS <- melt(apply(res_RATINGS_K2$mu[,big_RATINGS,(T1/2+1):T1],2,rank))
ranks_RANKINGS <- melt(apply(res_RANKINGS_K2$mu[,big_RANKINGS,(T1/2+1):T1],2,rank))
ranks_plot <- cbind(Model=rep(c("BTL-Binomial","BTL","Binomial"),each=nrow(ranks_ALL)),
                    rbind(ranks_ALL,ranks_RANKINGS,ranks_RATINGS))
ranks_plot$Model <- factor(ranks_plot$Model,levels=c("BTL-Binomial","BTL","Binomial"))
AIBS_rankcomp <- ggplot(ranks_plot,aes(x=factor(Var1,levels=order_ranks_ALL),color=Model,y=value))+
  geom_boxplot(outlier.alpha = 0.05)+theme_bw()+
  scale_y_continuous(breaks=seq(1,25,by=2))+
  labs(x="Proposal",y="Rank",color=element_blank())+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x=element_blank())
ggsave("ResultsPlots/Fig7.png",AIBS_rankcomp,width=10,height=5,units="in") 

## Appendix Figure 3 ####
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
ggsave("ResultsPlots/Appendix_Fig3.png",prior_plot,width=10,height=3,units="in") 

## Appendix Figure 4 ####
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
ggsave("ResultsPlots/Appendix_Fig4.png",AIBS_param2,width=10,height=3,units="in") 


## Appendix Figure 5 ####
g1<-ggplot(melt(apply(mfmm_mcmc$accept,2,cummean)),aes(Var1,value,group=Var2,color=Var2))+
  theme_bw()+scale_color_brewer(type="seq",palette="Dark2")+
  scale_y_continuous(limits=c(0,1))+scale_x_continuous(breaks=seq(0,12000,by=2000))+geom_line()+
  labs(color="Variable",x="Iteration",y="Cumulative Acceptance Probability",title="Metropolis-Hastings Acceptance Probabilities")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())
g2<-ggplot(melt(mfmm_mcmc$timing),aes(Var1,value,group=Var2,color=Var2))+
  theme_bw()+scale_color_brewer(type="seq",palette="Dark2")+geom_line(alpha=0.5)+
  scale_x_continuous(breaks=seq(0,12000,by=2000))+
  scale_y_log10(breaks = trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x)))+
  labs(color="Algorithm Step",x="Iteration",y="Time (seconds)",title="Algorithm Speed by Iteration")+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())
AIBS_alg <- grid.arrange(g1,g2,nrow=1)
ggsave("ResultsPlots/Appendix_Fig5.png",AIBS_alg,width=10,height=4,units="in") 

## Appendix Figure 6 ####
iters <- seq(1,12000,by=10)
postpred_ratings_mean <- matrix(NA,nrow=length(iters),ncol=J)
for(j in 1:J){
  ratings_j<-c()
  for(iter in iters){
    mus <- mfmm_mcmc$mu[j,mfmm_mcmc$Z[iter,],iter] 
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
    mus <- mfmm_mcmc$mu[j,mfmm_mcmc$Z[iter,],iter] 
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
    class <- mfmm_mcmc$Z[iter,i]
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
ggsave("ResultsPlots/Appendix_Fig6.png",AIBS_postpred,width=12,height=3,units="in") 

## Appendix Figure 7 ####
melt_K <- melt(mfmm_mcmc$K)
melt_K$Var1 <- rep(1:(nrow(melt_K)/2),2)
traceK<-ggplot(melt_K,aes(Var1,value,color=factor(Var2,levels=c("K","Kplus"),labels=c("K","K+"))))+
  geom_line(alpha=0.75)+theme_bw()+
  scale_x_continuous(breaks=seq(0,12000,by=2000))+
  scale_y_continuous(limits=c(1,max(mfmm_mcmc$K)),breaks=1:max(mfmm_mcmc$K))+
  scale_color_manual(values=colors[c(3,1)])+
  labs(x="Iteration",y="Value",color="",title="Trace Plot: K, K+")+
  theme(legend.position = "right",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
tracealpha<-ggplot(melt(mfmm_mcmc$alpha[,1:5]),aes(Var1,value,color=factor(Var2)))+
  geom_line()+theme_bw()+
  scale_color_manual(values=colors)+
  scale_x_continuous(breaks=seq(0,12000,by=2000))+
  labs(x="Iteration",y="Value",color="k",title=expression("Trace Plot:"~alpha))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
melt_sigma <- melt(mfmm_mcmc$sigma2)
melt_sigma$Var1 <- rep(1:(nrow(melt_sigma)/2),2)
tracesigma<-ggplot(melt_sigma,aes(Var1,value,color=factor(Var2)))+
  geom_line()+theme_bw()+
  scale_x_continuous(breaks=seq(0,12000,by=2000))+
  scale_color_manual(values=colors[c(3,1)],labels=c(expression(sigma[mu]^2),expression(sigma[nu]^2)))+
  labs(x="Iteration",y="Value",color="",title="Trace Plot: Variance Hyperparameters")+
  theme(legend.position = "right",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
tracegamma<-ggplot(melt(mfmm_mcmc$gamma),aes(Var1,value))+
  geom_line()+theme_bw()+
  scale_x_continuous(breaks=seq(0,12000,by=2000))+
  labs(x="Iteration",y="Value",color="",title=expression("Trace Plot:"~gamma))+
  theme(legend.position = "right",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
AIBS_trace1<-grid.arrange(tracegamma,tracesigma,tracealpha,traceK)
AIBS_trace2<-ggplot(melt(mfmm_mcmc$theta[,1:2]),aes(x=Var1,color=factor(Var2),y=value))+
  geom_line(alpha=c(rep(1,12000),rep(.5,12000)))+theme_bw()+
  scale_color_manual(values=colors)+
  scale_x_continuous(breaks=seq(0,12000,by=2000))+
  labs(x="Iteration",y="Value",title=expression("Trace Plot: "~theta),color="Class")
AIBS_trace3 <- grid.arrange(AIBS_trace1,AIBS_trace2,nrow=2,heights=c(.667,.333))
ggsave("ResultsPlots/Appendix_Fig7.png",AIBS_trace3,width=10,height=7.5,units="in") 

## Appendix Figure 8 ####
AIBS_trace1<-ggplot(melt(mfmm_mcmc$mu[,1:2,seq(1,12000,by=10)]),aes(x=Var3,color=factor(Var2),y=value))+
  geom_line(alpha=0.8)+theme_bw()+
  scale_x_continuous(breaks=seq(0,1200,by=300))+
  scale_color_manual(values=colors)+
  facet_wrap(~Var1)+labs(x="Iteration",y="Value",title=expression("Trace Plot: "~mu),color="Class")
AIBS_trace2<-ggplot(melt(mfmm_mcmc$nu[seq(1,12000,by=10),]),aes(x=Var1,y=value))+
  scale_x_continuous(breaks=seq(0,1200,by=300))+
  geom_line(alpha=0.6)+theme_bw()+
  facet_wrap(~Var2,nrow=3)+labs(x="Iteration",y="Value",title=expression("Trace Plot: "~nu))
AIBS_trace3<-grid.arrange(AIBS_trace1,AIBS_trace2,nrow=2,heights=c(0.5,0.5))
ggsave("ResultsPlots/Appendix_Fig8.png",AIBS_trace3,width=10,height=10,units="in") 

## Appendix Figure 9 ####
T1<-1000;T2<-5; K<-1
set.seed(20240502)
res_ALL_K1<-btlb_lcmm(rankings,ratings,M,K,delta_theta,delta_gamma,beta_mu,beta_nu,mu_0j,
                      T1,T2,sigma_mu,sigma_nu,sigma_theta,sigma_gamma)
set.seed(20240502)
res_RANKINGS_K1<-btl_lcmm(rankings,M,K,delta_gamma,beta_mu,mu_0j,
                          T1,T2,sigma_mu,sigma_gamma)
set.seed(20240502)
res_RATINGS_K1<-binom_lcmm(ratings,M,K,delta_gamma,beta_mu,beta_nu,mu_0j,
                           T1,T2,sigma_mu,sigma_nu,sigma_gamma)
ranks_ALL <- melt(apply(res_ALL_K1$mu[,1,(T1/2+1):T1],2,rank))
ranks_RATINGS <- melt(apply(res_RATINGS_K1$mu[,1,(T1/2+1):T1],2,rank))
ranks_RANKINGS <- melt(apply(res_RANKINGS_K1$mu[,1,(T1/2+1):T1],2,rank))
order_ranks_ALL <- ranks_ALL %>% group_by(Var1) %>% summarize(meanrank=mean(value)) %>% arrange(meanrank) %>%pull(Var1)
ranks_plot <- cbind(Model=rep(c("BTL-Binomial","BTL","Binomial"),each=nrow(ranks_ALL)),
                    rbind(ranks_ALL,ranks_RANKINGS,ranks_RATINGS))
ranks_plot$Model <- factor(ranks_plot$Model,levels=c("BTL-Binomial","BTL","Binomial"))
AIBS_rankcomp1 <- ggplot(ranks_plot,aes(x=factor(Var1,levels=order_ranks_ALL),color=Model,y=value))+
  geom_boxplot(outlier.alpha = 0.05)+theme_bw()+
  scale_y_continuous(breaks=seq(1,25,by=2))+
  labs(x="Proposal",y="Rank",color=element_blank())+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x=element_blank())
ggsave("ResultsPlots/Appendix_Fig9.png",AIBS_rankcomp1,width=10,height=5,units="in") 

