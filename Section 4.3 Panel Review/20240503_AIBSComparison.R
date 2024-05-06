### Functions ####
source("20240429_source.R")
load("Section 4.3 Panel Review/Data/AIBS_Clean.RData")

### General Constants and Hyperparameters ####
I <- nrow(rankings)
J <- ncol(rankings)

delta_gamma <- c(3,1)
beta_mu <- c(5,10)
beta_nu <- c(5,1)
mu_0j <- rep(0,J)
delta_theta <- c(10,0.5)
sigma_mu=0.5;sigma_nu=0.2;sigma_theta=2;sigma_gamma=1


## K=1 Case ####
T1<-500;T2<-5; K<-1
set.seed(20240502)
res_ALL_K2<-btlb_lcmm(rankings,ratings,M,K,delta_theta,delta_gamma,beta_mu,beta_nu,mu_0j,
                      T1,T2,sigma_mu,sigma_nu,sigma_theta,sigma_gamma)

set.seed(20240502)
res_RANKINGS_K2<-btl_lcmm(rankings,M,K,delta_gamma,beta_mu,mu_0j,
                          T1,T2,sigma_mu,sigma_gamma)
set.seed(20240502)
res_RATINGS_K2<-binom_lcmm(ratings,M,K,delta_gamma,beta_mu,beta_nu,mu_0j,
                           T1,T2,sigma_mu,sigma_nu,sigma_gamma)

ggplot(melt(res_ALL_K2$mu),aes(x=Var3,y=value,color=factor(Var1)))+geom_line()+facet_grid(~Var2)
ggplot(melt(res_RATINGS_K2$mu),aes(x=Var3,y=value,color=factor(Var1)))+geom_line()+facet_grid(~Var2)
ggplot(melt(res_RANKINGS_K2$mu),aes(x=Var3,y=value,color=factor(Var1)))+geom_line()+facet_grid(~Var2)

ranks_ALL <- melt(apply(res_ALL_K2$mu[,1,(T1/2+1):T1],2,rank))
ranks_RATINGS <- melt(apply(res_RATINGS_K2$mu[,1,(T1/2+1):T1],2,rank))
ranks_RANKINGS <- melt(apply(res_RANKINGS_K2$mu[,1,(T1/2+1):T1],2,rank))
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
ggsave("ResultsPlots/AIBS_rankcomparisonK1.png",AIBS_rankcomp1,width=10,height=5,units="in") 

## K=2 Case ####
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

ggplot(melt(res_ALL_K2$Z[(T1/2+1):T1,]),aes(x=factor(Var2),fill=factor(value)))+geom_bar(position="fill")
ggplot(melt(res_RATINGS_K2$Z[(T1/2+1):T1,]),aes(x=factor(Var2),fill=factor(value)))+geom_bar(position="fill")
ggplot(melt(res_RANKINGS_K2$Z[(T1/2+1):T1,]),aes(x=factor(Var2),fill=factor(value)))+geom_bar(position="fill")

ggplot(melt(res_ALL_K2$mu),aes(x=Var3,y=value,color=factor(Var1)))+geom_line()+facet_grid(~Var2)
ggplot(melt(res_RATINGS_K2$mu),aes(x=Var3,y=value,color=factor(Var1)))+geom_line()+facet_grid(~Var2)
ggplot(melt(res_RANKINGS_K2$mu),aes(x=Var3,y=value,color=factor(Var1)))+geom_line()+facet_grid(~Var2)

big_ALL <- which.max(apply(res_ALL_K2$alpha,2,mean))
big_RATINGS <- which.max(apply(res_RATINGS_K2$alpha,2,mean))
big_RANKINGS <- which.max(apply(res_RANKINGS_K2$alpha,2,mean))

ranks_ALL <- melt(apply(res_ALL_K2$mu[,big_ALL,(T1/2+1):T1],2,rank))
order_ranks_ALL <- ranks_ALL %>% group_by(Var1) %>% summarize(meanrank=mean(value)) %>% arrange(meanrank) %>%pull(Var1)
ggplot(ranks_ALL,aes(x=factor(Var1,levels=order_ranks_ALL),y=value))+geom_boxplot(outlier.alpha = 0.1)

ranks_RATINGS <- melt(apply(res_RATINGS_K2$mu[,big_RATINGS,(T1/2+1):T1],2,rank))
order_ranks_RATINGS <- ranks_RATINGS %>% group_by(Var1) %>% summarize(meanrank=mean(value)) %>% arrange(meanrank) %>%pull(Var1)
ggplot(ranks_RATINGS,aes(x=factor(Var1,levels=order_ranks_RATINGS),y=value))+geom_boxplot(outlier.alpha = 0.1)

ranks_RANKINGS <- melt(apply(res_RANKINGS_K2$mu[,big_RANKINGS,(T1/2+1):T1],2,rank))
order_ranks_RANKINGS <- ranks_RANKINGS %>% group_by(Var1) %>% summarize(meanrank=mean(value)) %>% arrange(meanrank) %>%pull(Var1)
ggplot(ranks_RANKINGS,aes(x=factor(Var1,levels=order_ranks_RANKINGS),y=value))+geom_boxplot(outlier.alpha = 0.1)

ranks_plot <- cbind(Model=rep(c("BTL-Binomial","BTL","Binomial"),each=nrow(ranks_ALL)),
                    rbind(ranks_ALL,ranks_RANKINGS,ranks_RATINGS))
ranks_plot$Model <- factor(ranks_plot$Model,levels=c("BTL-Binomial","BTL","Binomial"))
AIBS_rankcomp <- ggplot(ranks_plot,aes(x=factor(Var1,levels=order_ranks_ALL),color=Model,y=value))+
  geom_boxplot(outlier.alpha = 0.05)+theme_bw()+
  scale_y_continuous(breaks=seq(1,25,by=2))+
  labs(x="Proposal",y="Rank",color=element_blank())+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x=element_blank())
ggsave("ResultsPlots/AIBS_rankcomparisonK2.png",AIBS_rankcomp,width=10,height=5,units="in") 


## Saving Analysis Results .RData
save.image("Section 4.3 Panel Review/20240503_AIBSComparison.RData")