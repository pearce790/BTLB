source("20240429_source.R")

tmp <- matrix(NA,nrow=0,ncol=255)
for(iter in 1:100){
  tryCatch({
    load(paste0("Section 4.2 Conference/20240429_ConferenceResults/res_ConfSim_20240429_",iter,".Rdata"))
    tmp <- rbind(tmp,results)
  }, error = function(e) e)
}
results <- tmp
rm(tmp)


## Plot Results ####
J <- 8
names(results)[4:103] <- names(results[5:104])
names(results)[104] <- "theta"
names(results)[3] <- "iter"

accuracy_mu <- cbind(melt(results[,c(1,2,104,4:53)],id.vars = c(1:3)),
                     melt(results[,c(1,2,104,54:103)],id.vars = c(1:3))[,5])
names(accuracy_mu) <- c("M","R","theta","name","mu","mu_hat")
accuracy_theta <- results[,c(1,2,104,105)]
accuracy_theta$theta <- factor(accuracy_theta$theta)
accuracy_theta_plot <- data.frame(theta=c(1,10,20,40),theta_hat=c(1,10,20,40))
accuracy_theta_plot$theta <- factor(accuracy_theta_plot$theta)

accuracy_nu <- cbind(melt(results[,c(1,2,104,156:205)],id.vars = c(1:3)),
                     melt(results[,c(1,2,104,206:255)],id.vars = c(1:3))[,5])
names(accuracy_nu) <- c("M","R","theta","name","nu","nu_hat")


p1<-ggplot(accuracy_mu,aes(mu,mu_hat))+geom_point(alpha=0.02)+
  geom_abline(slope=1,intercept = 0,color="red")+
  facet_grid(cols=vars(M),rows=vars(R),labeller=label_both)+
  labs(x=expression(mu[0]),y=expression(hat(mu)))+
  theme_bw()+theme(panel.grid.minor = element_blank())
p2<-ggplot(accuracy_theta,aes(factor(theta),theta_map))+geom_violin()+
  geom_point(data=accuracy_theta_plot,aes(theta,theta_hat),color="red",alpha=.5)+
  facet_grid(cols=vars(M),rows=vars(R),labeller = label_both)+
  labs(x=expression(theta[0]),y=expression(hat(theta)))+
  theme_bw()+theme(panel.grid.minor = element_blank())
p3<-ggplot(accuracy_nu,aes(nu,nu_hat))+geom_point(alpha=0.02)+
  geom_abline(slope=1,intercept = 0,color="red")+
  facet_grid(cols=vars(M),rows=vars(R),labeller=label_both)+
  labs(x=expression(nu[0]),y=expression(hat(nu)))+
  theme_bw()+theme(panel.grid.minor = element_blank())
ggsave("ResultsPlots/Conference_res2.png",p2,width=6,height=5)
ggsave("ResultsPlots/Conference_res3.png",grid.arrange(p1,p3,nrow=1),width=12,height=5)

accuracy_mu$error <- accuracy_mu$mu_hat-accuracy_mu$mu
p4<-ggplot(accuracy_mu,aes(factor(theta),error))+geom_violin()+
  facet_grid(cols=vars(M),rows=vars(R),labeller = label_both)+
  geom_abline(slope=0,intercept = 0,color="red",linetype=2)+
  labs(x=expression(theta[0]),y=expression(hat(mu)-mu[0]))+
  theme_bw()+theme(panel.grid.minor = element_blank())
accuracy_nu$error <- accuracy_nu$nu_hat-accuracy_nu$nu
p5<-ggplot(accuracy_nu,aes(factor(theta),error))+geom_violin()+
  facet_grid(cols=vars(M),rows=vars(R),labeller = label_both)+
  geom_abline(slope=0,intercept = 0,color="red",linetype=2)+
  labs(x=expression(theta[0]),y=expression(hat(nu)-nu[0]))+
  theme_bw()+theme(panel.grid.minor = element_blank())
accuracy_theta$error <- accuracy_theta$theta_map-as.numeric(as.character(accuracy_theta$theta))
p6<-ggplot(accuracy_theta,aes(theta,error))+geom_violin()+
  facet_grid(cols=vars(M),rows=vars(R),labeller = label_both)+
  geom_abline(slope=0,intercept = 0,color="red",linetype=2)+
  labs(x=expression(theta[0]),y=expression(hat(theta)-theta[0]))+
  theme_bw()+theme(panel.grid.minor = element_blank())
ggsave("ResultsPlots/Conference_res4.png",p6,width=6,height=5)
ggsave("ResultsPlots/Conference_res5.png",grid.arrange(p4,p5,nrow=1),width=12,height=5)

btlb_distance <- apply(results,1,function(res){
  rankrate::kendall(order(rank(res[54:103],ties.method="random")),order(res[4:53]))})
standard_distance <- apply(results,1,function(res){
  rankrate::kendall(order(rank(res[106:155],ties.method="random")),order(res[4:53]))})
diff_pi <- cbind(results[,c(1,2,3,104)],btlb_distance,standard_distance)
diff_pi$btlb_distance <- diff_pi$btlb_distance/choose(50,2)
diff_pi$standard_distance <- diff_pi$standard_distance/choose(50,2)

ggplot(diff_pi,aes(standard_distance,btlb_distance,color=factor(theta)))+
  geom_point(size=1,alpha=.3)+
  facet_grid(cols=vars(R),rows=vars(M),labeller=label_both)+
  geom_abline(slope=1,intercept=0,color="red",linetype=2)+
  labs(x=expression(paste("Inaccuracy of ",hat(pi)[0]^X," to ",pi[0])),
       y=expression(paste("Inaccuracy of ",hat(pi)[0]^BTLB," to ",pi[0])),
       color=expression(theta))+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes=list(alpha=1,size=3)))
ggsave("ResultsPlots/Conference_res1.png",last_plot(),width=12,height=5)





