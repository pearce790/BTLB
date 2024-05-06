### Initalization ####
source("20240429_source.R")
library(RColorBrewer)

### Results Aggregation: Mechanisms I, II, III ####
iterations <- unlist(lapply(paste0("00",as.character(1:180)),function(val){str_sub(val,-3,-1)}))
results <- matrix(nrow=0,ncol=9)
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/20240429_SimulationResults/res_SimBTLB_20240429_Case1_",iter,".Rdata"))
    
    prob_K <- mean(res$results$K[,1]==res$K)
    MAE_nu <- mean(abs(apply(res$results$nu,2,mean)-res$nu))
    
    reorderings <- permutations(res$K,res$K)
    postmean_mu <- apply(res$results$mu[,1:res$K,res$results$K[,1]==res$K,drop=FALSE],c(1,2),mean)
    order <- reorderings[which.min(apply(reorderings,1,function(order){median(abs(postmean_mu[,order]-res$mu))})),]
    MAE_mu <- mean(abs(postmean_mu[,order]-res$mu))
    
    postmean_theta <- apply(res$results$theta[res$results$K[,1]==res$K,1:res$K,drop=FALSE],2,mean)
    MAE_theta <- mean(abs(postmean_theta[order]-res$theta))

    results <- rbind(results,c(res$model,res$rater_effects,res$I,res$K,res$seed,prob_K,MAE_nu,MAE_mu,MAE_theta))
  }, error = function(e) e)
}
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/20240429_SimulationResults/res_SimBTLB_20240429_Case2_",iter,".Rdata"))
    
    prob_K <- mean(res$results$K[,1]==res$K)
    MAE_nu <- mean(abs(apply(res$results$nu,2,mean)-res$nu))
    
    reorderings <- permutations(res$K,res$K)
    postmean_mu <- apply(res$results$mu[,1:res$K,res$results$K[,1]==res$K,drop=FALSE],c(1,2),mean)
    order <- reorderings[which.min(apply(reorderings,1,function(order){median(abs(postmean_mu[,order]-res$mu))})),]
    MAE_mu <- mean(abs(postmean_mu[,order]-res$mu))
    
    postmean_theta <- apply(res$results$theta[res$results$K[,1]==res$K,1:res$K,drop=FALSE],2,mean)
    MAE_theta <- mean(abs(postmean_theta[order]-res$theta))
    
    results <- rbind(results,c("Misspecified\nRatings",res$rater_effects,res$I,res$K,res$seed,prob_K,MAE_nu,MAE_mu,MAE_theta))
  }, error = function(e) e)
}
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/20240429_SimulationResults/res_SimBTLB_20240429_Case3_",iter,".Rdata"))
    
    prob_K <- mean(res$results$K[,1]==res$K)
    MAE_nu <- mean(abs(apply(res$results$nu,2,mean)-res$nu))
    
    reorderings <- permutations(res$K,res$K)
    postmean_mu <- apply(res$results$mu[,1:res$K,res$results$K[,1]==res$K,drop=FALSE],c(1,2),mean)
    order <- reorderings[which.min(apply(reorderings,1,function(order){median(abs(postmean_mu[,order]-res$mu))})),]
    MAE_mu <- mean(abs(postmean_mu[,order]-res$mu))
    
    postmean_theta <- apply(res$results$theta[res$results$K[,1]==res$K,1:res$K,drop=FALSE],2,mean)
    MAE_theta <- mean(abs(postmean_theta[order]-res$theta))
    
    results <- rbind(results,c("Misspecified\nRankings",res$rater_effects,res$I,res$K,res$seed,prob_K,MAE_nu,MAE_mu,MAE_theta))
  }, error = function(e) e)
}
rm(res)

results <- as.data.frame(results)
names(results) <- c("model","raters","I","K","seed","probK","MAE_nu","MAE_mu","MAE_theta")
results$model <- factor(results$model,levels=c("BTL-Binomial","Misspecified\nRatings","Misspecified\nRankings"))
results$raters <- as.logical(results$raters)
results$I <- factor(as.numeric(results$I))
results$K <- factor(as.numeric(results$K),labels=c(expression(K[0]~"= 1"),expression(K[0]~"= 2"),expression(K[0]~"= 4")))
results$seed <- as.numeric(results$seed)
results$probK <- as.numeric(results$probK)
results$MAE_nu <- as.numeric(results$MAE_nu)
results$MAE_mu <- as.numeric(results$MAE_mu)
results$MAE_theta <- as.numeric(results$MAE_theta)

results[results$model=="Misspecified\nRankings","MAE_theta"] <- NA
plot_results <- melt(results,id.vars=c(1,2,3,4,5))
plot_results$variable<- factor(plot_results$variable,levels=c("probK","MAE_nu","MAE_mu","MAE_theta"),
                               labels=c(expression("P["~K^+~"="~K[0]~"]"),
                                        expression("MAE("~nu~")"),expression("MAE("~mu~")"),expression("MAE("~theta~")")))

ggplot(plot_results,aes(color=I,y=value,x=model))+geom_boxplot(outlier.size = .5,outlier.alpha=0.5)+
  facet_grid(variable~K,scales = "free_y",labeller = label_parsed)+
  scale_color_manual(values=c("skyblue","#3182BD","#08519C"))+
  theme_bw()+theme(legend.position = "bottom",panel.grid.major.x = element_blank(),
                   panel.grid.minor = element_blank())+
  labs(color="Number of Judges, I",y=element_blank(),x=element_blank())
ggsave("ResultsPlots/MFM_Simulation1.pdf",plot=last_plot(),width=9,height=4,units="in") 


### Results Aggregation: Additional Model Fit to Mechanism I####
results <- matrix(nrow=0,ncol=4)
for(iter in iterations){
  tryCatch({
    load(paste0("Section 4.1 Simulation Study/20240429_SimulationResults/res_SimBTLB_20240429_Case4_",iter,".Rdata"))
    MAE_K <- mean(abs(res$results$K[,1]-res$K))
    results <- rbind(results,c(res$I,res$K,res$seed,MAE_K))
  }, error = function(e) e)
}

results <- as.data.frame(results)
names(results) <- c("I","K","seed","MAE_K")
results$I <- factor(as.numeric(results$I))
results$K <- factor(as.numeric(results$K),labels=c(expression(K[0]~"= 1"),expression(K[0]~"= 2"),expression(K[0]~"= 4")))
results$seed <- as.numeric(results$seed)
results$MAE_K <- as.numeric(results$MAE_K)

ggplot(results,aes(x=I,y=MAE_K))+geom_boxplot(outlier.size = .5,outlier.alpha=0.5)+
  facet_grid(~K,labeller = label_parsed)+
  theme_bw()+theme(legend.position = "bottom")+scale_y_continuous(breaks=seq(0,20,by=2))+
  labs(x="Number of Judges, I",y="MAE(K)")
ggsave("ResultsPlots/MFM_Simulation2.pdf",plot=last_plot(),width=9,height=2,units="in") 

