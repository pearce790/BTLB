#### Instructions ####
# This .R script may be run in its entirety to replicate the results presented in 
# Section 4.2 of the manuscript and associated Appendix B.2. Note that running this 
# file involves estimating 3200 unique models to data, which may be time-consuming 
# and slow. Specifically, if the script is run without parallelization, it would 
# take approximately 1--2 days to run on a standard but high-powered 2023 MacBook 
# Pro. Parallelization along the outer loop in the code that follows can reduce 
# total computation time to approximately 1 hour. Associated slurm files are 
# included in the "slurm files" folder. 

# Additionally, note that the output files from our 100 simulation iterations are
# included in the provided "results_ConferenceSim" folder, which allows the user 
# to skip model fitting entirely and replicate the creation of plots by skipping 
# directly to the "Figures" section of the .R script, if desired.

#### Load Source Files ####
source("source.R")

#### Code to Run Simulations ####
for(iter in 1:100){

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
  
  # iterating over M, R, theta 
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
      }}}
  
  results <- as.data.frame(results)
  names(results) <- c("M","R","iter",paste0("mu",1:50),paste0("mu_map",1:50),
                      "theta","theta_map",paste0("Xbar",1:50),paste0("nu",1:50),paste0("nu_map",1:50))
  save(results,file=paste0("Section 4.2 Conference/simulation_results/res_ConfSim_20240429_",iter,".Rdata"))
}

#### Figures ####

## load the stored MAP estimates from simulations above
tmp <- matrix(NA,nrow=0,ncol=255)
for(iter in 1:100){
  load(paste0("Section 4.2 Conference/results_ConferenceSim/res_",iter,".Rdata"))
  tmp <- rbind(tmp,results)
}
tmp <- as.data.frame(tmp)
names(tmp) <- attr(results,"names")[1:255]
results <- tmp
J <- 8

## Main Paper Figure 3 ####
btlb_distance <- apply(results,1,function(res){
  rankrate::kendall(order(rank(res[54:103],ties.method="random")),order(res[4:53]))})
standard_distance <- apply(results,1,function(res){
  rankrate::kendall(order(rank(res[106:155],ties.method="random")),order(res[4:53]))})
diff_pi <- cbind(results[,c(1,2,3,104)],btlb_distance,standard_distance)
diff_pi$btlb_distance <- diff_pi$btlb_distance/choose(50,2)
diff_pi$standard_distance <- diff_pi$standard_distance/choose(50,2)
p7<-ggplot(diff_pi,aes(standard_distance,btlb_distance,color=factor(theta)))+
  geom_point(size=1,alpha=.3)+
  facet_grid(cols=vars(R),rows=vars(M),labeller=label_both)+
  geom_abline(slope=1,intercept=0,color="red",linetype=2)+
  labs(x=expression(paste("Inaccuracy of ",hat(pi)[0]^X," to ",pi[0])),
       y=expression(paste("Inaccuracy of ",hat(pi)[0]^BTLB," to ",pi[0])),
       color=expression(theta))+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes=list(alpha=1,size=3)))
ggsave("ResultsPlots/Fig3.png",p7,width=12,height=5)

## Appendix Figure 1 ####
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

ggsave("ResultsPlots/Appendix_Fig1.png",
       grid.arrange(grid.arrange(p1,p3,nrow=1),
                    grid.arrange(NULL,p2,NULL,nrow=1,widths=c(.2,.6,.2)),
                    nrow=2),
       width=12,height=10)

## Appendix Figure 2 ####
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
ggsave("ResultsPlots/Appendix_Fig2.png",
       grid.arrange(grid.arrange(p4,p5,nrow=1),
                    grid.arrange(NULL,p6,NULL,nrow=1,widths=c(.2,.6,.2)),
                    nrow=2),
       width=12,height=10)



