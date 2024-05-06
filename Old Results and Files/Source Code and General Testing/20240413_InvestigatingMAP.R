library(tidyverse)
source("20240412_source.R")

results <- data.frame(iter=1:80,K=NA,loglik=NA)
for(iter in 1:80){
  tryCatch({
    load(paste0("~/BTLBinomial/Sushi Results/res_Sushi20240412_MAP_iter",iter,".Rdata"))
    results[iter,] <- c(iter,res$map_res$K,res$map_res$loglikelihood)
    rm(res,iter)
  },error=function(e){print(iter)})
}
results %>% group_by(K) %>% summarize(count=n(),max_loglik=max(loglik))


iter<-12
load(paste0("~/BTLBinomial/Sushi Results/res_Sushi20240412_MAP_iter",iter,".Rdata"))
Z_6 <- apply(res$map_res$z_hat,1,which.max)

iter<-40
load(paste0("~/BTLBinomial/Sushi Results/res_Sushi20240412_MAP_iter",iter,".Rdata"))
Z_8 <- apply(res$map_res$z_hat,1,which.max)

iter<-43
load(paste0("~/BTLBinomial/Sushi Results/res_Sushi20240412_MAP_iter",iter,".Rdata"))
Z_12 <- apply(res$map_res$z_hat,1,which.max)

iter<-77
load(paste0("~/BTLBinomial/Sushi Results/res_Sushi20240412_MAP_iter",iter,".Rdata"))
Z_16 <- apply(res$map_res$z_hat,1,which.max)

table(Z_16,Z_12)

table(apply(res$map_res$z_hat,1,which.max))
res$map_res$alpha
res$map_res$gamma
res$map_res$sigma2_mu
res$map_res$theta
round(expit(res$map_res$mu),3)[,c(9,11)]
min(dist(t(expit(res$map_res$mu))))


load("/Users/michaelpearce/BTLBinomial/Sushi Results/res_Sushi20240412_MAP_iter78.Rdata")