library(tidyverse)
source("20240429_source.R")

results <- data.frame(iter=1:80,K=NA,loglik=NA)
for(iter in 1:80){
  tryCatch({
    load(paste0("~/BTLBinomial/Sushi Results 20240429/res_Sushi20240429_MAP_iter",iter,".Rdata"))
    results[iter,] <- c(iter,res$map_res$K,res$map_res$loglikelihood)
    rm(res,iter)
  },error=function(e){print(iter)})
}
results %>% group_by(K) %>% summarize(count=n(),max_loglik=max(loglik))


