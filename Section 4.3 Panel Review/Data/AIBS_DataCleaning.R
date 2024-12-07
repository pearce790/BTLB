## Data Import/Cleaning ####
library(dplyr)
library(readxl)

rawdata <- read_excel("Section 4.3 Panel Review/Data/GrantPeerReview_Scoring_Ranking_DataA2.xlsx", sheet = "Pnl2")
raw_assignments <- read_excel("Section 4.3 Panel Review/Data/GrantPeerReview_Scoring_Ranking_DataA2.xlsx", sheet = "Pnl2_Assignments")
merged_data <- as.data.frame(rawdata[2:18,2:27]) %>% left_join(as.data.frame(rawdata[28:40,2:27]), by="RandomPropID")
judge_IDS <- merged_data$RandomPropID
paper_IDS <- names(rawdata[,3:27])
ratings <- merged_data[,2:26]
orderings <- merged_data[,27:51]
for(j in 1:ncol(ratings)){ratings[,j] <- as.numeric(ratings[,j])}
ratings <- as.matrix(10*(ratings-1))
colnames(ratings) <- NULL
rankings <- matrix(NA,nrow=nrow(ratings),ncol=ncol(ratings))
rankings[,1] <- unlist(apply(orderings,1,function(x){if("First Place" %in% x){return(which(x=="First Place"))}else{return(NA)}}))
rankings[,2] <- unlist(apply(orderings,1,function(x){if("Second Place" %in% x){return(which(x=="Second Place"))}else{return(NA)}}))
rankings[,3] <- unlist(apply(orderings,1,function(x){if("Third Place" %in% x){return(which(x=="Third Place"))}else{return(NA)}}))
rankings[,4] <- unlist(apply(orderings,1,function(x){if("Fourth Place" %in% x){return(which(x=="Fourth Place"))}else{return(NA)}}))
rankings[,5] <- unlist(apply(orderings,1,function(x){if("Fifth Place" %in% x){return(which(x=="Fifth Place"))}else{return(NA)}}))
rankings[,6] <- unlist(apply(orderings,1,function(x){if("Sixth Place" %in% x){return(which(x=="Sixth Place"))}else{return(NA)}}))
assignments <- matrix(data=TRUE,nrow=nrow(ratings),ncol=ncol(ratings))
for(paper in paper_IDS){
 judge_COIS <- unlist(raw_assignments[which(raw_assignments$RandomPropID == paper),c(4)])
 judge_COIS <- unlist(strsplit(judge_COIS,", "))
 if(length(judge_COIS)>0 & all(!is.na(judge_COIS))){
 for(judge in judge_COIS){
   assignments[which(judge_IDS == judge),which(paper_IDS == paper)] <- FALSE
 }}
}
assignments[1,setdiff(1:25,21)] <- FALSE
assignments[2,setdiff(1:25,24)] <- FALSE
assignments[16,setdiff(1:25,21)] <- FALSE
attr(rankings,"assignments") <- assignments
M <- 40
rm(merged_data,orderings,raw_assignments,rawdata,j,judge,judge_COIS,judge_IDS,paper,paper_IDS,assignments)
save.image("Section 4.3 Panel Review/Data/AIBS_Clean.RData")
