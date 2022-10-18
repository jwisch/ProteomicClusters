library(data.table)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(tidymodels)
library(xgboost)
library(pROC)
library(ggROC)
library(caret)
library(supclust)
library(imputeMissings)
library(glmnet)
library(heatmaply)
library(ppsr) #predictive power score

#Bonkers. I can classify pretty well the difference between Emerging and AD Pathology (AUC = 0.76),
#I can classify very well the difference between AD Pathology and healthy controls (AUC > 0.8)
#I can't classify unhealthy (Emerging + AD) vs healthy (AUC = 0.66)


IDs <- read.csv("./Data/Clustered_Participants.csv")
csf_df <- read.csv("./Data/20190925_CSF_afterQC_proteomic_exprs.csv")

#Reading in proteomic datasets
df_pro_csf <- read.csv("./data/v4_CSF_afterQC_BasicPheno-All.csv")
df_pro_csf$MAPID <- as.numeric(as.character(df_pro_csf$MAPID))
df_pro_csf <- df_pro_csf[, c("MAPID", "ExtIdentifier", "TimePoint",  "CDR.Score")]


df_pro_csf$TimePoint <-as.Date(format(as.Date(as.character(df_pro_csf$TimePoint), "%m/%d/%y"), "20%y-%m-%d"))


csf_df <- gather(csf_df, ex_ID,measurement,  EXID40000001884654:EXID40000001735507)

csf_df <- merge(csf_df, df_pro_csf, by.x = "ex_ID", by.y = "ExtIdentifier", all = FALSE)

csf_df <- merge(csf_df, IDs, by.x = "MAPID", by.y = "ID", all = FALSE)

csf_check <- csf_df[!duplicated(csf_df[, c("MAPID", "upvoted_cluster")]), c("MAPID", "upvoted_cluster")]
table(IDs$upvoted_cluster)
table(csf_check$upvoted_cluster)

csf_df <- csf_df[!duplicated(csf_df),]
csf_df <- data.table(csf_df)[, measurement := mean(measurement, na.rm = TRUE), by = list(MAPID, TargetFullName, TimePoint)]
csf_df <- unique(csf_df)
csf_df <- spread(csf_df[, c("MAPID", "TargetFullName", "measurement", "TimePoint",
                            "upvoted_cluster")], TargetFullName, measurement)

table(data.frame(data.table(csf_df)[, .N, by = MAPID])$N)

csf_df <- unique(data.table(csf_df)[order(-TimePoint)], by="MAPID")

#median imputation for missingness
values <- compute(csf_df[,4:716])
csf_df[,4:716] <- impute(csf_df[,4:716], method = "median", object = values)

csf_df$upvoted_cluster <- as.factor(csf_df$upvoted_cluster)
levels(csf_df$upvoted_cluster) <- c(NA, 1, 0) #classifying both ad and emerging as 1's, healthy as 0's
csf_df$upvoted_cluster <- as.numeric(as.character(csf_df$upvoted_cluster))
csf_df <- csf_df[!is.na(csf_df$upvoted_cluster),]
NumClust <- 10


####################################################################################################
##Single Run
####################################################################################################
IDs <- read.csv("./Data/Clustered_Participants.csv")
csf_df <- read.csv("./Data/20190925_CSF_afterQC_proteomic_exprs.csv")

#Reading in proteomic datasets
df_pro_csf <- read.csv("./Data/v4_CSF_afterQC_BasicPheno-All.csv")
df_pro_csf$MAPID <- as.numeric(as.character(df_pro_csf$MAPID))
df_pro_csf <- df_pro_csf[, c("MAPID", "ExtIdentifier", "TimePoint",  "CDR.Score")]


df_pro_csf$TimePoint <-as.Date(format(as.Date(as.character(df_pro_csf$TimePoint), "%m/%d/%y"), "20%y-%m-%d"))


csf_df <- gather(csf_df, ex_ID,measurement,  EXID40000001884654:EXID40000001735507)

csf_df <- merge(csf_df, df_pro_csf, by.x = "ex_ID", by.y = "ExtIdentifier", all = FALSE)

csf_df <- merge(csf_df, IDs, by.x = "MAPID", by.y = "ID", all = FALSE)

csf_check <- IDs[!duplicated(IDs[, c("ID", "upvoted_cluster")]), c("ID", "upvoted_cluster")]


csf_df <- csf_df[!duplicated(csf_df),]
csf_df <- data.table(csf_df)[, measurement := mean(measurement, na.rm = TRUE), by = list(MAPID, TargetFullName, TimePoint)]
csf_df <- unique(csf_df)
csf_df <- spread(csf_df[, c("MAPID", "TargetFullName", "measurement", "TimePoint")], TargetFullName, measurement)
csf_df_1 <- unique(data.table(csf_df)[order(-TimePoint)], by="MAPID")
csf_df <- merge(csf_df, csf_check, by.x = "MAPID", by.y = "ID", all = FALSE)

####################################################################################################
##################Predictive Power Score###########################################
####################################################################################################
# csf_pps <- csf_df
# names(csf_pps) <- c("MAPID", "TimePoint", paste0("x", 0:(ncol(csf_pps)-4)), "upvoted_cluster")
# protein_names <- names(csf_pps[, 3:715])
# csf_pps$upvoted_cluster <- as.factor(csf_pps$upvoted_cluster)
# pps_df <- data.frame("Protein" = rep(NA, length(protein_names)),
#                      "tree_pps" = rep(NA, length(protein_names)),
#                      "glm_pps" = rep(NA, length(protein_names)))
# for(i in 1:length(protein_names)){
#   pps_df$Protein[i] <- names(csf_df)[i + 2]
#   pps_df$tree_pps[i] <- score(csf_pps, x = protein_names[i], y = 'upvoted_cluster', algorithm = 'tree')[['pps']]
#   pps_df$glm_pps[i] <- score(csf_pps, x = protein_names[i], y = 'upvoted_cluster', algorithm = 'glm')[['pps']]
# }
# 
# write.csv(pps_df, "./Data/PredictivePowerScores.csv", row.names = FALSE)

pps_df <- read.csv("./Data/PredictivePowerScores.csv")
pps_df <- pps_df[with(pps_df, order(-tree_pps)),]

  p1 <- ggplot(pps_df[1:10,], aes(x = reorder(Protein, tree_pps), y = tree_pps)) + geom_col() + 
    theme_bw() + xlab("") + ylab("Predictive Power Score - Tree Based Modeling") + coord_flip() 
  
  pps_df <- pps_df[with(pps_df, order(-glm_pps)),]
  p2 <- ggplot(pps_df[1:10,], aes(x = reorder(Protein, glm_pps), y = glm_pps)) + geom_col() + 
    theme_bw() + xlab("") + ylab("Predictive Power Score - Linear Modeling") + coord_flip() 

grid.arrange(p1, p2, nrow = 1, top = "Predictive Power Score - Top 10 Proteins")
####################################################################################################
###################Getting correlation map###########################################
####################################################################################################

library(dendextend)
library(RColorBrewer)
# https://bio723-class.github.io/Bio723-book/clustering-in-r.html
corr_results <- cor(as.matrix(csf_df[, 3 : 715]), use = "pairwise.complete.obs", method = "spearman")
#corrplot(corr_results, method = "square")



x <- corr_results

row_dend  <- x %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 9) %>% 
  #set("branches_lwd", c(1, 3)) %>%
  ladderize
# rotate_DendSer(ser_weight = dist(x))
col_dend  <- x %>% 
  t %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 9) %>% 
 # set("branches_lwd", c(1, 2)) %>%
  ladderize
#    rotate_DendSer(ser_weight = dist(t(x)))

ggheatmap(corr_results, fontsize_row = 4, fontsize_col = 4,
          Rowv = row_dend, Colv = col_dend)




clusters <- cutree(col_dend, k=15) #look at N biggest clusters
table(clusters)
plot(color_branches(col_dend, k=9),leaflab="none")
clusters.df <- data.frame(protein = names(clusters), cluster = clusters) #making it to a thing I can read

#Generate heat map from cluster
cluster_select <- clusters.df[clusters.df$cluster == 14,]

  ggheatmap(corr_results[(colnames(corr_results) %in% cluster_select$protein &
                           rownames(corr_results) %in% cluster_select$protein),])
####################################################################################################


#median imputation for missingness
values <- compute(csf_df[,3:715])
csf_df[,3:715] <- impute(csf_df[,3:715], method = "median", object = values)

csf_df$upvoted_cluster <- as.factor(csf_df$upvoted_cluster)
# csf_df$upvoted_cluster <- as.numeric(as.character(csf_df$upvoted_cluster))
# csf_df <- csf_df[!is.na(csf_df$upvoted_cluster),]
NumClust <- 10

#need to run below for each comparison
csf_analysis <- csf_df
levels(csf_analysis$upvoted_cluster) <- c(1, NA, 0) #classifying AD = 1, healthy = 0
csf_analysis <- csf_analysis[!is.na(csf_analysis$upvoted_cluster)]
csf_analysis <- csf_analysis[with(csf_analysis, order(upvoted_cluster))]

set.seed(6646)
mod1 <- pelora(as.matrix(csf_analysis[,3:715]), as.numeric(as.character(csf_analysis$upvoted_cluster)), noc = NumClust,
               standardize = TRUE)
 
tmp <- data.frame(mod1$values)
tmp$Classification <- as.factor(csf_analysis$upvoted_cluster)
X_names <- paste0("X", seq(from = 1, to = NumClust))
#mod_tmp <- glm(as.formula(paste("Classification ~ ", paste(X_names, collapse= "+"))), family = "binomial", data = tmp)

cvfit <- cv.glmnet(as.matrix(tmp[, X_names]), tmp$Classification, family = "binomial", type.measure = "class", nfolds = 10)
# plot(cvfit)
retained_nums <- data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
retained_nums$keeps <- row.names(retained_nums)
retained_nums <- retained_nums[!(retained_nums$X1 == 0), "keeps"]
retained_nums <- retained_nums[2:length(retained_nums)]
retained_nums <- tidyr::extract_numeric(retained_nums)

test_matrix <- predict(mod1, csf_analysis[,3:715])
colnames(test_matrix) <- X_names

standardized_test_matrix <- standardize.genes(test_matrix[, X_names])

results <- data.frame("Probability" = predict(cvfit, newx = standardized_test_matrix$x, s = "lambda.min", type = "response"),
                      "Reality" = csf_analysis$upvoted_cluster)
results$Prediction <- ifelse(results$X1 < 0.5, 1, 0)
roc_obj <- roc(results$Reality, results$X1)

plot(roc_obj)

print(auc(roc_obj))

tmp_list <- list()


for(j in 1:length(retained_nums)){
  tmp_list[[j]] <- mod1$gene.names[mod1$genes[[retained_nums[j]]]]}

print("First Cluster of Proteomics used for classification:")
tmp_list[[1]]

print("Second Cluster of Proteomics used for classification:")
tmp_list[[2]]

print("Third Cluster of Proteomics used for classification:")
tmp_list[[3]]

protein_names <- unlist(tmp_list)
protein_names <- c( unique(protein_names))
protein_names <- gsub(" ", ".", protein_names)
protein_names <- gsub("-", ".", protein_names)
protein_names <- gsub("/", ".", protein_names)
protein_names <- gsub(",", ".", protein_names)
protein_names <- gsub("(", ".", protein_names, fixed = TRUE)
protein_names <- gsub(")", ".", protein_names)
csf_analysis <- data.frame(csf_analysis)
colnames(csf_analysis)[3:8] <- c("14.3.3.protein.beta.alpha", "14.3.3.protein.epsilon", "14.3.3.protein.family",    
                           "14.3.3.protein.sigma",  "14.3.3.protein.theta","14.3.3.protein.zeta.delta")

heatmap_df <-as.matrix(csf_analysis[, protein_names])
heatmap_df_scaled <- normalize(heatmap_df)
heatmaply(heatmap_df_scaled)
levels(csf_analysis$upvoted_cluster) <- c("AD", "No Path")
#https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
ggheatmap(heatmap_df_scaled, row_side_colors = csf_analysis[, "upvoted_cluster"], 
        Rowv = NULL)


####################################################################################################
####################################################################################################
####################################################################################################







# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2002-3-12-research0069
auc_list <- list()
results_list <- list()
grouped_genes <- list()
seed_list <- list()
for(i in 1:10){
  RandSeed <- sample(1:10000, 1)
  set.seed(RandSeed)
  seed_list[[i]] <- RandSeed
  csf_split <- initial_split(csf_df, strata = "upvoted_cluster", prop = 0.8)
  csf_train <- training(csf_split)
  csf_test <- testing(csf_split)
  # https://genomebiology.biomedcentral.com/articles/10.1186/gb-2002-3-12-research0069
  mod1 <- pelora(as.matrix(csf_train[,4:716]), as.numeric(csf_train$upvoted_cluster), noc = NumClust,
                 standardize = TRUE)
  # summary(mod1)
  # plot(mod1)
  # fitted(mod1)
  
  tmp <- data.frame(mod1$values)
  tmp$Classification <- as.factor(csf_train$upvoted_cluster)
  X_names <- paste0("X", seq(from = 1, to = NumClust))
  #mod_tmp <- glm(as.formula(paste("Classification ~ ", paste(X_names, collapse= "+"))), family = "binomial", data = tmp)
  
  cvfit <- cv.glmnet(as.matrix(tmp[, X_names]), tmp$Classification, family = "binomial", type.measure = "class", nfolds = 10)
  # plot(cvfit)
  retained_nums <- data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
  retained_nums$keeps <- row.names(retained_nums)
  retained_nums <- retained_nums[!(retained_nums$X1 == 0), "keeps"]
  retained_nums <- retained_nums[2:length(retained_nums)]
  retained_nums <- tidyr::extract_numeric(retained_nums)
  
  test_matrix <- predict(mod1, csf_test[,4:716])
  colnames(test_matrix) <- X_names
  
  standardized_test_matrix <- standardize.genes(test_matrix[, X_names])
  
  results <- data.frame("Probability" = predict(cvfit, newx = standardized_test_matrix$x, s = "lambda.min", type = "response"),
                        "Reality" = csf_test$upvoted_cluster)
  results$Prediction <- ifelse(results$X1 > 0.5, 1, 0)
  results_list[[i]] <- results
  tmp_list <- list()
  for(j in 1:length(retained_nums)){
    tmp_list[[j]] <- mod1$gene.names[mod1$genes[[retained_nums[j]]]]}
  grouped_genes[[i]] <- list(tmp_list)
  # mod1$gene.names[mod1$genes[[4]]]
  # mod1$gene.names[mod1$genes[[9]]]
  # mod1$gene.names[mod1$genes[[5]]]
  # mod1$gene.names[mod1$genes[[10]]]
  
  roc_obj <- roc(results$Reality, results$X1)
  auc_list[[i]] <- as.numeric(auc(roc_obj)) 
}

saveRDS(seed_list, "./Results/SuccessfulProteomicsClassificationofEmergingvsNonAD2_seeds.RDS")
saveRDS(grouped_genes, "./Results/SuccessfulProteomicsClassificationofEmergingvsNonAD2_genes.RDS")
saveRDS(auc_list, "./Results/SuccessfulProteomicsClassificationofEmergingvsNonAD2_auc.RDS")
saveRDS(results_list, "./Results/SuccessfulProteomicsClassificationofEmergingvsNonAD2_results.RDS")


mean(unlist(readRDS("./Results/SuccessfulProteomicsClassificationofEmergingvsAD_auc.RDS")))
mean(unlist(readRDS("./Results/SuccessfulProteomicsClassificationofADvsNonAD_auc.RDS")))
mean(unlist(readRDS("./Results/SuccessfulProteomicsClassificationofEmergingvsNonAD2_auc.RDS")))
mean(unlist(readRDS("./Results/UnSuccessfulProteomicsClassificationofUnhealthyvsHealthy_auc.RDS")))
