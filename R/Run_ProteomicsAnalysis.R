

library(heatmaply)
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
set.seed(6646)

source("./R/functions.R")
demogs <- read_xlsx("./Data/mod_demographics.xlsx")

IDs <- read.csv("./Results/ID_UpvotedClusters.csv")
csf_df <- read.csv("./Data/20190925_CSF_afterQC_proteomic_exprs.csv")

#Reading in proteomic datasets
df_pro_csf <- read.csv("./data/v4_CSF_afterQC_BasicPheno-All.csv")
df_pro_csf$MAPID <- as.numeric(as.character(df_pro_csf$MAPID))
df_pro_csf <- df_pro_csf[, c("MAPID", "ExtIdentifier", "TimePoint",  "CDR.Score")]

labels_for_Fig4 <- read_xlsx("./Results/ProteinRoles_BrainRevision.xlsx")
labels_for_Fig4 <- labels_for_Fig4[1:19,] #THese are predefined labels to color the Figure 4 barplot

df_pro_csf$TimePoint <-as.Date(format(as.Date(as.character(df_pro_csf$TimePoint), "%m/%d/%y"), "20%y-%m-%d"))


csf_df <- gather(csf_df, ex_ID,measurement,  EXID40000001884654:EXID40000001735507)

csf_df <- merge(csf_df, df_pro_csf, by.x = "ex_ID", by.y = "ExtIdentifier", all = FALSE)

csf_df <- merge(csf_df, IDs, by.x = "MAPID", by.y = "ID", all = FALSE)

csf_df <- csf_df[!duplicated(csf_df),]
csf_df <- data.table(csf_df)[, measurement := mean(measurement, na.rm = TRUE), by = list(MAPID, TargetFullName, TimePoint)]
csf_df <- unique(csf_df)
csf_df <- spread(csf_df[, c("MAPID", "TargetFullName", "measurement", "TimePoint")], TargetFullName, measurement)
#csf_df <- unique(data.table(csf_df)[order(-TimePoint)], by="MAPID")
csf_df <- merge(csf_df, IDs, by.x = "MAPID", by.y = "ID", all = FALSE)

Proteome_counts <- data.table(csf_df)[, .N, by = MAPID]
write.csv(Proteome_counts, "./Results/NumProcessedProteomes.csv", row.names = FALSE)

imputed_medians=data.frame(apply(csf_df,2,f))

csf_df <- cbind(csf_df[,1:2], imputed_medians[,3:715], csf_df[,716])
colnames(csf_df)[716] <- "upvoted_cluster"

csf_df$upvoted_cluster <- as.factor(csf_df$upvoted_cluster)

csf_df <- merge(csf_df, demogs[, c("ID", "BIRTH", "GENDER")],
                by.x = "MAPID", by.y = "ID", all = FALSE) #2 = female
csf_df$Age <- as.numeric(csf_df$TimePoint - as.Date(csf_df$BIRTH, format = "%Y-%m-%d"))/365.25
csf_df$Sex <- ifelse(csf_df$GENDER == 2, "female", "male")


#saving a holdout set
#Doing it by ID, since that takes out the random effect
#Also stratifying for equal balance across classes
tmp <- unique(csf_df[, c("MAPID", "upvoted_cluster")])
train_set <- caret::createDataPartition(tmp$upvoted_cluster, p = 0.8, list = FALSE)
train_IDs <- tmp[train_set, "MAPID"]
holdout_IDs <- tmp[-train_set, "MAPID"]



#Running through pelora for 1 - 10 clusters for each combo
AD_vs_Emerging_set <- get_MultipleClustersAUC(data.frame(csf_df), train_IDs$MAPID, 10, c(1, 0, NA), FALSE)
AD_vs_No_set <- get_MultipleClustersAUC(data.frame(csf_df), train_IDs$MAPID, 10, c(1, NA, 0), FALSE)
Emerging_vs_No_set <- get_MultipleClustersAUC(data.frame(csf_df), train_IDs$MAPID, 10, c(NA, 1, 0), FALSE)

#Running through pelora for 1 - 10 clusters for each combo, rerunning with age and sex as covariates
AD_vs_Emerging_set_cov <- get_MultipleClustersAUC(data.frame(csf_df), train_IDs$MAPID, 10, c(1, 0, NA), TRUE)
AD_vs_No_set_cov <- get_MultipleClustersAUC(data.frame(csf_df), train_IDs$MAPID, 10, c(1, NA, 0), TRUE)
Emerging_vs_No_set_cov <- get_MultipleClustersAUC(data.frame(csf_df), train_IDs$MAPID, 10, c(NA, 1, 0), TRUE)


#Generate final forecast based on best AUC
plot(AD_vs_Emerging_set[[1]]$NumClust, AD_vs_Emerging_set[[1]]$AUC) #2 is best
plot(AD_vs_No_set[[1]]$NumClust, AD_vs_No_set[[1]]$AUC) #10 is best
plot(Emerging_vs_No_set[[1]]$NumClust, Emerging_vs_No_set[[1]]$AUC) #5 is best


#Generate final forecast based on best AUC, with covariates
plot(AD_vs_Emerging_set_cov[[1]]$NumClust, AD_vs_Emerging_set[[1]]$AUC) #2
plot(AD_vs_No_set_cov[[1]]$NumClust, AD_vs_No_set[[1]]$AUC) #9 is best
plot(Emerging_vs_No_set_cov[[1]]$NumClust, Emerging_vs_No_set[[1]]$AUC) #5 is best


forecast_AD_Emerging <- predict(AD_vs_Emerging_set[[2]][[2]], AD_vs_Emerging_set[[3]], type = "class")
probs_AD_Emerging <- predict(AD_vs_Emerging_set[[2]][[2]], AD_vs_Emerging_set[[3]], type = "prob")

y_holdout_AD_Emerging <- AD_vs_Emerging_set[[4]]
levels(y_holdout_AD_Emerging) <- c(1, 0, NA) #AD vs Emerging

forecast_AD_Emerging <- data.frame("forecast" = unlist(forecast_AD_Emerging),
                                   "probability" = unlist(probs_AD_Emerging),
                                   "actual" = as.numeric(as.character(y_holdout_AD_Emerging)))
colnames(forecast_AD_Emerging) <- c("forecast", "probability", "actual")
table(data.frame(forecast_AD_Emerging)$forecast, data.frame(forecast_AD_Emerging)$actual)

rocobj_AD_Emerging <- roc(forecast_AD_Emerging$actual, forecast_AD_Emerging$probability)
as.numeric(ci(rocobj_AD_Emerging))
Fig4B <-pROC::ggroc(rocobj_AD_Emerging) + ggtitle( "B. Positive vs. Intermediate") + theme_bw() + theme(axis.text=element_text(size=12)) +
  annotate(geom = "text", x = 0.25, y = 0.25, label = paste0("AUC: ", round(as.numeric(ci(rocobj_AD_Emerging))[2], 3)), size = 5)



#Generate final forecast based on best AUC

forecast_AD_no <- predict(AD_vs_No_set[[2]][[10]], AD_vs_No_set[[3]], type = "class")
probs_AD_no <- predict(AD_vs_No_set[[2]][[10]], AD_vs_No_set[[3]], type = "prob")

y_holdout_AD_no <- AD_vs_No_set[[4]]
levels(y_holdout_AD_no) <- c(1, NA, 0) #AD vs no

forecast_AD_no <- data.frame("forecast" = unlist(forecast_AD_no),
                             "probability" = unlist(probs_AD_no),
                             "actual" = as.numeric(as.character(y_holdout_AD_no)))
colnames(forecast_AD_no) <- c("forecast", "probability", "actual")
table(data.frame(forecast_AD_no)$forecast, data.frame(forecast_AD_no)$actual)

rocobj_AD_no <- roc(forecast_AD_no$actual, forecast_AD_no$probability)
as.numeric(ci(rocobj_AD_no))
Fig4C <- pROC::ggroc(rocobj_AD_no) + ggtitle( "C. Positive vs. Negative") + theme_bw() + theme(axis.text=element_text(size=12)) +
  annotate(geom = "text", x = 0.25, y = 0.25, label = paste0("AUC: ", round(as.numeric(ci(rocobj_AD_no))[2], 3)), size = 5)




#Generate final forecast based on best AUC

forecast_Emerging_vs_No <- predict(Emerging_vs_No_set[[2]][[5]], Emerging_vs_No_set[[3]], type = "class")
probs_Emerging_vs_No <- predict(Emerging_vs_No_set[[2]][[5]], Emerging_vs_No_set[[3]], type = "prob")

y_holdout_Emerging_vs_No_set <- Emerging_vs_No_set[[4]]
levels(y_holdout_Emerging_vs_No_set) <- c(NA, 1, 0) #AD vs Emerging

forecast_Emerging_vs_No_set <- data.frame("forecast" = unlist(forecast_Emerging_vs_No),
                                          "probability" = unlist(probs_Emerging_vs_No),
                                          "actual" = as.numeric(as.character(y_holdout_Emerging_vs_No_set)))
colnames(forecast_Emerging_vs_No_set) <- c("forecast", "probability", "actual")
table(data.frame(forecast_Emerging_vs_No_set)$forecast, data.frame(forecast_Emerging_vs_No_set)$actual)

rocobj_Emerging_vs_No_set <- roc(forecast_Emerging_vs_No_set$actual, forecast_Emerging_vs_No_set$probability)
as.numeric(ci(rocobj_Emerging_vs_No_set))
Fig4A <- pROC::ggroc(rocobj_Emerging_vs_No_set) + ggtitle( "A. Intermediate vs. Negative") + theme_bw() + theme(axis.text=element_text(size=12)) + 
  annotate(geom = "text", x = 0.25, y = 0.25, label = paste0("AUC: ", round(as.numeric(ci(rocobj_Emerging_vs_No_set))[2], 3)), size = 5)

#Can tell AD vs no and Emerging vs No, but not really AD vs. Emerging, which is consistent with original submission

#Figuring out genes important for classification
importantGenes_Emerging_No <- get_importantGenes(Emerging_vs_No_set[[2]][[5]]$genes, Emerging_vs_No_set[[2]][[5]]$crit, Emerging_vs_No_set[[2]][[5]]$gene.names)
importantGenes_AD_No <- get_importantGenes(AD_vs_No_set[[2]][[10]]$genes, AD_vs_No_set[[2]][[10]]$crit, AD_vs_No_set[[2]][[10]]$gene.names)
importantGenes_AD_Emerging <- get_importantGenes(AD_vs_Emerging_set[[2]][[2]]$genes, AD_vs_Emerging_set[[2]][[2]]$crit, AD_vs_Emerging_set[[2]][[2]]$gene.names)

importantGenes_Emerging_No <- merge(importantGenes_Emerging_No, labels_for_Fig4[, c("Default_Name", "Nice_Name", "Category")],
                                    by.x = "gene", by.y = "Default_Name", all.x = TRUE, all.y = FALSE)
importantGenes_Emerging_No <- importantGenes_Emerging_No[with(importantGenes_Emerging_No, order(-score)),]
importantGenes_Emerging_No$Category <- as.factor(importantGenes_Emerging_No$Category)

importantGenes_AD_No <-  merge(importantGenes_AD_No, labels_for_Fig4[, c("Default_Name", "Nice_Name", "Category")],
                               by.x = "gene", by.y = "Default_Name", all.x = TRUE, all.y = FALSE)
importantGenes_AD_No <- importantGenes_AD_No[with(importantGenes_AD_No, order(-score)),]
importantGenes_AD_No$Category <- as.factor(importantGenes_AD_No$Category)

Fig4D <- 
  ggplot(importantGenes_Emerging_No[1:10,], aes(x=reorder(Nice_Name, score), y = score, label=score)) + 
  geom_bar(stat='identity', width=.9, aes(fill = Category))  +
  coord_flip() + ylab("Log Likelihood Criterion") + xlab("") +
  theme_bw() + ggtitle("D. Intermediate vs. Negative") +
  theme(legend.position = "bottom", axis.text=element_text(size=12)) +
  scale_fill_manual(values = c("#fde725", "#5ec962", "#440154", "#21918c"))

Fig4E <- ggplot(importantGenes_AD_No[1:10,], aes(x=reorder(Nice_Name, score), y = score, label=score)) + 
  geom_bar(stat='identity', width=.9, aes(fill= Category))  +
  coord_flip() + ylab("Log Likelihood Criterion") + xlab("") +
  theme_bw()+ ggtitle("E. Positive vs. Negative")+
  theme(legend.position = "bottom", axis.text=element_text(size=12)) +
  scale_fill_manual(values = c("#fde725", "#5ec962", "#3b528b", "#440154", "#21918c"))

lay <- rbind(c(1, 1, 2, 2, 3, 3), c(4, 4, 4, 5, 5, 5), c(4, 4, 4, 5, 5, 5), c(4, 4, 4, 5, 5, 5))
Cairo(file="./Results/Figure4.png",
      type="png",
      units="px", 
      width=1800, 
      height=1000, 
      pointsize=12, 
      dpi="auto")

grid.arrange(Fig4A, Fig4B, Fig4C, Fig4D, Fig4E, layout_matrix = lay)
dev.off()
#Next Step would be to fill the bars based on function

#Here I'm figuring out what the AUC is with age + sex for each model
forecast_AD_Emerging_cov <- predict(AD_vs_Emerging_set_cov[[2]][[2]], AD_vs_Emerging_set_cov[[3]], type = "class")
probs_AD_Emerging_cov <- predict(AD_vs_Emerging_set_cov[[2]][[2]], AD_vs_Emerging_set_cov[[3]], type = "prob")

y_holdout_AD_Emerging_cov <- AD_vs_Emerging_set_cov[[4]]
levels(y_holdout_AD_Emerging_cov) <- c(1, 0, NA) #AD vs Emerging

forecast_AD_Emerging_cov <- data.frame("forecast" = unlist(forecast_AD_Emerging_cov),
                                       "probability" = unlist(probs_AD_Emerging_cov),
                                       "actual" = as.numeric(as.character(y_holdout_AD_Emerging_cov)))
colnames(forecast_AD_Emerging_cov) <- c("forecast", "probability", "actual")
table(data.frame(forecast_AD_Emerging_cov)$forecast, data.frame(forecast_AD_Emerging_cov)$actual)

rocobj_AD_Emerging_cov <- roc(forecast_AD_Emerging_cov$actual, forecast_AD_Emerging_cov$probability)
as.numeric(ci(rocobj_AD_Emerging_cov)) #AUC = 0.625, 95%CI: 0.179, 1.000


forecast_AD_no_cov <- predict(AD_vs_No_set_cov[[2]][[9]], AD_vs_No_set_cov[[3]], AD_vs_No_set_cov[[5]],type = "class")
probs_AD_no_cov <- predict(AD_vs_No_set_cov[[2]][[9]], AD_vs_No_set_cov[[3]], AD_vs_No_set_cov[[5]], type = "prob")

y_holdout_AD_no_cov <- AD_vs_No_set_cov[[4]]
levels(y_holdout_AD_no_cov) <- c(1, NA, 0) #AD vs no

forecast_AD_no_cov <- data.frame("forecast" = unlist(forecast_AD_no_cov),
                                 "probability" = unlist(probs_AD_no_cov),
                                 "actual" = as.numeric(as.character(y_holdout_AD_no_cov)))
colnames(forecast_AD_no_cov) <- c("forecast", "probability", "actual")
table(data.frame(forecast_AD_no_cov)$forecast, data.frame(forecast_AD_no_cov)$actual)

rocobj_AD_no_cov <- roc(forecast_AD_no_cov$actual, forecast_AD_no_cov$probability)
as.numeric(ci(rocobj_AD_no_cov)) #AUC = 0.811, 95%CI = 0.669, 0.954


forecast_Emerging_vs_No_cov <- predict(Emerging_vs_No_set_cov [[2]][[5]], Emerging_vs_No_set_cov[[3]], type = "class")
probs_Emerging_vs_No_cov  <- predict(Emerging_vs_No_set_cov [[2]][[5]], Emerging_vs_No_set_cov[[3]], type = "prob")

y_holdout_Emerging_vs_No_set_cov  <- Emerging_vs_No_set_cov [[4]]
levels(y_holdout_Emerging_vs_No_set_cov ) <- c(NA, 1, 0) #AD vs Emerging

forecast_Emerging_vs_No_set_cov  <- data.frame("forecast" = unlist(forecast_Emerging_vs_No_cov ),
                                               "probability" = unlist(probs_Emerging_vs_No_cov ),
                                               "actual" = as.numeric(as.character(y_holdout_Emerging_vs_No_set_cov )))
colnames(forecast_Emerging_vs_No_set_cov ) <- c("forecast", "probability", "actual")
table(data.frame(forecast_Emerging_vs_No_set_cov )$forecast, data.frame(forecast_Emerging_vs_No_set_cov )$actual)

rocobj_Emerging_vs_No_set_cov  <- roc(forecast_Emerging_vs_No_set_cov $actual, forecast_Emerging_vs_No_set_cov $probability)
as.numeric(ci(rocobj_Emerging_vs_No_set_cov )) #0.865, 95%CI = 0.640, 1.000



#Applying lasso regression to try to assess predictive power of age & sex alone




Int_Pos_rocObj <- get_LassoAUC(c(1, 0, NA), csf_df, holdout_IDs) #Int vs Pos
Int_Neg_rocObj <- get_LassoAUC(c(1, NA, 0), csf_df, holdout_IDs) #Int vs Neg
Pos_Neg_rocObj <- get_LassoAUC(c(NA, 1, 0), csf_df, holdout_IDs) #Pos vs Neg

as.numeric(ci(Int_Pos_rocObj))
as.numeric(ci(Int_Neg_rocObj))
as.numeric(ci(Pos_Neg_rocObj))
