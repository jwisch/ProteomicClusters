library(tidyr)
library(data.table)
library(imputeMissings)
library(ppsr)
library(pROC)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggcorrplot)
library('Cairo')
CairoWin() #dealing with aliasing


clusters <- read.csv("./Data/csf_processed_for_analysis.csv")
clusters <- clusters[, c("ID", "class")]
clusters <- unique(clusters)

csf_df <- read.csv("./Data/20190925_CSF_afterQC_proteomic_exprs.csv")

# "important" proteins come from Fig 4 in manuscript
labels_for_Fig4 <- read_xlsx("./Results/ProteinRoles_BrainRevision.xlsx")
labels_for_Fig4 <- labels_for_Fig4[1:19,] #THese are predefined labels to color the Figure 4 barplot

#Reading in proteomic datasets
df_pro_csf <- read.csv("./Data/v4_CSF_afterQC_BasicPheno-All.csv")
df_pro_csf$MAPID <- as.numeric(as.character(df_pro_csf$MAPID))
df_pro_csf <- df_pro_csf[, c("MAPID", "ExtIdentifier", "TimePoint",  "CDR.Score")]


df_pro_csf$TimePoint <-as.Date(format(as.Date(as.character(df_pro_csf$TimePoint), "%m/%d/%y"), "20%y-%m-%d"))


csf_df <- gather(csf_df, ex_ID,measurement,  EXID40000001884654:EXID40000001735507)

csf_df <- merge(csf_df, df_pro_csf, by.x = "ex_ID", by.y = "ExtIdentifier", all = FALSE)

csf_df <- merge(csf_df, clusters, by.x = "MAPID", by.y = "ID", all = FALSE)

csf_df <- csf_df[!duplicated(csf_df),]
csf_df <- data.table(csf_df)[, measurement := mean(measurement, na.rm = TRUE), by = list(MAPID, TargetFullName, TimePoint)]
csf_df <- unique(csf_df)
csf_df <- spread(csf_df[, c("MAPID", "TargetFullName", "measurement", "TimePoint")], TargetFullName, measurement)
#csf_df <- unique(data.table(csf_df)[order(-TimePoint)], by="MAPID")
csf_df <- merge(csf_df, clusters, by.x = "MAPID", by.y = "ID", all = FALSE)

#median imputation for missingness
#median imputation for missingness
f=function(x){
  x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
  x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
  x #display the column
}
imputed_medians=data.frame(apply(csf_df,2,f))

csf_df <- cbind(csf_df[,1:2], imputed_medians[,3:715], csf_df[,716])
colnames(csf_df)[716] <- "upvoted_cluster"

csf_df$upvoted_cluster <- as.factor(csf_df$upvoted_cluster)
names(csf_df) <- make.names(names(csf_df), unique = TRUE)

csf_df <- data.frame(csf_df)
csf_df[3:715] <- sapply(csf_df[3:715],as.numeric)




corr <- round(cor(csf_df[,names(csf_df) %in% labels_for_Fig4$Default_Name]), 3)
g <- ggcorrplot(corr, hc.order = TRUE) + 
  scale_fill_viridis(limits = c(-0.345, 1)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1, size = 12),
                                                                 axis.text.y = element_text(size = 12)) + xlab("") + ylab("") + labs(fill = "") +
  ggtitle("Correlation Matrix of Identified Genes")

ggsave(g, filename = 'SuppFig13.png', dpi = 300, type = 'cairo',
       width = 16, height = 16, units = 'in')

#Calculating individual predictive power score for each protein
ppsr_list <- list()
for(i in 3:715){
  ppsr_list[[i-2]] <- ppsr::score(csf_df, x = names(csf_df)[i], y = 'upvoted_cluster', algorithm = 'glm')[['pps']]
}

ppsr_scores <- data.frame("Protein" = names(csf_df)[3:715],
                          "PredictivePowerScore" = unlist(ppsr_list))


#Calculating AUC for each individual classification

csf_df <- merge(csf_df, clusters[, c("ID", "class")], by.x = "MAPID", by.y = "ID", all = FALSE)
AD_Emerging <- csf_df[!(csf_df$upvoted_cluster == "AD Biomarker Negative"),]
AD_Emerging$Class <- ifelse(AD_Emerging$upvoted_cluster == "AD Biomarker Positive", 1, 0)

AD_No <- csf_df[!(csf_df$upvoted_cluster == "Intermediate AD Biomarkers"),]
AD_No$Class <- ifelse(AD_No$upvoted_cluster == "AD Biomarker Positive", 1, 0)

Emerging_No <- csf_df[!(csf_df$upvoted_cluster == "AD Biomarker Positive"),]
Emerging_No$Class <- ifelse(Emerging_No$upvoted_cluster == "Intermediate AD Biomarkers", 1, 0)

get_AUC_and_P <- function(DF){

AD_Emerging_AUC <- list()
AD_Emerging_p <- list()
for(i in 3:715){
logit_formula <- as.formula(paste0("Class ~ ", names(DF)[i]))
logit_model <- glm(logit_formula, family = "binomial", data = DF)
AD_Emerging_p[[i - 2]]  <- summary(logit_model)$coefficients[2, 4]
forecast_comparison <- data.frame("Prediction" = predict(logit_model, DF, type = "response"),
                                  "Actual" = DF[,"Class"])
ROC <- roc(forecast_comparison$Actual, forecast_comparison$Prediction)
AD_Emerging_AUC[[i - 2]] <- auc(ROC)} 

comparison_DF <- data.frame("Protein" = names(DF)[3:715],
                            "AUC" = unlist(AD_Emerging_AUC),
                            "p_value" = unlist(AD_Emerging_p))

return(comparison_DF)}

AD_Emerging_DF <- get_AUC_and_P(AD_Emerging)

AD_No_DF <- get_AUC_and_P(AD_No)

Emerging_No_DF <- get_AUC_and_P(Emerging_No)


write.csv(ppsr_scores, "./Results/ppsr_individual_proteins.csv", row.names = FALSE)
write.csv(AD_Emerging_DF, "./Results/IndividualProteinAUC_ADvsEmerging.csv", row.names = FALSE)
write.csv(AD_No_DF, "./Results/IndividualProteinAUC_ADvsNo.csv", row.names = FALSE)
write.csv(Emerging_No_DF, "./Results/IndividualProteinAUC_EmergingvsNo.csv", row.names = FALSE)


################################################################################

ppsr_scores <- read.csv("./Results/ppsr_individual_proteins.csv")

#VIZ
ppsr_scores <- arrange(ppsr_scores, PredictivePowerScore)
ppsr_scores$Class <- as.factor(ifelse(ppsr_scores$Protein %in% labels_for_Fig4$Default_Name, 
                                         "Appears in pelora classification", 
                                         "Does not appear in pelora classification"))
ggplot(ppsr_scores[ppsr_scores$PredictivePowerScore > 0.245,], 
       aes(x = reorder(Protein, PredictivePowerScore), y = PredictivePowerScore, fill = Class)) +
  geom_col() + coord_flip() + theme_bw() + xlab(" ") + ylab("Predictive Power Score") + scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "bottom")


AD_No_DF$Class <- as.factor(ifelse(AD_No_DF$Protein %in% labels_for_Fig4$Default_Name, 
                                      "Appears in pelora classification", 
                                      "Does not appear in pelora classification"))
ggplot(AD_No_DF[AD_No_DF$AUC > 0.75,], 
       aes(x = reorder(Protein, AUC), y = AUC, fill = Class)) +
  geom_col() + coord_flip() + theme_bw() + xlab(" ") + ylab("AUC") + scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "bottom", legend.title = element_blank())


Emerging_No_DF$Class <- as.factor(ifelse(Emerging_No_DF$Protein %in% labels_for_Fig4$Default_Name, 
                                   "Appears in pelora classification", 
                                   "Does not appear in pelora classification"))
ggplot(Emerging_No_DF[Emerging_No_DF$AUC > 0.75,], 
       aes(x = reorder(Protein, AUC), y = AUC, fill = Class)) +
  geom_col() + coord_flip() + theme_bw() + xlab(" ") + ylab("AUC") + scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "bottom", legend.title = element_blank())




csf_df$class <- factor(csf_df$class.x, levels=c("Intermediate AD Biomarkers", "AD Biomarker Negative", "AD Biomarker Positive"))

labels_for_Fig4[labels_for_Fig4$Category == "Neurodegeneration","Default_Name"]


my_comparisons <- list( c("AD Biomarker Positive", "AD Biomarker Negative"), 
                        c("AD Biomarker Positive", "Intermediate AD Biomarkers"), 
                        c("Intermediate AD Biomarkers", "AD Biomarker Negative") )

p1 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "X14.3.3.protein.zeta.delta",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("14-3-3 Protein Zeta Delta") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p2 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Neuronal.growth.regulator.1",
                  color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Neuronal Growth Regulator 1") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p3 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "X14.3.3.protein.epsilon",
                  color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("14-3-3 Protein Epsilon") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p4 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Semaphorin.6A",
                  color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Semaphorin 6A") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p5 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "X14.3.3.protein.family",
                  color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("14-3-3 Protein Family") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p6 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Acidic.leucine.rich.nuclear.phosphoprotein.32.family.member.B",
                  color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Acidic Leucine Rich Nuclear Phosphoprotein 32 Family Member B") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)


labels_for_Fig4[labels_for_Fig4$Category == "BBB/Vascular","Default_Name"]

p1 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "SPARC.related.modular.calcium.binding.protein.1",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("SPARC Reulated Modular Calcium Binding Protein 1") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p2 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Vascular.endothelial.growth.factor.C",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Vascular Endothelial Growth Factor C") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p3 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Discoidin.domain.containing.receptor.2",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Discoidin Domain Containing Receptor 2") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p4 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Endostatin",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Endostatin") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p5 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "ADP.ribosyl.cyclase.cyclic.ADP.ribose.hydrolase.2",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("ADP Ribosyl Cyclase Cyclic ADP Ribose Hydrolase 2") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 2)

labels_for_Fig4[labels_for_Fig4$Category == "Immune Function","Default_Name"]

p1 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Antileukoproteinase",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Antileukoproteinase") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p2 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Immunoglobulin.G",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Immunoglobulin G") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p3 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Hepatitis.A.virus.cellular.receptor.2",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Hepatitis A Virus Cellular Receptor 2") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)


labels_for_Fig4[labels_for_Fig4$Category == "Other","Default_Name"]

p1 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Interleukin.20.receptor.subunit.alpha",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("IL-20") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p2 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Hepcidin",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Hepcidin") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p3 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Superoxide.dismutase..Mn...mitochondrial",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Superoxide Dismutase MN") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)

labels_for_Fig4[labels_for_Fig4$Category == "Inflamation","Default_Name"]

p1 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Thrombospondin.4",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Thrombospondin 4") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

p2 <- ggpubr::ggboxplot(csf_df, x = "class.x", y = "Growth.differentiation.factor.15",
                        color = "class.x", add = "jitter")+ 
  ggpubr::stat_compare_means(comparisons = my_comparisons) + xlab("") + ylab("Growth Differentiation Factor 15") +
  scale_colour_manual(values = c( "#21918c", "#fde725", "#440154")) + theme(legend.position = "none")

grid.arrange(p1, p2, nrow = 2, ncol = 1)


