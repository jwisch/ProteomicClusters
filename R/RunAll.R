

library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(lme4)
library(table1)
library(plotly)
library(data.table)
library(RobMixReg)
library(lcmm)
library(tidyr)
library(pROC)
library(ggROC)
library(caret)
library(supclust)
library(imputeMissings)
library(glmnet)
library(heatmaply)
library(ape)
library(survival)
library('Cairo')
library(tableone)
library(ggfortify)
library(tableone)
library(MonoInc)


#CairoWin() #dealing with aliasing

source("./R/functions.R")


#######################PARAMS############################
MinNumVisits <- 4
############################



#######################READ IN DATA############################
pet <- read_xlsx("./Data/mod_pet.xlsx")
mmse <- read_xlsx("./Data/mod_b4_cdr.xlsx")
demogs <- read_xlsx("./Data/mod_demographics.xlsx")
csf <- read_xlsx("./Data/mod_lumipulse.xlsx")
np <- read.csv("./Data/mod_psychometrics.csv")
mri <- read_xlsx("./Data/mod_mri.xlsx")
tau <- read_xlsx("./Data/mod_tau.xlsx")
wmh <- read_xlsx("./Data/mod_wmh.xlsx")
nfl <- read.csv("./Data/csf_nfl_20210301.csv")

csf_long <- setDT(csf)[, .N, by = ID]
csf_long <- csf_long[csf_long$N >= MinNumVisits,]

csf <- csf[, c("ID", "CSF_LP_DATE", "LUMIPULSE_CSF_AB42", "LUMIPULSE_CSF_AB40",
               "LUMIPULSE_CSF_pTau", "LUMIPULSE_CSF_tTau")]
csf$LUMIPULSE_CSF_AB42_AB40 <- csf$LUMIPULSE_CSF_AB42 / csf$LUMIPULSE_CSF_AB40
csf$CSF_LP_DATE <- as.Date(csf$CSF_LP_DATE, format = "%Y-%m-%d")
csf <- csf[csf$ID %in% csf_long$ID,]
########################################################

#####################BASIC DATA CLEANING##################
demogs$BIRTH <- as.Date(demogs$BIRTH, format = "%Y-%m-%d")
demogs <- demogs[, c("ID", "BIRTH", "GENDER", "EDUC", "apoe", "race2")]
demogs$apoe4 <- ifelse(demogs$apoe == 34 | demogs$apoe == 44, 1, 0)
demogs <- demogs[demogs$ID %in% unique(csf_long$ID),]
demogs <- merge(demogs, cluster_results, by = "ID")



csf <- merge(csf, demogs, by = "ID")
csf$Age_at_LP <- round(as.numeric(csf$CSF_LP_DATE - csf$BIRTH) / 365.25, 1)
nfl$LP_DATE <- as.Date(nfl$LP_DATE, format = "%m/%d/%Y")
csf <- merge(csf, nfl, by.x = c("ID", "CSF_LP_DATE"), by.y = c("ID", "LP_DATE"), all.x = TRUE, all.y = FALSE)

pet<- pet[!(pet$PUP_QC_Status == "Failed" | pet$PUP_QC_Status == "Quarantined"), ]
pet <- pet[pet$Tracer == "PIB",]
pet <- pet[, c("ID", "PET_Date", "pib_fsuvr_rsf_tot_cortmean")]
pet$PET_Date <- as.Date(pet$PET_Date, format = "%Y-%m-%d")
pet <- pet[!is.na(pet$pib_fsuvr_rsf_tot_cortmean),]
pet <- merge(pet, demogs, by = "ID", all = FALSE)
pet$Age_at_Scan <- as.numeric(pet$PET_Date - pet$BIRTH) / 365.25


mmse <- mmse[, c("ID", "testdate", "cdr", "MMSE")]
mmse$testdate <- as.Date(mmse$testdate, format = "%Y-%m-%d")
mmse <- mmse[!is.na(mmse$MMSE),]
mmse <- merge(mmse, demogs, by = "ID", all = FALSE)
mmse$Age_at_Visit <- as.numeric(mmse$testdate - mmse$BIRTH) / 365.25


np <- np[, c("ID", "psy_date", "srttotal", "lmdelay", "digsym", "MEMUNITS")]
np$psy_date <- as.Date(np$psy_date, format = "%d-%b-%y")
np <- merge(np, demogs, by = "ID")
np$Age_at_Visit <- as.numeric(np$psy_date - np$BIRTH) / 365.25


wmh <- wmh[, c("ID", "session_date", "WMH_volume")]
wmh$session_date <- as.Date(wmh$session_date, format = "%Y-%m-%d")
wmh <- merge(wmh, demogs, by = "ID", all = FALSE)
wmh$Age_at_Scan <- as.numeric(wmh$session_date - wmh$BIRTH) / 365.25


tau <- tau[, c("ID", "PET_Date", "Tauopathy")]
tau$PET_Date <- as.Date(tau$PET_Date, format = "%Y-%m-%d")
tau <- tau[!is.na(tau$Tauopathy),]
tau <- merge(tau, demogs, by = "ID", all = FALSE)
tau$Age_at_Scan <- as.numeric(tau$PET_Date - tau$BIRTH) / 365.25

mri <- mri[mri$FS_QC_Status == "Passed" | mri$FS_QC_Status == "Passed with edits",]
mri <- mri[,c("ID", "session_date", "CortSig_Thickness", "MR_TOTV_INTRACRANIAL", 
              "MR_LV_HIPPOCAMPUS", "MR_RV_HIPPOCAMPUS")]
mri$session_date <- as.Date(mri$session_date, format = "%Y-%m-%d")
mri <- merge(mri, demogs, by = "ID", all = FALSE)
mri$Age_at_Scan <- as.numeric(mri$session_date - mri$BIRTH) / 365.25
############################################################################


#######################UNSUPERVISED CLUSTERING############################
set.seed(2002)
gmm1_2 <- hlme(LUMIPULSE_CSF_pTau ~ LUMIPULSE_CSF_AB42, subject = "ID", random=~1+LUMIPULSE_CSF_AB42, ng = 1, data =
                 csf, verbose = FALSE)
gmm2_2 <- gridsearch(rep = 100, maxiter = 10, minit = gmm1_2,
                     hlme(LUMIPULSE_CSF_pTau ~ LUMIPULSE_CSF_AB42, subject = "ID", random=~1+LUMIPULSE_CSF_AB42,
                          ng = 2, data = csf, mixture = ~ LUMIPULSE_CSF_AB42,
                          nwg=T, verbose = FALSE))
gmm3_2 <- gridsearch(rep = 100, maxiter = 10, minit = gmm1_2,
                     hlme(LUMIPULSE_CSF_pTau ~ LUMIPULSE_CSF_AB42, subject = "ID", random=~1+LUMIPULSE_CSF_AB42,
                          ng = 3, data = csf, mixture = ~ LUMIPULSE_CSF_AB42,
                          nwg=T, verbose = FALSE))
gmm4_2 <- gridsearch(rep = 100, maxiter = 10, minit = gmm1_2,
                     hlme(LUMIPULSE_CSF_pTau ~ LUMIPULSE_CSF_AB42, subject = "ID", random=~1+LUMIPULSE_CSF_AB42,
                          ng = 4, data = csf, mixture = ~ LUMIPULSE_CSF_AB42,
                          nwg=T, verbose = FALSE))

#summarytable(gmm1_2, gmm2_2, gmm3_2, gmm4_2) #figure out which one is best - smallest BIC
########################################################

#########PLOTTING TRAJECTORIES OF CLASSES################
#Proposed by Reviewer 1...not sure if necessary
#Included in response to reviewers but not in manuscript
cluster_proj <- data.frame("LUMIPULSE_CSF_AB42" = seq(from = min(csf$LUMIPULSE_CSF_AB42),
                                                      to = max(csf$LUMIPULSE_CSF_AB42), length(csf$LUMIPULSE_CSF_AB42)))

plot(lcmm::predictY(gmm3_2, cluster_proj, var.time = "LUMIPULSE_CSF_AB42"))
########################################################


#########TESTING MONOTONICITY################
#Provided in response to reviewers but not in manuscript
csf_surv$pTau_AB42 <- csf_surv$LUMIPULSE_CSF_pTau / csf_surv$LUMIPULSE_CSF_AB42

get_proportion_monotonic(csf_surv, "ID", "LUMIPULSE_CSF_pTau", 'inc')
get_proportion_monotonic(csf_surv[csf_surv$class == "AD Biomarker Negative",], "ID", "LUMIPULSE_CSF_pTau", 'inc')
get_proportion_monotonic(csf_surv[csf_surv$class == "Intermediate AD Biomarkers",], "ID", "LUMIPULSE_CSF_pTau", 'inc')
get_proportion_monotonic(csf_surv[csf_surv$class == "AD Biomarker Positive",], "ID", "LUMIPULSE_CSF_pTau", 'inc')

get_proportion_monotonic(csf_surv, "ID", "LUMIPULSE_CSF_AB42_AB40", 'dec')
get_proportion_monotonic(csf_surv[csf_surv$class == "AD Biomarker Negative",], "ID", "LUMIPULSE_CSF_AB42_AB40", 'dec')
get_proportion_monotonic(csf_surv[csf_surv$class == "Intermediate AD Biomarkers",], "ID", "LUMIPULSE_CSF_AB42_AB40", 'dec')
get_proportion_monotonic(csf_surv[csf_surv$class == "AD Biomarker Positive",], "ID", "LUMIPULSE_CSF_AB42_AB40", 'dec')

get_proportion_monotonic(csf_surv, "ID", "pTau_AB42", 'inc')
get_proportion_monotonic(csf_surv[csf_surv$class == "AD Biomarker Negative",], "ID", "pTau_AB42", 'inc')
get_proportion_monotonic(csf_surv[csf_surv$class == "Intermediate AD Biomarkers",], "ID", "pTau_AB42", 'inc')
get_proportion_monotonic(csf_surv[csf_surv$class == "AD Biomarker Positive",], "ID", "pTau_AB42", 'inc')

########################################################




#######################SURVIVAL ANALYSIS############################
csf_surv <- merge(csf[, c("ID", "CSF_LP_DATE", "LUMIPULSE_CSF_AB42", 
                          "LUMIPULSE_CSF_pTau", "LUMIPULSE_CSF_AB42_AB40")], demogs[, c("ID", "BIRTH", "GENDER", "EDUC", "apoe")],
                  by = "ID", all.x = TRUE, all.y = FALSE)
csf_surv$Age <- as.numeric(csf_surv$CSF_LP_DATE - as.Date(csf_surv$BIRTH, format = "%Y-m-%d")) / 365.25

csf_surv$Apos <- ifelse(csf_surv$LUMIPULSE_CSF_AB42_AB40 < 0.0673, 1, 0)
csf_surv$Tpos <- ifelse(csf_surv$LUMIPULSE_CSF_pTau > 42.5, 1, 0)
csf_surv$CSF_LP_DATE <- as.Date(csf_surv$CSF_LP_DATE, format = "%Y-%m-%d")
csf_surv <- merge(csf_surv, gmm3_2$pprob[, c("ID", "class")], by = "ID")
pet$ID <- as.factor(pet$ID)

write.csv(gmm3_2$pprob[, c("ID", "class")], "./Results/ID_UpvotedClusters.csv", row.names = FALSE)

csf_surv$class <- as.factor(csf_surv$class)
levels(csf_surv$class) <- c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative")
csf_surv$class <- factor(csf_surv$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

#Amyloid
csf_surv_A <- get_modality_survival_df(csf_surv, "CSF_LP_DATE", "Apos")
km_obj_A <- get_run_survival_analysis(csf_surv_A)

#Tau
csf_surv_T <- get_modality_survival_df(csf_surv, "CSF_LP_DATE","Tpos")
km_obj_T <- get_run_survival_analysis(csf_surv_T)

#Neurodegeneration
mmse <- data.frame(mmse)
mmse$Npos <- ifelse(mmse$cdr > 0, 1, 0)
names(mmse)[names(mmse) == 'Age_at_Visit'] <- 'Age'
mmse <- merge(mmse, gmm3_2$pprob[, c("ID", "class")], by = "ID", all = FALSE)

mmse$class <- as.factor(mmse$class)
levels(mmse$class) <- c("Intermediate AD Biomarkers", "AD Biomarker Positive", "AD Biomarker Negative")
mmse$class <- factor(mmse$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))


csf_surv_N <- get_modality_survival_df(mmse, "testdate","Npos")
km_obj_N <- get_run_survival_analysis(csf_surv_N)

########################################################


#######################GAMM ANALYSIS############################
#https://astrostatistics.psu.edu/su07/R/html/mgcv/html/gamm.html
#https://stats.stackexchange.com/questions/536167/how-to-fit-a-repeated-measure-gam-model-with-mgcv
library(mgcv)
library(itsadug)
source("./gamm_hacks.R") #https://github.com/soskuthy/gamm_intro/blob/master/gamm_hacks.r


#Random smooths like: https://arxiv.org/pdf/1703.05339.pdf
csf_surv$ID <- as.factor(csf_surv$ID)

BaseModel <- bam(LUMIPULSE_CSF_pTau ~  s(LUMIPULSE_CSF_AB42, bs = "cr") +
                   s(LUMIPULSE_CSF_AB42, class, bs = "fs",  m = 1, k = 3),
                 data = csf_surv, method = "fREML")
plot_smooth(BaseModel, view="LUMIPULSE_CSF_AB42", plot_all="class",
            rug=F, rm.ranef=T)

plot(BaseModel)
summary.coefs(BaseModel)
plot_smooth.cont(BaseModel, view="LUMIPULSE_CSF_AB42", plot_all.c="class",
                 rug=F, rm.ranef=T)

df_forecast <- data.frame("LUMIPULSE_CSF_AB42" = c(seq(from = min(csf_surv[csf_surv$class == "AD Biomarker Negative",]$LUMIPULSE_CSF_AB42),
                                                       to = max(csf_surv[csf_surv$class == "AD Biomarker Negative",]$LUMIPULSE_CSF_AB42), 
                                                       length.out = 1000),
                                                   seq(from = min(csf_surv[csf_surv$class == "AD Biomarker Positive",]$LUMIPULSE_CSF_AB42),
                                                       to = max(csf_surv[csf_surv$class == "AD Biomarker Positive",]$LUMIPULSE_CSF_AB42), 
                                                       length.out = 1000),
                                                   seq(from = min(csf_surv[csf_surv$class == "Intermediate AD Biomarkers",]$LUMIPULSE_CSF_AB42),
                                                       to = max(csf_surv[csf_surv$class == "Intermediate AD Biomarkers",]$LUMIPULSE_CSF_AB42), 
                                                       length.out = 1000)),
                          "class" = c(rep("AD Biomarker Negative", 1000), rep("AD Biomarker Positive", 1000),
                                      rep("Intermediate AD Biomarkers", 1000)))

df_forecast$class <- as.factor(df_forecast$class)
df_forecast$class <- factor(df_forecast$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

forecast <- predict(BaseModel, df_forecast, type = "link", se.fit = TRUE)
df_forecast$forecast <- forecast$fit
df_forecast$upper <- forecast$fit + 1.96 * forecast$se.fit
df_forecast$lower <- forecast$fit - 1.96 * forecast$se.fit
Fig1A <- ggplot(df_forecast, aes(x = LUMIPULSE_CSF_AB42, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() + #geom_vline(xintercept = 0.0673, linetype = "dashed", colour = "grey67") +
  geom_hline(yintercept = 42.5, linetype = "dashed", colour = "grey67") +
  geom_vline(xintercept = 1098, linetype = "dashed", colour = "grey67") +
  ylab("CSF pTau") + xlab("CSF AB42") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster") + scale_x_reverse()  +
  geom_point(csf_surv, mapping=aes(x = LUMIPULSE_CSF_AB42, y = LUMIPULSE_CSF_pTau, colour = class, fill = class)) +  
  geom_line(csf_surv, mapping=aes(x = LUMIPULSE_CSF_AB42, y = LUMIPULSE_CSF_pTau, colour = class, fill = class, group = ID), alpha = 0.7) +
  theme(legend.position = "bottom", legend.text = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.title = element_text(size = 20), 
        legend.title = element_blank(), title = element_text(size = 24))  + ggtitle(label = "A.")



RatioModel <- bam(LUMIPULSE_CSF_pTau ~  s(LUMIPULSE_CSF_AB42_AB40, bs = "cr") +
                    s(LUMIPULSE_CSF_AB42_AB40, class, bs = "fs",  m = 1, k = 3),
                  data = csf_surv, method = "fREML")


df_forecast_ratio <- data.frame("LUMIPULSE_CSF_AB42_AB40" = c(seq(from = min(csf_surv[csf_surv$class == "AD Biomarker Negative",]$LUMIPULSE_CSF_AB42_AB40),
                                                                  to = max(csf_surv[csf_surv$class == "AD Biomarker Negative",]$LUMIPULSE_CSF_AB42_AB40), 
                                                                  length.out = 1000),
                                                              seq(from = min(csf_surv[csf_surv$class == "AD Biomarker Positive",]$LUMIPULSE_CSF_AB42_AB40),
                                                                  to = max(csf_surv[csf_surv$class == "AD Biomarker Positive",]$LUMIPULSE_CSF_AB42_AB40), 
                                                                  length.out = 1000),
                                                              seq(from = min(csf_surv[csf_surv$class == "Intermediate AD Biomarkers",]$LUMIPULSE_CSF_AB42_AB40),
                                                                  to = max(csf_surv[csf_surv$class == "Intermediate AD Biomarkers",]$LUMIPULSE_CSF_AB42_AB40), 
                                                                  length.out = 1000)),
                                "class" = c(rep("AD Biomarker Negative", 1000), rep("AD Biomarker Positive", 1000),
                                            rep("Intermediate AD Biomarkers", 1000)))

df_forecast_ratio$class <- as.factor(df_forecast_ratio$class)
df_forecast_ratio$class <- factor(df_forecast_ratio$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

forecast_ratio <- predict(RatioModel, df_forecast_ratio, type = "link", se.fit = TRUE)
df_forecast_ratio$forecast <- forecast_ratio$fit
df_forecast_ratio$upper <- forecast_ratio$fit + 1.96 * forecast_ratio$se.fit
df_forecast_ratio$lower <- forecast_ratio$fit - 1.96 * forecast_ratio$se.fit
Fig1B <- ggplot(df_forecast_ratio, aes(x = LUMIPULSE_CSF_AB42_AB40, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() + geom_vline(xintercept = 0.0673, linetype = "dashed", colour = "grey67") +
  geom_hline(yintercept = 42.5, linetype = "dashed", colour = "grey67") +
  ylab("CSF pTau") + xlab("CSF AB42/AB40") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster") + scale_x_reverse()  +
  geom_point(csf_surv, mapping=aes(x = LUMIPULSE_CSF_AB42_AB40, y = LUMIPULSE_CSF_pTau, colour = class, fill = class)) +  
  geom_line(csf_surv, mapping=aes(x = LUMIPULSE_CSF_AB42_AB40, y = LUMIPULSE_CSF_pTau, colour = class, fill = class, group = ID), alpha = 0.7) +
  theme(legend.position = "bottom", legend.text = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.title = element_text(size = 20), 
        legend.title = element_blank(), title = element_text(size = 24)) + ggtitle(label = "B.")


Cairo(file="./Results/Figure1.png",
      type="png",
      units="px", 
      width=1800, 
      height=1000, 
      pointsize=12, 
      dpi="auto")
grid.arrange(Fig1A, Fig1B, nrow = 2, ncol = 1)
dev.off()


Cairo(file="./Results/ManuscriptFigures/ThumbNail.png",
      type="png",
      units="px", 
      width=1800, 
      height=1400, 
      pointsize=12, 
      dpi="auto")
Fig1A + ggtitle(label = "")
dev.off()


###########AMYLOID
##CSF Abeta
AbetaModel <- bam(LUMIPULSE_CSF_AB42_AB40 ~ s(Age, bs = "cr") + 
                    s(Age, class, bs = "fs", m = 1, k = 3),
                  data = csf_surv, method = "fREML")
plot_smooth(AbetaModel, view="Age", plot_all="class",
            rug=F, rm.ranef=T)

plot(AbetaModel)
summary.coefs(AbetaModel)
plot_smooth.cont(AbetaModel, view="Age", plot_all.c="class",
                 rug=F, rm.ranef=T)

df_forecast <- get_makeForecast(data.frame(csf_surv), 
                                ClassVec = c("AD Biomarker Negative", "Intermediate AD Biomarkers", "AD Biomarker Positive"), 
                                Model = AbetaModel)
df_forecast$class <- factor(df_forecast$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

Plot3A <- ggplot(df_forecast, aes(x = Age, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() + geom_hline(yintercept = 0.0673, linetype = "dashed", colour = "grey67") +
  xlab("Age at LP") + ylab("CSF AB42 / AB40 Ratio") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster") +
  geom_point(csf_surv, mapping=aes(x = Age, y = LUMIPULSE_CSF_AB42_AB40, colour = class, fill = class), alpha = 0.5) +  
  geom_line(csf_surv, mapping=aes(x = Age, y = LUMIPULSE_CSF_AB42_AB40, colour = class, fill = class, group = ID), alpha = 0.3) +
  ggtitle("A.")

get_AgeAtPositive(df_forecast, "AD Biomarker Negative", 0.0673, "<")  
get_AgeAtPositive(df_forecast, "Intermediate AD Biomarkers", 0.0673, "<")  
get_AgeAtPositive(df_forecast, "AD Biomarker Positive", 0.0673, "<")  
#Age at Abetqa positive = 72.3 (95%CI: 70.1, 75.4)

## PIB PET
pet <- merge(pet[, c("ID", "PET_Date", "pib_fsuvr_rsf_tot_cortmean")], csf_surv[, c("ID", "class", "BIRTH")],
             by = "ID", all = FALSE)
pet <- pet[!duplicated(pet),]
pet$Age <- as.numeric(as.Date(pet$PET_Date, format = "%Y-%d-%m") - pet$BIRTH) / 365.25
pet <- pet[!is.na(pet$pib_fsuvr_rsf_tot_cortmean),]


PIBModel <- bam(pib_fsuvr_rsf_tot_cortmean ~ s(Age, bs = "cr") + 
                  s(Age, class, bs = "fs", m = 1, k = 3),
                data = pet, method = "fREML")
plot_smooth(PIBModel, view="Age", plot_all="class",
            rug=F, rm.ranef=T)

summary.coefs(PIBModel)

df_forecast_pib <- get_makeForecast(pet, ClassVec = c("AD Biomarker Negative", "Intermediate AD Biomarkers", "AD Biomarker Positive"), 
                                    Model = PIBModel)
pet$class <- factor(pet$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))
df_forecast_pib$class <- factor(df_forecast_pib$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

Plot3B <- ggplot(df_forecast_pib, aes(x = Age, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() + geom_hline(yintercept =1.42, linetype = "dashed", colour = "grey67") +
  xlab("Age at Scan") + ylab("PET PIB Summary Value") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster")  +
  geom_point(pet, mapping=aes(x = Age, y = pib_fsuvr_rsf_tot_cortmean, colour = class, fill = class), alpha = 0.5) +  
  geom_line(pet, mapping=aes(x = Age, y = pib_fsuvr_rsf_tot_cortmean, colour = class, fill = class, group = ID), alpha = 0.3) +
  ggtitle("B.")

get_AgeAtPositive(df_forecast_pib, "AD Biomarker Negative", 1.41999, ">")  
get_AgeAtPositive(df_forecast_pib,  "Intermediate AD Biomarkers", 1.41999, ">")  
get_AgeAtPositive(df_forecast_pib, "AD Biomarker Positive", 1.41999, ">")  


#pTau
pTauModel <- bam(LUMIPULSE_CSF_pTau ~ s(Age, bs = "cr") + 
                   s(Age, class, bs = "fs", m = 1, k = 3),
                 data = csf_surv, method = "fREML")
plot_smooth(pTauModel, view="Age", plot_all="class",
            rug=F, rm.ranef=T)

summary.coefs(pTauModel)

df_forecast_pTau <- get_makeForecast(data.frame(csf_surv), ClassVec = c("AD Biomarker Negative", "Intermediate AD Biomarkers", "AD Biomarker Positive"), 
                                     Model = pTauModel)
df_forecast_pTau$class <- factor(df_forecast_pTau$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

Plot3C <- ggplot(df_forecast_pTau, aes(x = Age, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() + geom_hline(yintercept = 42.5, linetype = "dashed", colour = "grey67") +
  xlab("Age at LP") + ylab("CSF pTau") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster") +
  geom_point(csf_surv, mapping=aes(x = Age, y = LUMIPULSE_CSF_pTau, colour = class, fill = class), alpha = 0.5) +  
  geom_line(csf_surv, mapping=aes(x = Age, y = LUMIPULSE_CSF_pTau, colour = class, fill = class, group = ID), alpha = 0.3) +
  ggtitle("C.")

get_AgeAtPositive(df_forecast_pTau, "AD Biomarker Negative", 42.4999, ">")  
get_AgeAtPositive(df_forecast_pTau, "Intermediate AD Biomarkers", 42.4999, ">")  
get_AgeAtPositive(df_forecast_pTau, "AD Biomarker Positive", 42.4999, ">")  

#PET Tau
tau <- merge(tau[, c("ID", "PET_Date", "Tauopathy", "Age_at_Scan")], csf_surv[, c("ID", "class")])
tau <- tau[!duplicated(tau),]
tau <- tau[!is.na(tau$Tauopathy),]
colnames(tau)[4] <- "Age"
tauModel <- bam(Tauopathy ~ s(Age, bs = "cr") + 
                  s(Age, class, bs = "fs", m = 1, k = 3),
                data = tau, method = "fREML")
plot_smooth(tauModel, view="Age", plot_all="class",
            rug=F, rm.ranef=T)

summary.coefs(tauModel)

df_forecast_tau <- get_makeForecast(data.frame(tau), ClassVec = c("AD Biomarker Negative", "Intermediate AD Biomarkers", "AD Biomarker Positive"), 
                                    Model = tauModel)
tau$class <- factor(tau$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))
df_forecast_tau$class <- factor(df_forecast_tau$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

Plot3D <- ggplot(df_forecast_tau, aes(x = Age, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() + geom_hline(yintercept = 1.22, linetype = "dashed", colour = "grey67") +
  xlab("Age at Scan") + ylab("PET AV1451 Summary Value") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster") +
  geom_point(tau, mapping=aes(x = Age, y = Tauopathy, colour = class, fill = class), alpha = 0.5) +  
  geom_line(tau, mapping=aes(x = Age, y = Tauopathy, colour = class, fill = class, group = ID), alpha = 0.3) +
  ggtitle("D.")

get_AgeAtPositive(df_forecast_tau, "AD Biomarker Negative", 1.22, ">")  
get_AgeAtPositive(df_forecast_tau, "Intermediate AD Biomarkers", 1.22, ">")  
get_AgeAtPositive(df_forecast_tau, "AD Biomarker Positive", 1.22, ">")  


mri <- merge(mri[, c("ID", "CortSig_Thickness", "Age_at_Scan", "session_date")], csf_surv[, c("ID", "class")], by = "ID", all = FALSE)
mri <- mri[!duplicated(mri),]
colnames(mri)[3] <- "Age"
mriModel <- bam(CortSig_Thickness ~ s(Age, bs = "cr") + 
                  s(Age, class, bs = "fs", m = 1, k = 3),
                data = mri, method = "fREML")

df_forecast_mri <- get_makeForecast(data.frame(mri), ClassVec = c("AD Biomarker Negative", "Intermediate AD Biomarkers", "AD Biomarker Positive"), 
                                    Model = mriModel)
mri$class <- factor(mri$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))
df_forecast_mri$class <- factor(df_forecast_mri$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

Plot3E <- ggplot(df_forecast_mri, aes(x = Age, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() +
  xlab("Age at Scan") + ylab("Cortical Thickness") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster")  +
  geom_point(mri, mapping=aes(x = Age, y = CortSig_Thickness, colour = class, fill = class), alpha = 0.5) +  
  geom_line(mri, mapping=aes(x = Age, y = CortSig_Thickness, colour = class, fill = class, group = ID), alpha = 0.3) +
  ggtitle("E.")



wmh <- merge(wmh[, c("ID", "WMH_volume", "Age_at_Scan", "session_date")], csf_surv[, c("ID", "class")], by = "ID", all = FALSE)
wmh <- wmh[!duplicated(wmh),]
colnames(wmh)[3] <- "Age"
wmhModel <- bam(WMH_volume ~ s(Age, bs = "cr") + 
                  s(Age, class, bs = "fs", m = 1, k = 3),
                data = wmh, method = "fREML")

df_forecast_wmh <- get_makeForecast(data.frame(wmh), ClassVec = c("AD Biomarker Negative", "Intermediate AD Biomarkers", "AD Biomarker Positive"), 
                                    Model = wmhModel)
wmh$class <- factor(wmh$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))
df_forecast_wmh$class <- factor(df_forecast_wmh$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

Plot3F <- ggplot(df_forecast_wmh, aes(x = Age, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() +
  xlab("Age at Scan") + ylab("WMH Volume") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster")  +
  geom_point(wmh, mapping=aes(x = Age, y = WMH_volume, colour = class, fill = class), alpha = 0.5) +  
  geom_line(wmh, mapping=aes(x = Age, y = WMH_volume, colour = class, fill = class, group = ID), alpha = 0.3) +
  ggtitle("F.")


nfl <- merge(nfl, csf_surv[, c("ID", "BIRTH", "class")], by = "ID", all = FALSE)
nfl <- nfl[!duplicated(nfl),]
nfl$Age <- as.numeric(as.Date(nfl$LP_DATE, format = "%Y-%m-%d") - nfl$BIRTH)/365.25
nflModel <- bam(LN_CSF_NfL ~ s(Age, bs = "cr") +
                  s(Age, class, bs = "fs", m = 1, k = 3),
                data = nfl, method = "fREML")
df_forecast_nfl <- get_makeForecast(data.frame(nfl), ClassVec = c("AD Biomarker Negative", "Intermediate AD Biomarkers", "AD Biomarker Positive"), 
                                    Model = nflModel)
nfl$class <- factor(nfl$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))
df_forecast_nfl$class <- factor(df_forecast_nfl$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))

Plot3G <-
  ggplot(df_forecast_nfl, aes(x = Age, y = forecast, colour = class, fill = class)) + geom_smooth()  +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4)  + scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme_bw() +
  xlab("Age at Scan") + ylab("LN(CSF NfL)") + labs(colour = "Assigned Cluster", fill = "Assigned Cluster") +
  geom_point(nfl, mapping=aes(x = Age, y = LN_CSF_NfL, colour = class, fill = class), alpha = 0.5) +  
  geom_line(nfl, mapping=aes(x = Age, y = LN_CSF_NfL, colour = class, fill = class, group = ID), alpha = 0.3) +
  ggtitle("G.")

nfl_compare <- data.frame("Age" = df_forecast_nfl[df_forecast_nfl$class == "AD Biomarker Negative", "Age"],
                          "Class1_lower" = df_forecast_nfl[df_forecast_nfl$class == "AD Biomarker Negative", "lower"],
                          "Class2_lower" = df_forecast_nfl[df_forecast_nfl$class == "Intermediate AD Biomarkers", "lower"],
                          "Class3_upper" = df_forecast_nfl[df_forecast_nfl$class == "AD Biomarker Positive", "upper"])

################################################################################
##########PRESENTED FIGURES/TABLES########################################
################################################################################


#############TABLE1########################################
#Stuff for the tableone table
csf_comparison <- get_early_late(csf, "CSF_LP_DATE")
mmse_comparison <- get_early_late(mmse, "testdate")
mri_comparison <- get_early_late(mri, "session_date")
np_comparison <- get_early_late(np, "psy_date")
pet_comparison <- get_early_late(pet, "PET_Date")
tau_comparison <- get_early_late(tau, "PET_Date")
wmh_comparison <- get_early_late(wmh, "session_date")

df_tableone <- merge(csf_comparison, mmse_comparison[, c("ID", "cdr_baseline", "MMSE_baseline", 
                                                         "cdr_recent", "MMSE_recent")], by = "ID", all = TRUE)
df_tableone <- merge(df_tableone, mri_comparison[, c("ID", "CortSig_Thickness_baseline", "CortSig_Thickness_recent")], by = "ID", all = TRUE)
df_tableone <- merge(df_tableone, np_comparison[, c("ID", "lmdelay_baseline", "digsym_baseline", "lmdelay_recent", "digsym_recent")], by = "ID", all = TRUE)
df_tableone <- merge(df_tableone, pet_comparison[, c("ID", "pib_fsuvr_rsf_tot_cortmean_baseline", "pib_fsuvr_rsf_tot_cortmean_recent")], by = "ID", all = TRUE)
df_tableone <- merge(df_tableone, tau_comparison[, c("ID", "Tauopathy_baseline", "Tauopathy_recent")], by = "ID", all = TRUE)
df_tableone <- merge(df_tableone, wmh_comparison[, c("ID", "WMH_volume_baseline", "WMH_volume_recent")], by = "ID", all = TRUE)

tmp <- df_tableone[,c("ID", "CSF_LP_DATE_baseline")]
mmse_tmp <- mmse[, c("ID", "testdate", "cdr", "MMSE")]

setDT(tmp)
setDT(mmse_tmp)

setkey(mmse_tmp, ID, testdate)[, CSF_LP_DATE_baseline := testdate]
nearest_np <- mmse_tmp[tmp, roll = 'nearest']
colnames(nearest_np)[ 3:5] <- c("cdr_baseline", "MMSE_baseline", "sumbox_baseline")

df_tableone <- df_tableone[, !(names(df_tableone) %in% c("cdr_baseline", "MMSE_baseline",
                                                         "sumbox_baseline"))]

df_tableone <- merge(df_tableone, nearest_np, by = c("ID"))
df_tableone <- merge(df_tableone, gmm3_2$pprob[, c("ID", "class")], by = "ID") #merging wtih labels from clustering for stratification

df_tableone$enrollment_duration <- df_tableone$Age_at_LP_recent - df_tableone$Age_at_LP_baseline

#number of LPs
get_counts_LP <- data.table(csf)[, .N, by = ID]
colnames(get_counts_LP)[2] <- "NumberLPs"
get_counts_CDR <- data.table(mmse)[, .N, by = ID]
colnames(get_counts_CDR)[2] <-"NumberCDRs"

get_counts_Proteome <- read.csv("./Results/NumProcessedProteomes.csv")
colnames(get_counts_Proteome) <- c("ID", "NumberProteomes")


df_tableone <- merge(df_tableone, get_counts_LP, by = "ID", all.x = TRUE, all.y = FALSE)
df_tableone <- merge(df_tableone, get_counts_CDR, by = "ID", all.x = TRUE, all.y = FALSE)
df_tableone <- merge(df_tableone, get_counts_Proteome, by = "ID", all.x = TRUE, all.y = FALSE)

myVars = c( "Age_at_LP_baseline", "Age_at_LP_recent",
            "enrollment_duration", "GENDER_baseline", 
            "apoe_baseline", "race2_baseline", "EDUC_baseline",
            "cdr_baseline", "cdr_recent", 
            "MMSE_baseline", "MMSE_recent", 
            "NumberLPs", "NumberCDRs", "NumberProteomes")
catVars = c( "GENDER_baseline", 
             "apoe_baseline", "race2_baseline",
             "cdr_baseline", "cdr_recent")
CreateTableOne(vars = myVars, strata = "class", data = df_tableone, factorVars = catVars)
##################################################################################

#############FIGURE1########################################

##################################################################################

#############FIGURE2########################################
km_obj_A[[2]]$class <- factor(km_obj_A[[2]]$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))
km_obj_T[[2]]$class <- factor(km_obj_T[[2]]$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))
km_obj_N[[2]]$class <- factor(km_obj_N[[2]]$class, levels=c("AD Biomarker Positive", "Intermediate AD Biomarkers", "AD Biomarker Negative"))


p1 <-   ggplot(km_obj_A[[2]], aes(x = Age, y = Probability, group = class, colour = class, fill = class)) +  
  geom_line() + theme_bw() + xlab("Age") + ylab("Percent Amyloid Negative") +
  xlim(c(47, 87))+
  scale_colour_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, colour = NA) + ggtitle("A.") +
  theme(legend.text = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.title = element_blank())

p2 <- ggplot(km_obj_T[[2]], aes(x = Age, y = Probability, group = class, colour = class, fill = class)) +  
  geom_line() + theme_bw() + xlab("Age") + ylab("Percent Tau Negative") +
  xlim(c(47, 87))+
  scale_colour_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, colour = NA) + ggtitle("B.")+
  theme(legend.text = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.title = element_blank())

p3 <- ggplot(km_obj_N[[2]], aes(x = Age, y = Probability, group = class, colour = class, fill = class)) +  
  geom_line() + theme_bw() + xlab("Age") + ylab("Percent CDR = 0") +
  xlim(c(47, 87))+
  scale_colour_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, colour = NA) + ggtitle("C.")+
  theme(legend.text = element_text(size = 16), 
        axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.title = element_blank())

Cairo(file="./Results/Figure2.png",
      type="png",
      units="px", 
      width=1800, 
      height=1200, 
      pointsize=12, 
      dpi="auto")
lemon::grid_arrange_shared_legend(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()
##################################################################################

######FIGURE 3###########################################
lay <- rbind(c(1, 1, 1, 2, 2, 2), c(3, 3, 3, 4, 4, 4), c(5, 5, 6, 6, 7, 7))
Cairo(file="./Results/Figure3.png",
      type="png",
      units="px", 
      width=1800, 
      height=1000, 
      pointsize=12, 
      dpi="auto")
grid.arrange(Plot3A + theme(legend.position="none"), Plot3B + theme(legend.position="none"), 
             Plot3C + theme(legend.position="none"), Plot3D + theme(legend.position="none"),
             Plot3E + theme(legend.position="none"), Plot3F + theme(legend.position = "bottom"), 
             Plot3G + theme(legend.position="none"), layout_matrix = lay)
dev.off()
##################################################################################
write.csv(csf_surv, "./Data/csf_processed_for_analysis.csv", row.names = FALSE) #saving for use in IndividualProteinAnalysis
