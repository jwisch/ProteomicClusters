get_early_late <- function(DF, DateName){
  early <- setDT(DF)[order(get(DateName)), head(.SD, 1L), by = ID]
  late <- setDT(DF)[order(get(DateName)), tail(.SD, 1L), by = ID]
  
  early <- data.frame(early)
  colnames(early) <- paste0(names(early), "_baseline")
  colnames(early)[1] <- "ID"
  
  late <- data.frame(late)
  colnames(late) <- paste0(names(late), "_recent")
  colnames(late)[1] <- "ID"
  
  df_return <- merge(early, late, by = "ID")
  return(df_return)
}

#Function for date matching
MatchbyNearestDate<-function(df1, df2, ID, Date1, Date2){
  z <- lapply(intersect(df1[,ID],df2[,ID]),function(id) {
    df1 <- df1[df1[,ID] == id,]
    df2 <- df2[df2[,ID] == id,]
    
    df1[,"indices"] <- sapply(df1[,Date1],function(d) which.min(abs(df2[,Date2] - d)))
    df2[,"indices"] <- 1:nrow(df2)
    
    merge(df1,df2,by=c(ID,'indices'))
  })
  
  df_matched <- do.call(rbind,z)
  df_matched$indices <- NULL
  return(df_matched)
  
}

#This calculates PACC change in an equivalent manner to the Donohue et al paper
#I don't think it's actually appropriate for this project
get_PACC_change <- function(np, np_baseline){
  np_baseline <- np_baseline[!duplicated(np_baseline$ID),]
  np_change_from_base <- np %>% 
    group_by(ID) %>% 
    arrange(ID, psy_date) %>% 
    mutate(srttotal_changefrombaseline = (srttotal - first(na.omit(srttotal))) / sd(np[,"srttotal"], na.rm = TRUE),
           lmdelay_changefrombaseline = (lmdelay - first(na.omit(lmdelay))) / sd(np[,"lmdelay"], na.rm = TRUE),
           digsym_changefrombaseline = (digsym - first(na.omit(digsym))) / sd(np[,"digsym"], na.rm = TRUE),
           MEMUNITS_changefrombaseline = (MEMUNITS - first(na.omit(MEMUNITS))) / sd(np[,"MEMUNITS"], na.rm = TRUE),
           MMSE_changefrombaseline = (MMSE - first(na.omit(MMSE))) / sd(np[,"MMSE"], na.rm = TRUE))
  np_change_from_base$PACC <- as.vector(ifelse(is.na(np_change_from_base[,"MEMUNITS"]),
                                     rowMeans(np_change_from_base[, c("srttotal_changefrombaseline",
                                                                      "lmdelay_changefrombaseline",
                                                                      "digsym_changefrombaseline",
                                                                      "MMSE_changefrombaseline")], na.rm = TRUE) * 4,
                                     rowMeans(np_change_from_base[, c("srttotal_changefrombaseline",
                                                                      "MEMUNITS_changefrombaseline",
                                                                      "digsym_changefrombaseline",
                                                                      "MMSE_changefrombaseline")], na.rm = TRUE) * 4))
  return(np_change_from_base)
  
}


get_barplot <- function(DF){
  DF[, "Group"] <- as.factor(DF[, "Group"])
  
  level_names <- paste0( "Protein Cluster ", seq(from = 1, to = length(unique(DF[, "Group"]))))
  levels(DF[, "Group"]) <- c(level_names)
  
  DF <- data.frame(data.table(DF)[, .(OrderingCrit = sum(Crit),
                                      Crit = Crit,
                                      Group = Group), by = Gene])
  
  p <- ggplot(DF, aes(x = reorder(Gene, OrderingCrit), y = Crit)) + 
    geom_col(aes(fill = Group)) + xlab("Gene") + ylab("Criterion") +
    coord_flip() + theme_bw() + theme(legend.position = "bottom") + 
    scale_fill_viridis(name = " ", discrete = TRUE) 
  return(p)}


get_survival_obj <- function(csf_surv, status_text){
  #Amyloid positive before study starts:
  csf_time_to_amyloid_right_censor <- setDT(csf_surv)[order(-CSF_LP_DATE), head(.SD, 1L), by = ID]
  csf_time_to_amyloid_right_censor <- data.frame(csf_time_to_amyloid_right_censor)[data.frame(csf_time_to_amyloid_right_censor)[, status_text] == 1 & 
                                                                                     !(csf_time_to_amyloid_right_censor$ID %in% 
                                                                                         csf_time_to_amyloid_left_censor$ID),]
  #Getting average of all visits...if all 0 - stay negative whole time, if all 1 - stay positive whole time
  tmp <- setDT(csf_surv)[, mean(get(status_text)), by = ID]
  
  
  #These guys start negative and finish negative
  left_censors <-  tmp[tmp$V1 == 0, "ID"]
  left_censors$ID <- as.numeric(as.character(left_censors$ID))
  csf_surv <- data.frame(csf_surv)
  csf_surv$ID <- as.numeric(as.character(csf_surv$ID))
  csf_time_to_amyloid_left_censor <- csf_surv[(csf_surv$ID %in% left_censors$ID),]
  csf_time_to_amyloid_left_censor <- setDT(csf_time_to_amyloid_left_censor)[order(-CSF_LP_DATE), head(.SD, 1L), by = ID] 
  
  
  #Identify converters and their first converted visit...
  csf_converters <- csf_surv[!(csf_surv$ID %in% csf_time_to_amyloid_left_censor$ID) &
                               !(csf_surv$ID %in% csf_time_to_amyloid_right_censor$ID),]
  #Getting all Apos visits for each converter....
  csf_converters <- csf_converters[csf_converters[, status_text] == 1,]
  #Then grabbing earliest date for each of these guys
  csf_converters <- setDT(csf_converters)[order(CSF_LP_DATE), head(.SD, 1L), by = ID]
  
  
  #TODO: Combine the three dataframes in a way that matches with the stack overflow post, then run survival analysis
  csf_time_to_amyloid_left_censor$Converts <- 0 #they never converted to Apos
  csf_time_to_amyloid_left_censor$Followup <- csf_time_to_amyloid_left_censor$Time_from_Baseline #how long they were followed
  csf_time_to_amyloid_left_censor$Followup2 <- NA #they never experienced event
  
  csf_time_to_amyloid_right_censor$Converts <- 1 #They started out Apos
  csf_time_to_amyloid_right_censor$Followup <- NA #experienced event before study
  csf_time_to_amyloid_right_censor$Followup2 <- csf_time_to_amyloid_right_censor$Time_from_Baseline #how long they were followed
  
  csf_converters$Converts <- 1 #They convert during study
  csf_converters$Followup <- csf_converters$Time_from_Baseline
  csf_converters$Followup2 <- csf_converters$Time_from_Baseline
  
  surv_frame <- rbind(csf_time_to_amyloid_left_censor[, c("Converts", "Followup", "Followup2", "upvoted_cluster")],
                      csf_time_to_amyloid_right_censor[, c("Converts", "Followup", "Followup2", "upvoted_cluster")],
                      csf_converters[, c("Converts", "Followup", "Followup2", "upvoted_cluster")])
  
  
  Surv.Obj <- Surv(surv_frame$Followup, surv_frame$Followup2, type = 'interval2')
  
  km <- survfit(Surv.Obj ~ surv_frame$upvoted_cluster, conf.type = "none")
  
  library(icenReg)
  #https://www.jstatsoft.org/article/view/v081i12
  surv_frame$Followup[is.na(surv_frame$Followup)] <- -Inf
  surv_frame$Followup2[is.na(surv_frame$Followup2)] <- Inf
  fit_A <- ic_sp(cbind(Followup, Followup2) ~ upvoted_cluster, data = surv_frame,
                 model = 'ph', bs_samples = 500)
  
  return(list(km, fit_A))}

get_plottable_survobj <- function(km_list){
  res <- summary(km_list)
  cols <- lapply(c(2:6, 8:11) , function(x) res[x])
  tbl <- do.call(data.frame, cols)
  tbl$upper <- tbl$surv + 1.96*tbl$std.err
  tbl$lower <- tbl$surv - 1.96*tbl$std.err
  return(tbl)
}


#Survival Analysis Functions


get_modality_survival_df <- function(DF, Date, Status){
  csf_surv <- data.table(DF)
  max <- csf_surv[csf_surv[, .I[which.max(get(Status))], by=ID]$V1]#earliest date of the maximum amyloid status attained
  max <- data.frame(max)[, c("ID", "Age", Status, "class")]
  min <- csf_surv[csf_surv[, .I[which.min(get(Status))], by=ID]$V1]#earliest date of the minimum amyloid status attained
  min <- data.frame(min)[, c("ID", "Age", Status, "class")]
  
  csf_surv_X <- csf_surv[order(get(Date)), head(.SD, 1L), by = ID]
  csf_surv_final <- csf_surv[order(-get(Date)), head(Age, 1L), by = ID]
  colnames(csf_surv_final) <- c("ID", "MaxAge")
  csf_surv_X <- merge(csf_surv_X, csf_surv_final, by = "ID")
  csf_surv_X$TotalEnrollment <- csf_surv_X$MaxAge - csf_surv_X$Age
  
  
  #borrowing this example https://stackoverflow.com/questions/41968606/left-censoring-for-survival-data-in-r
  csf_surv_X$ChangeTime <- max$Age - min$Age
  
  csf_surv_X<-csf_surv_X%>%mutate(Status = case_when(
    ChangeTime == 0 & get(Status)==1 ~ "LeftCensor",
    ChangeTime == 0 & get(Status)==0 ~ "RightCensor",
    TRUE~"Converts"
  ))
  
  csf_surv_X$FollowupTime <- ifelse(csf_surv_X$Status == "LeftCensor", NA, #If left censored, use NA
                                    ifelse(csf_surv_X$Status == "RightCensor", csf_surv_X$MaxAge, #if right censored, use oldest age of enrollment
                                           csf_surv_X$Age + csf_surv_X$ChangeTime)) #if they actually change, use age at conversion
  
  csf_surv_X$FollowupTime2 <- ifelse(csf_surv_X$Status == "LeftCensor", csf_surv_X$MaxAge, #If left censored, use oldest age of enrollment
                                     ifelse(csf_surv_X$Status == "RightCensor", NA, #if right censored, use NA
                                            csf_surv_X$Age + csf_surv_X$ChangeTime)) #if they actually change, use age at conversion
  
  return(csf_surv_X)
}


get_run_survival_analysis <- function(csf_surv_X){
  
  Followup1 <- csf_surv_X$FollowupTime
  Followup2 <- csf_surv_X$FollowupTime2
  
  Surv.Obj <- Surv(Followup1, Followup2, type = 'interval2') #- indicates left censor, + indicates right censor
  km <- survfit(Surv.Obj ~ csf_surv_A$class, conf.type = "log-log")
  
  plot_df <- data.frame("Age" = km$time, 
                        "Probability" = km$surv,
                        "min" = km$lower,
                        "max" = km$upper,
                        "class" = as.factor(c(rep("AD Biomarker Positive", as.numeric(km$strata[1])),
                                              rep("Intermediate AD Biomarkers", as.numeric(km$strata[2])),
                                              rep("AD Biomarker Negative", as.numeric(km$strata[3])))))
  #hardcoding confidence intervals at the extremes
  plot_df$min[ plot_df$Probability == 0] <- 0
  plot_df$max[ plot_df$Probability == 0] <- 0
  plot_df$min[ plot_df$Probability == 1] <- 1
  plot_df$max[ plot_df$Probability == 1] <- 1
  
  return(list(km, plot_df))
}


get_proportion_monotonic <- function(DF, IDCOL, YCOL, DIRECTION){
  tmp <- data.frame(monotonic(DF, id.col = IDCOL, y.col = YCOL, direction = DIRECTION))
  return(sum(tmp$Montonic) / nrow(tmp))
  
}

get_AgeAtPositive <- function(DF, CLASS, CUTOFF, DIRECTION){
  ifelse(DIRECTION == ">",
         print(paste0("Mean Age = ", 
                      head(DF[DF[, "class"]== CLASS & DF[,"forecast"] > CUTOFF,"Age"])[1],
                      ", 95% CI = (",
                      head(DF[DF[, "class"]== CLASS & DF[,"upper"] > CUTOFF,"Age"])[1],
                      ", ",
                      head(DF[DF[, "class"]== CLASS & DF[,"lower"] > CUTOFF,"Age"])[1],
                      ")")),
         print(paste0("Mean Age = ", 
                      head(DF[DF[, "class"]== CLASS & DF[,"forecast"] < CUTOFF,"Age"])[1],
                      ", 95% CI = (",
                      head(DF[DF[, "class"]== CLASS & DF[,"lower"] < CUTOFF,"Age"])[1],
                      ", ",
                      head(DF[DF[, "class"]== CLASS & DF[,"upper"] < CUTOFF,"Age"])[1],
                      ")"))         
         
  )   }

get_makeForecast <- function(DF, AgeCol = "Age", ClassVec = c(1, 2, 3), Model){
  df_forecast <- expand.grid("Age" = seq(min(DF[,AgeCol]), max(DF[,AgeCol]), length.out = 50),
                             "class" = ClassVec)
  forecast <- predict(Model, df_forecast, type = "link", se.fit = TRUE)
  df_forecast$forecast <- forecast$fit
  df_forecast$upper <- forecast$fit + 1.96 * forecast$se.fit
  df_forecast$lower <- forecast$fit - 1.96 * forecast$se.fit
  df_forecast$class <- as.factor(df_forecast$class)
  return(df_forecast)
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#median imputation for missingness

f=function(x){
  x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
  x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
  x #display the column
}

get_MultipleClustersAUC <- function(DF, train_IDs, NumberClust, LEVELS, IncCovariates = FALSE){
  
  DF <- data.frame(DF)
  csf_df_available <- data.frame(DF[DF[,"MAPID"] %in% train_IDs,])
  csf_df_holdout <- data.frame(DF[!(DF[,"MAPID"]%in% train_IDs),])

  y_holdout <- csf_df_holdout$upvoted_cluster
 # csf_analysis_holdout <- csf_df_holdout[, 3:715]
  
  
  #need to run below for each comparison
  result_list<- list()
  storeMod_list <- list()
  
  csf_analysis <- csf_df_available
  
  levels(csf_analysis$upvoted_cluster) <- LEVELS
  csf_analysis <- csf_analysis[!is.na(csf_analysis$upvoted_cluster),]
  
  
  
  #split to train/test
  df_split <- initial_split(csf_analysis, strata = upvoted_cluster, prop = 0.75)
  csf_analysis_train <- data.frame(training(df_split))
  
 if(IncCovariates == FALSE){
    csf_analysis_test <- data.frame(testing(df_split))[,3:715]
    csf_analysis_holdout <- data.frame(csf_df_holdout)[,c(3:715)]
  }
  if(IncCovariates == TRUE){
    csf_analysis_test <- data.frame(testing(df_split))[,c(3:715)]
    covariates_test <- data.frame(testing(df_split))[,718:719]
    csf_analysis_holdout <- data.frame(csf_df_holdout)[,c(3:715)]
    covariates_holdout <- data.frame(csf_df_holdout[, 718:719])
  }

 y_test <- data.frame(testing(df_split)$upvoted_cluster)


  y <- csf_analysis_train$upvoted_cluster
  csf_analysis_train <- standardize.genes(csf_analysis_train[,3:715]) 
  covariates <- standardize.genes(data.frame(training(df_split))[,c("GENDER", "Age")])

  #applying standardization to test & holdout sets
  for (j in 1:(dim(csf_analysis_test)[2])) {
    csf_analysis_test[, j] <- (csf_analysis_test[, j] - csf_analysis_train$means[j])/csf_analysis_train$sdevs[j]
    csf_analysis_holdout[, j] <- (csf_analysis_holdout[,j] - csf_analysis_train$means[j]) / csf_analysis_train$sdevs[j]
  }

for(j in 1:dim(covariates_test)[2]){
  covariates_test[, j] <- (covariates_test[,j] - covariates$means[j]) / covariates$sdevs[j]
  covariates_holdout[, j] <- (covariates_holdout[,j] - covariates$means[j]) / covariates$sdevs[j]
  
}

  #AD+ vs Emerging
  for(i in 1:NumberClust){
    
    NumClust <- i
    df_split$data$GENDER <- as.numeric(as.character(df_split$data$GENDER))

   if(IncCovariates == FALSE){
      mod_noCovariates <- pelora(x = as.matrix(csf_analysis_train$x), y = as.numeric(as.character(y)), 
                                 noc = NumClust,
                                 standardize = FALSE)
      forecast <- predict(mod_noCovariates, csf_analysis_test, type = "class")
      probs <- predict(mod_noCovariates, csf_analysis_test,  type = "prob")
     }
    if(IncCovariates == TRUE){
      mod_noCovariates <- pelora(x = as.matrix(csf_analysis_train$x), y = as.numeric(as.character(y)),
                                 u = as.matrix(covariates$x),
                                 noc = NumClust,
                                 standardize = FALSE)
      forecast <- predict(mod_noCovariates, csf_analysis_test, covariates_test, type = "class")
      probs <- predict(mod_noCovariates, csf_analysis_test, covariates_test, type = "prob")
    }


    
    forecast <- data.frame("forecast" = unlist(forecast),
                           "probability" = unlist(probs),
                           "actual" = y_test)
    colnames(forecast) <- c("forecast", "probability", "actual")
    table(forecast$forecast, forecast$actual)
    
    rocobj <- roc(forecast$actual, forecast$probability)
    
    result <- data.frame("NumClust" = NumClust,
                         "AUC" = as.numeric(ci(rocobj))[2],
                         "AUC_min" = as.numeric(ci(rocobj))[1],
                         "AUC_max" = as.numeric(ci(rocobj))[3])
    result_list[[i]] <- result 
    storeMod_list[[i]] <- mod_noCovariates 
  }
  
  return(list(rbindlist(result_list), storeMod_list, csf_analysis_holdout, y_holdout, covariates_holdout))}


get_importantGenes <- function(GENES, CRIT, GENE_NAMES){
  df <- data.frame("gene_no" = unlist(GENES),
                   "score" = unlist(CRIT))
  df$gene <- unlist(GENE_NAMES)[df$gene_no]
  df_sum <- data.table(df)[, sum(score), by = gene]
  df_sum <- data.frame(df_sum)
  colnames(df_sum)[2] <- "score"
  return(df)
  
}

get_LassoAUC <- function(LEVELS, DF, holdout_IDs){
  csf_binom <- DF[!(DF[,"MAPID"] %in% holdout_IDs[,"MAPID"]),]
  levels(csf_binom$upvoted_cluster) <- LEVELS
  cvfit <- cv.glmnet(as.matrix(csf_binom[!is.na(csf_binom$upvoted_cluster),c("Age", "GENDER")]), 
                     as.numeric(as.character(csf_binom[!is.na(csf_binom$upvoted_cluster),]$upvoted_cluster)))
  coef(cvfit, s = "lambda.min")
  
  csf_binom_holdout <- csf_df[csf_df$MAPID %in% holdout_IDs$MAPID,]
  levels(csf_binom_holdout$upvoted_cluster) <- LEVELS
  
  holdout_results <- data.frame("Actual" = csf_binom_holdout[!is.na(csf_binom_holdout$upvoted_cluster),c("upvoted_cluster")],
                                "Forecast" = 
                                  predict(cvfit, as.matrix(csf_binom_holdout[!is.na(csf_binom_holdout$upvoted_cluster),c("Age", "GENDER")]),
                                          s = "lambda.min", type = "response"))
  colnames(holdout_results) <- c("Actual", "Forecast")
  rocobj <- roc(holdout_results$Actual, holdout_results$Forecast)
  return(rocobj)
}
