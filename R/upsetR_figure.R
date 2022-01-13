library(UpSetR)
library(data.table)

csf <- read.csv("./Data/csf_processed_for_analysis.csv")
mri <- read_xlsx("././Data/mod_mri.xlsx")
tau <- read_xlsx("././Data/mod_tau.xlsx")
wmh <- read_xlsx("././Data/mod_wmh.xlsx")
pet <- read_xlsx("././Data/mod_pet.xlsx")
mmse <- read_xlsx("././Data/mod_b4_cdr.xlsx")
demogs <- read_xlsx("././Data/mod_demographics.xlsx")

visits <- data.table(csf)[, .(count = length(upvoted_cluster)), by = ID]
print(paste0("Mean number of CSF visits = ", mean(visits$count)))
print(paste0("SD number of CSF visits = ", sd(visits$count)))
csf <- csf[!duplicated(csf$ID), c("ID", "upvoted_cluster")]

mmse <- mmse[mmse$ID %in% csf$ID,]
visits <- data.table(mmse)[, .(count = length(cdr)), by = ID]
print(paste0("Mean number of visits = ", mean(visits$count)))
print(paste0("SD number of visits = ", sd(visits$count)))


get_modality <- function(DF, ID_name = "ID", QC_name = "FS_QC_Status", Measure_name = "MR_TOTV_INTRACRANIAL"){
  DF <- data.frame(DF)
  DF <- DF[DF[, QC_name] != "Failed", ]
  DF <- DF[DF[, ID_name] %in% csf$ID, c(ID_name, Measure_name)]
  DF <- DF[!is.na(DF[,Measure_name]),]
  
  return(DF)
}

mri <- get_modality(mri, "ID", "FS_QC_Status", "MR_TOTV_INTRACRANIAL")
visits <- data.table(mri)[, .(count = length(MR_TOTV_INTRACRANIAL)), by = ID]
print(paste0("Mean number of visits = ", mean(visits$count)))
print(paste0("SD number of visits = ", sd(visits$count)))
mri <- mri[!duplicated(mri$ID),]

pet <- get_modality(pet, "ID", "PUP_QC_Status", "pib_fsuvr_tot_cortmean")
visits <- data.table(pet)[, .(count = length(pib_fsuvr_tot_cortmean)), by = ID]
print(paste0("Mean number of visits = ", mean(visits$count)))
print(paste0("SD number of visits = ", sd(visits$count)))
pet <- pet[!duplicated(pet$ID),]

tau <- get_modality(tau, "ID", "PUP_QC_Status", "Tauopathy")
visits <- data.table(tau)[, .(count = length(Tauopathy)), by = ID]
print(paste0("Mean number of visits = ", mean(visits$count)))
print(paste0("SD number of visits = ", sd(visits$count)))
tau <- tau[!duplicated(tau$ID),]


wmh <- wmh[wmh$ID %in% csf$ID, c("ID", "WMH_volume")]
visits <- data.table(wmh)[, .(count = length(WMH_volume)), by = ID]
print(paste0("Mean number of visits = ", mean(visits$count)))
print(paste0("SD number of visits = ", sd(visits$count)))
wmh <- wmh[!duplicated(wmh$ID),]

upset_df <- list(`Longitudinal CSF` = csf$ID,
                 `Cross-Sectional Proteomics` = csf$ID,
                `Longitudinal Cognitive Screening` = csf$ID,
                 `PET - PIB` = pet$ID,
                 `PET - AV1451` = tau$ID,
                 `Structural MRI` = mri$ID,
                 WMH = wmh$ID)
upset(fromList(upset_df), order.by = "freq",
      sets = c("Cross-Sectional Proteomics", "Longitudinal Cognitive Screening", "WMH", 
               "Structural MRI","PET - AV1451", "PET - PIB", "Longitudinal CSF"),
      keep.order = TRUE)
