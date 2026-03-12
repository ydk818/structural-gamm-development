
library(gamm4)  #
library(readxl)  
library(ggplot2)  
library(readr)
library(dplyr)
library(mgcv)
library(haven)
library(sjPlot)
library(caret)


### Load structural data
qc_data <- read_excel("/gpfs2/scratch/melmarsr/DK/data/datasetupdate.xlsx", sheet = "QC_DAIRC")
names(qc_data)[names(qc_data) == "age...10"] <- "age"
names(qc_data)[names(qc_data) == "PGUID...2"] <- "PGUID"
label_map <- read_csv("/gpfs2/scratch/melmarsr/DK/data/label_map.csv")
shn_data <- read_csv("/gpfs2/scratch/melmarsr/ABCD7.0/tabulated/abcd_imgincl01.csv")
shn_data $Visit_ID <- gsub(paste0("_year", "$"), "year", shn_data $Visit_ID)
shn_data <- subset(shn_data, select = c("Visit_ID", "SNH"))
names(shn_data)[names(shn_data) == "Visit_ID"] <- "ID_event"
names(shn_data)[names(shn_data) == "SNH"] <- "SHN"
qc_data <- merge(qc_data,shn_data, by='ID_event')

#thres<-15
#shn_id <- shn_data[shn_data$SHN < thres,]$ID_event
qc_id<- read_csv("/gpfs2/scratch/melmarsr/DK/data/abcd_imgincl01.csv")
#qc_data <- qc_data %>% filter(ID_event%in% shn_id)
qc_index=qc_id[qc_id$imgincl_t1w_include==1,]$VisitID

### Filter the flagged subjects and load covariates 
df_cov= subset(qc_data, select = c("VisitID", "sex", "age","PGUID","S_S","ID_event","SHN"))
df_cov= df_cov %>% filter(VisitID %in% qc_index)
#df_cov= df_cov %>% filter(ID_event%in% shn_id)

### Find labels for structural data
structural_mapper=na.omit(label_map[label_map$table_name_nda =="smrip102",])
structural_list=structural_mapper$var_name
structural_dk <- read_csv("/gpfs2/scratch/melmarsr/DK/data/ABCD_sync/abcd_smrip102.csv")

### winsorized
winsorize_3sd <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  lower_bound <- mean_x - 3 * sd_x
  upper_bound <- mean_x + 3 * sd_x
  x[x < lower_bound] <- lower_bound
  x[x > upper_bound] <- upper_bound
  return(x)
}
### filter outliers
filter_outliers <- function(df) {
  # Identify numeric columns
  numeric_cols <- sapply(df, is.numeric)
  
  # Subset numeric data
  num_df <- df[, numeric_cols]
  
  # Calculate mean and SD
  mu <- sapply(num_df, mean, na.rm = TRUE)
  sigma <- sapply(num_df, sd, na.rm = TRUE)
  
  # Logical matrix: TRUE if within 3 SDs
  within_bounds <- mapply(function(x, m, s) {
    abs(x - m) <= 5 * s
  }, num_df, mu, sigma, SIMPLIFY = FALSE)
  
  # Combine all TRUEs row-wise
  keep_rows <- Reduce(`&`, within_bounds)
  
  # Return filtered full data frame (incl. categorical columns)
  return(df[keep_rows, ])
}

structural_dk_bl <- filter(structural_dk , eventname == "baseline_year_1_arm_1")
structural_dk_2y <- filter(structural_dk , eventname == "2_year_follow_up_y_arm_1")
structural_dk_4y <- filter(structural_dk , eventname == "4_year_follow_up_y_arm_1")
structural_dk_6y <- filter(structural_dk , eventname == "6_year_follow_up_y_arm_1")
structural_dk_bl <- filter_outliers(structural_dk_bl[,0:483])
structural_dk_2y <- filter_outliers(structural_dk_2y[,0:483])
structural_dk_4y <- filter_outliers(structural_dk_4y[,0:483])
structural_dk_6y <- filter_outliers(structural_dk_6y[,0:483])

# Apply the winsorization function to each column
#df_winsorized <- as.data.frame(lapply(task, winsorize_3sd))
structural_dk <-rbind(structural_dk_bl, structural_dk_2y, structural_dk_4y, structural_dk_6y)

# Display result
head(structural_dk)
filter_mapper= label_map %>%  filter(grepl("Cortical thickness in mm",var_label))
roi_list=filter_mapper$var_name


### Load ICV variable
icv_data=read_csv("/gpfs2/scratch/melmarsr/DK/data/ABCD_sync/abcd_smrip102.csv")
icv=subset(icv_data, select = c("VisitID", "smri_vol_scs_intracranialv","subjectkey"))

### Load CBCL data 
phenotype_fix="cbcl_scr_syn_internal_r"
cbcl_data=read_csv("/gpfs2/scratch/melmarsr/DK/data/cbcl_nih_sma_diff.csv")
cbcl_adhd=subset(cbcl_data, select = c("subjectkey", phenotype_fix))

### Load Polygenic Risk Scores
#genetic_data= read_sav("/gpfs1/home/d/y/dyuan/matt_gamm/mds_and_MDD_prs.sav")
#genetic=subset(genetic_data, select = c("subjectid", "S8","C1","C2","C3","C4"))
#names(genetic)[names(genetic) == 'subjectid'] <- 'subjectkey'

### Merge Fix Factors into one Dataframe
fix_factor=merge(icv,cbcl_adhd, by='subjectkey')


# Create a table of counts for each unique item
item_counts <- table(df_cov$PGUID)

# Filter items that appear more than once
visit=1
count=visit-1
filtered_list <- df_cov$PGUID[df_cov$PGUID%in% names(item_counts[item_counts > count])]
df_cov <- df_cov %>% filter(PGUID %in% filtered_list )
str(df_cov)
#unique(df_cov[,4])

roi= roi_list[152]
label_map[label_map$var_name==roi,]

#residualized icv from structural data
#structural_icv=merge(structural_destr,icv, by='VisitID')
#for (roi in roi_list) {
#  model <- lm(structural_icv[[roi]] ~ structural_icv$smri_vol_scs_intracranialv)
#  structural_icv[[paste0(roi, "_resid")]] <- resid(model)
#}


# desikan abcd_smrip102 abcd_smrip102.csv
# subcortical structural_list[429:444][445:459][460:464] data: abcd_smrip102
# surface area structural_list[303:453] data: abcd_mrisdp102
# cortical thickness structural_list[303:453] data: abcd_mrisdp102
# cortical thickness roi_list[0:68] Cortical thickness in mm of APARC ROI 
# cortical area roi_list[152:219]] Cortical area in mm
# cortical volume roi_list[152:219]] Cortical volume in mm

#SA roi_list[152:219]
#SV roi list[1:16][17:30][31:32]
### Residualized ROI from ICV
#structural_icv=merge(structural_dk,icv, by='VisitID')


#roi_names <- setdiff(names(df), "ICV")   # all columns except ICV

#for (roi in roi) {
#  model <- lm(structural_icv$mrisdp_453 ~ structural_icv$smri_vol_scs_intracranialv)
#  structural_icv[[paste0(roi, "_resid")]] <- resid(model)
#}

#roi='mrisdp_453_resid'
###
for (roi in roi_list[152:216]){
  #roi_res= paste0(roi,"_resid",sep ="")
  ### Merge ROI with fix factors and drops na's
  BOLD=subset(structural_dk, select = c("VisitID", roi))
  BOLD <- na.omit(BOLD)
  df0 <- merge(df_cov,fix_factor, by='VisitID')
  df1<-merge(BOLD,df0, by='VisitID')
  df2 <- df1[df1$sex!='NaN',]
  df3<- df2[df2$age!='NaN',]
  # df4 <- df3[df3$cbcl_scr_syn_attention_r!='NaN',]
  #df4 <- df3[df3$S8!='NaN',]
  df_clean <- df3[df3$age!=0,]
  df_clean <- na.omit(df_clean)
  
  length(unique(df_clean$PGUID))
  # Ensure that 'PGUID', 'Sex', 'scanner'and 'S_S' are factors
  
  df_clean$PGUID <- as.factor(df_clean$PGUID)
  df_clean$sex <- as.factor(df_clean$sex)
  df_clean$S_S <- as.factor(df_clean$S_S)

  
  # Make sure continuous variables are numeric
  df_clean$ICV <-  as.numeric(df_clean$smri_vol_scs_intracranialv)
  df_clean$age <- as.numeric(df_clean$age)
  #df_clean[,2] <- as.numeric(df_clean[,2] )
  #df_clean$S8 <- as.numeric(df_clean$S8)
  #df_clean$C1 <- as.numeric(df_clean$C1)
  #df_clean$C2 <- as.numeric(df_clean$C2)
  #df_clean$C3 <- as.numeric(df_clean$C3)
  #df_clean$C4 <- as.numeric(df_clean$C4)
  df_clean$ROI <- as.numeric(df_clean[,2])
  df_clean$SHN <- as.numeric(df_clean$SHN)
  df_clean$CBCL <- as.numeric(df_clean$cbcl_scr_syn_internal_r)

  # Run the GAMM using gamm4
  model1 <- gamm4(formula = ROI ~ s(age,SHN)+sex+ICV,   random=~(1|S_S/PGUID),data = df_clean)
  #model2 <- gamm4(formula = ROI ~ s(age, by=sex)+sex+ICV,   random=~(1|S_S/PGUID),data = df_clean)
  k_check<-k.check(model1$gam)  
  # Summary of the model
  sum_model=summary(model1$gam)
  
  # ICC
  subject_variance <- as.numeric(VarCorr(model1$mer)$PGUID[1, 1])
  
  # Residual variance (within-subject variance)
  residual_variance <- attr(VarCorr(model1$mer), "sc")^2
  
  # Calculate ICC
  ICC <- subject_variance / (subject_variance + residual_variance)
  
  # Predictions with ggplot2
  predictions = predict(model1$gam, newdata=df_clean, se.fit = TRUE)
  newdat = cbind(df_clean, predictions)
  #newdat = within(newdat, {
  #     lower = fit-1.96*se.fit
  #     upper = fit+1.96*se.fit
  # })
  
  # Consolidating new data and predictions
  label=label_map[label_map$var_name==roi,]
  
  
  measure = strsplit(label$var_label, split="in ")[[1]][1]
  measure = gsub("\\s","_",measure)
  region  = strsplit(label$var_label, split="ROI ")[[1]][2]
  region = gsub("\\s","_",region)
  print(measure)
  print(region)
  #region="Total Surface Area(mm^2)"
  p<-ggplot(newdat, aes(x = age, y = fit))  + geom_smooth() +
    scale_x_continuous(
      breaks = seq(96, 205, by = 12),  # Breaks at every 12 months
      labels = seq(8, 17, by = 1)       # Labels the x-axis as 9, 10, ..., 17
    )+
    #scale_y_continuous(limits=c(170000,200000)
    #)+
    labs(x = "Age", y =" Cortical Thickness", title = "Effect of Age on ROI") + theme_minimal()
  dir_path=paste("/gpfs2/scratch/melmarsr/DK/CT_DK_",measure,visit, "_SHNinteraction/",sep="")
  dir.create(dir_path)
  plot_data<-data.frame(layer_data(p,1))
  plot_data_name=paste("/gpfs2/scratch/melmarsr/DK/CT_DK_",measure,visit, "_SHNinteraction/", measure,region, ".csv", sep="")
  write.csv(plot_data, file=plot_data_name,row.names = FALSE)
  file_name <- paste("/gpfs2/scratch/melmarsr/DK/CT_DK_",measure,visit,"_SHNinteraction/", measure,region, ".png", sep="")
  ggsave(filename= file_name, plot=p, width =6, height = 4)
  
  print(paste("Plot saved as", file_name))
  sum_file_name=paste("/gpfs2/scratch/melmarsr/DK/CT_DK_",measure,visit,"_SHNinteraction/", measure,region, ".txt", sep="")
  capture.output(sum_model,file=sum_file_name,k_check,ICC)
  predict_data= data.frame(newdat$age,newdat$fit)
  predict_file_name=paste("/gpfs2/scratch/melmarsr/DK/CT_DK_",measure,visit,"_SHNinteraction/predicted_", measure,region, ".csv", sep="")
  write.csv(predict_data, file=predict_file_name,row.names = FALSE)
}
