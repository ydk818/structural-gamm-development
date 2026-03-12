library(gamm4)  #
library(readxl)  
library(ggplot2)  
library(readr)
library(dplyr)
library(mgcv)
library(haven)
library(sjPlot)
library(caret)


### Load covariates and structural data
qc_data <- read_excel("/gpfs2/scratch/melmarsr/DK/data/datasetupdate.xlsx", sheet = "QC_DAIRC")
names(qc_data)[names(qc_data) == "age...10"] <- "age"
names(qc_data)[names(qc_data) == "PGUID...2"] <- "PGUID"
names(qc_data)[names(qc_data) == "S_S"] <- "S_S"
qc_id<- read_csv("/gpfs2/scratch/melmarsr/DK/data/abcd_imgincl01.csv")
qc_index=qc_id[qc_id$imgincl_t1w_include==1,]$VisitID
label_map <- read_csv("/gpfs2/scratch/melmarsr/DK/data/label_map.csv")

### Find labels for structural data
structural_mapper=na.omit(label_map[label_map$table_name_nda =="abcd_smrip102",])
structural_list=structural_mapper$var_name
structural_destr <- read_csv("/gpfs2/scratch/melmarsr/DK/data/ABCD_sync/abcd_smrip102.csv")
filter_mapper= structural_mapper %>%  filter(grepl("Cortical area in mm of APARC ROI ",var_label))
roi_list=filter_mapper$var_name

### Filter the flagged subjects and load covariates 
df_cov= subset(qc_data, select = c("VisitID", "sex", "age","PGUID","S_S"))
df_cov= df_cov %>% filter(VisitID %in% qc_index)

## Winsorization for structural data
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

structural_destr_bl <- filter(structural_destr, eventname == "baseline_year_1_arm_1")
structural_destr_2y <- filter(structural_destr, eventname == "2_year_follow_up_y_arm_1")
structural_destr_4y <- filter(structural_destr, eventname == "4_year_follow_up_y_arm_1")
structural_destr_6y <- filter(structural_destr, eventname == "6_year_follow_up_y_arm_1")
structural_destr_bl <- filter_outliers(structural_destr_bl[,0:68])
structural_destr_2y <- filter_outliers(structural_destr_2y[,0:68])
structural_destr_4y <- filter_outliers(structural_destr_4y[,0:68])
structural_destr_6y <- filter_outliers(structural_destr_6y[,0:68])

# Apply the winsorization function to each column
structural_destr <-rbind(structural_destr_bl, structural_destr_2y, structural_destr_4y, structural_destr_6y)
#structural_destr <- as.data.frame(lapply(structural_destr, winsorize_3sd))

### Load ICV variable
icv_data=read_csv("/gpfs2/scratch/melmarsr/DK/data/ABCD_sync/abcd_smrip102.csv")
icv=subset(icv_data, select = c("VisitID", "smri_vol_scs_intracranialv","subjectkey"))

### Load CBCL data 
phenotype_fix="cbcl_scr_syn_attention_r"
cbcl_data=read_csv("/gpfs2/scratch/melmarsr/DK/data/cbcl_nih_sma_diff.csv")
cbcl_adhd=subset(cbcl_data, select = c("subjectkey", phenotype_fix))

### Load Polygenic Risk Scores
genetic_data= read_sav("/gpfs1/home/d/y/dyuan/matt_gamm/mds_and_MDD_prs.sav")
genetic=subset(genetic_data, select = c("subjectid", "S8","C1","C2","C3","C4"))
names(genetic)[names(genetic) == 'subjectid'] <- 'subjectkey'

### Merge Fix Factors into one Dataframe
fix_factor=merge(genetic,icv, by='subjectkey')

# Create a table of counts for each unique item
item_counts <- table(df_cov$PGUID)

# Filter items that appear more than once
visit=1
count=visit-1
filtered_list <- df_cov$PGUID[df_cov$PGUID%in% names(item_counts[item_counts > count])]
df_cov <- df_cov %>% filter(PGUID %in% filtered_list )
str(df_cov)
unique(df_cov[,4])

### Test one ROI before the loop 

roi= roi_list[0:1]
label_map[label_map$var_name==roi,]
label=label_map[label_map$var_name==roi,]
print(label)
### Different loops for different measures
for (roi in roi_list[0:68]){
  
  ### Merge ROI with fix factors and drops na's
  BOLD=subset(structural_destr, select = c("VisitID", roi))
  df0 <- merge(df_cov,icv, by='VisitID')
  df1<-merge(BOLD,df0, by='VisitID')
  df2 <- df1[df1$sex!='NaN',]
  df3<- df2[df2$age!='NaN',]
  #df4 <- df3[df3$cbcl_scr_syn_attention_r!='NaN',]
  #df4 <- df3[df3$S8!='NaN',]
  df_clean <- df3[df3$age!=0,]
  df_clean <- na.omit(df_clean)
  
  # Ensure that 'PGUID', 'Sex', 'scanner'and 'S_S' are factors
  
  df_clean$PGUID <- as.factor(df_clean$PGUID)
  df_clean$sex <- as.factor(df_clean$sex)
  df_clean$scanner <- as.factor(df_clean$S_S)
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
  
  # Run the GAMM interactive model
  #model<- gam(formula = ROI  ~ s(age) + S8 + ti(age, S8) + sex + ICV + C1 + C2 + C3 + C4)
  model1 <- gamm4(formula = ROI ~ s(age)+sex+ICV,   random=~(1|S_S/PGUID),data = df_clean)
  model2 <- gamm4(formula = ROI ~ s(age, by=sex)+sex+ICV,   random=~(1|S_S/PGUID),data = df_clean)
  
  
  # Summary of the model
  sum=anova(model1$mer,model2$mer)
  # Plot the model
  p = plot_model(model2, type = "pred", terms = c("age", "sex"), axis.title =c('age','ROI'), show.data = FALSE)
  

  # Save the plot as png and summary file as txt
  label=label_map[label_map$var_name==roi,]
  
  measure= strsplit(label$var_label, split="in mm")[[1]][1]
  measure = gsub("\\s","_",measure)
  region= strsplit(label$var_label, split=" ROI ")[[1]][2]
  region = gsub("\\s","_",region)
  
  file_name <- paste("/gpfs2/scratch/melmarsr/DK/Cortical_thickness_",visit,"_sex_interaction/ggplot_", measure,region, ".png", sep="")
  ggsave(filename= file_name, plot=p, width =6, height = 4)
  print(paste("Plot saved as", file_name))
  sum_file_name=paste("/gpfs2/scratch/melmarsr/DK/Cortical_thickness_",visit,"_sex_interaction/ggplot_", measure,region, ".txt", sep="")
  capture.output(sum,file=sum_file_name)
}

