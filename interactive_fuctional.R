#install.packages("gamm4")
library(gamm4)  # for gamm4 function
library(readxl)  # for reading Excel files
library(ggplot2)  # for enhanced plotting
library(readr)
library(dplyr)
library(mgcv)
library(haven)
library(sjPlot)
library(caret)
library(lme4)
library(lmerTest)

# Read data from Excel file
df_cov <- read_excel("/gpfs2/scratch/melmarsr/DK/data/datasetupdate.xlsx", sheet = "MID")
names(df_cov)[names(df_cov) == "age...6"] <- "age"
df_cov=subset(df_cov, select= c("VisitID","QC_DAIRC","PGUID",'age','site','scanner','S_S','sex'))
label_map <- read_csv("/gpfs2/scratch/melmarsr/DK/data/label_map.csv")

# Drop subjects with performance flags
#performance_pguid <- read_excel("/gpfs2/scratch/melmarsr/DK/data/nback_6.0_PGUIDs_QC.xlsx", sheet = "Sheet1")
#names(performance_pguid)[names(performance_pguid) == "5.0 PGUID format"] <- "PGUID"
#df_cov<-merge(performance_pguid,df_cov,by="PGUID")
#df_cov<-filter(df_cov, df_cov$PGUID %in% performance_pguid$PGUID)
df_cov<-na.omit(df_cov)

# Find brain data file
##abcd_midabwdp01--MID
##abcd_tfsstabwdp101--SST
##abcd_tfabwdp101--Nback
task_mapper=na.omit(label_map[label_map$table_name_nda =="abcd_midabwdp01",])
task<- read_csv("/gpfs2/scratch/melmarsr/DK/data/ABCD_sync/abcd_midabwdp01.csv")

# Filter contrasts and runs
filter_mapper= task_mapper %>%  filter(grepl("Beta weight for MID",var_label))
roi_list=filter_mapper$var_name

### Outliers detection
### Winsorization
winsorize_3sd <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  lower_bound <- mean_x - 3 * sd_x
  upper_bound <- mean_x + 3 * sd_x
  x[x < lower_bound] <- lower_bound
  x[x > upper_bound] <- upper_bound
  return(x)
}
### Filter outliers
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

functional_destr_bl <- filter(task, eventname == "baseline_year_1_arm_1")
functional_destr_2y <- filter(task, eventname == "2_year_follow_up_y_arm_1")
functional_destr_4y <- filter(task, eventname == "4_year_follow_up_y_arm_1")
functional_destr_6y <- filter(task, eventname == "6_year_follow_up_y_arm_1")
functional_destr_bl <- filter_outliers(functional_destr_bl[,0:531])
functional_destr_2y <- filter_outliers(functional_destr_2y[,0:531])
functional_destr_4y <- filter_outliers(functional_destr_4y[,0:531])
functional_destr_6y <- filter_outliers(functional_destr_6y[,0:531])

# Apply the outliers filtering function to each column
#df_winsorized <- as.data.frame(lapply(task, winsorize_3sd))
functional_destr <-rbind(functional_destr_bl, functional_destr_2y, functional_destr_4y, functional_destr_6y)
head(functional_destr)

# Check the covariates of your data
str(df_cov)
df_cov= subset(df_cov, select = c("VisitID", "sex", "age","PGUID","S_S",'QC_DAIRC'))
    df_cov <- df_cov[df_cov$QC_DAIRC==1,]

    
# Create a table of counts for each unique item
item_counts <- table(df_cov$PGUID)

# Filter items that appear more than once
visit=1
count=visit-1
filtered_list <- df_cov$PGUID[df_cov$PGUID%in% names(item_counts[item_counts > count])]


df_cov <- df_cov %>% filter(PGUID %in% filtered_list )
str(df_cov)

roi = roi_list[1]

for (roi in roi_list[1:148]){
  # merge covariates and brain data
  BOLD=subset(functional_destr, select = c("VisitID", roi))
  df1<-merge(BOLD,df_cov, by='VisitID')
  df2 <- df1[df1$sex!='NaN',]
  df_clean <- df2[df2$age!='NaN',]
  df_clean <- df_clean[df_clean$age!=0,]
  
  # Ensure that 'PGUID', 'Sex', and 'S_S' are factors
  df_clean$PGUID <- as.factor(df_clean$PGUID)
  df_clean$sex <- as.factor(df_clean$sex)
  df_clean$S_S<- as.factor(df_clean$S_S)
  # Make sure 'Age' is numeric
  df_clean$age <- as.numeric(df_clean$age)
  df_clean$ROI <- as.numeric(df_clean[,2] )
  
  # Run the GAMM using gamm4
  #y.gamm4 <- gam(formula = ROI ~ s(age*sex), random=~(1|S_S/PGUID),data = df_clean,REML=TRUE)
  model1 <- gamm4(formula = ROI ~ s(age)+sex,   random=~(1|S_S/PGUID),data = df_clean,ML=TRUE)
  model2 <- gamm4(formula = ROI ~ s(age, by=sex)+sex,   random=~(1|S_S/PGUID),data = df_clean,ML=TRUE)

  # Summary of the model
  sum=anova(model1$mer,model2$mer)
  # Plot the model
  p = plot_model(model2, type = "pred", terms = c("age", "sex"), axis.title =c('age','ROI'), show.data = FALSE)
  
  # Get the ROI labels
  label=label_map[label_map$var_name==roi,]
  contrast= strsplit(strsplit(label$var_label, split=" contrast in ")[[1]][1], split = "for ")[[1]][2]
  contrast= gsub("\\s","_",contrast)
  region =paste(gsub("\\s", "_", strsplit(strsplit(label$var_label, split=" contrast in ")[[1]][2], split = "cortical Destrieux ROI")[[1]]),collapse="")
  print(contrast)
  print(region)
  
  file_name <- paste("/gpfs2/scratch/melmarsr/DK/",contrast,visit,"_sex_interaction/", contrast,region, ".png", sep="")
  ggsave(filename= file_name, plot=p, width =6, height = 4)
  print(paste("Plot saved as", file_name))
  sum_file_name=paste("/gpfs2/scratch/melmarsr/DK/",contrast,visit,"_sex_interaction/summary_", contrast,region, ".txt", sep="")
  capture.output(paste(summary(model1$gam),summary(model1$gam),sum),file=sum_file_name)
  
}

