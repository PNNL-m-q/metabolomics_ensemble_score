library(tidyverse)

################################################################################
## FULL MODEL PREDICTIONS------------------------------------------------------- 
################################################################################

setwd("~/Git_Repos/refined_ss_score/")

# Load test data 
test_data <- fread("~/Downloads/MQ_RF/Test_Sum.txt") 

# Load test predictions
test_pred <- readRDS("Result_Data/test_pred.RDS")

# Add labels
test_pred <- test_pred %>%
  mutate(
    Label = "",
    Label = ifelse(pred_class == Truth.Annotation & pred_class == 1, "True Positive", Label),
    Label = ifelse(pred_class == Truth.Annotation & pred_class == 0, "True Negative", Label),
    Label = ifelse(pred_class == 0 & Truth.Annotation == 1, "False Negative", Label),
    Label = ifelse(pred_class == 1 & Truth.Annotation == 0, "False Positive", Label)
  )

# Double check: 27 FN, 14501 TN, 2443 TP, 558 FP
table(test_pred$Label)

# Add labels and prediction accuracy
test_data$Label <- test_pred$Label
test_data$Prob <- test_pred$.pred_1
test_data$Bin <- paste0(test_data$`Sample name`, " & ", test_data$`Peak Index`)

# Select Bins, Labels, Probabilities, and calculate ranks
test_ranks <- test_data %>%
  select(Bin, Label, Prob, `Compound Name`) %>% 
  group_by(Bin) %>%
  mutate(Rank = rank(-1 * Prob, ties.method = "min"),
         Model = "Full Model")

################################################################################
## REDUCED MODEL PREDICTIONS----------------------------------------------------
################################################################################

# Load test data 
test_data <- fread("~/Downloads/MQ_RF/Test_Sum.txt") 

# Load reduced
red_pred <- readRDS("Result_Data/reduced_test_pred.RDS")

# Add labels
red_pred <- red_pred %>%
  mutate(
    Label = "",
    Label = ifelse(pred_class == Truth.Annotation & pred_class == 1, "True Positive", Label),
    Label = ifelse(pred_class == Truth.Annotation & pred_class == 0, "True Negative", Label),
    Label = ifelse(pred_class == 0 & Truth.Annotation == 1, "False Negative", Label),
    Label = ifelse(pred_class == 1 & Truth.Annotation == 0, "False Positive", Label)
  )

# Double check: 41 FN, 14404 TN, 2429 TP, 655 FP
table(red_pred$Label)

# Add labels and prediction accuracy
test_data$Label <- red_pred$Label
test_data$Prob <- red_pred$.pred_1
test_data$Bin <- paste0(test_data$`Sample name`,  " & ", test_data$`Peak Index`)

# Select Bins, Labels, Probabilities, and calculate ranks
red_ranks <- test_data %>%
  select(Bin, Label, Prob, `Compound Name`) %>% 
  group_by(Bin) %>%
  mutate(Rank = rank(-1 * Prob, ties.method = "min"),
         Model = "Reduced Model")

################################################################################
## BIND TEST AND REDUCED MODEL RANKS FOR FP and FN------------------------------
################################################################################

FP_FN_Ranks <- rbind(test_ranks, red_ranks)
fwrite(FP_FN_Ranks, "Result_Data/FP_FN_Ranks.txt", quote = F, row.names = F, sep = "\t")

################################################################################
## GET BIN SIZES----------------------------------------------------------------
################################################################################

# Set the working directory
# setwd()

# List all files
files <- list.files(".")

bin_sizes <- do.call(rbind, lapply(files, function(file) {
  
  # Read data 
  fread(file) %>%
    mutate(Bin = paste0(`Sample name`, " & ", `Peak Index`)) %>%
    group_by(Bin) %>%
    summarise(Count = n())
  
}))

fwrite(bin_sizes, "~/Git_Repos/refined_ss_score/Result_Data/BinSizes.csv", quote = F, row.names = F)






