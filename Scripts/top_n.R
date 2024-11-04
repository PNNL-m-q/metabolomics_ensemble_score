# Calculate the top N 
library(tidyverse)
library(data.table)
library(randomForest)

# Set working directory
setwd("~/Git_Repos/refined_ss_score/")

# List all files
files <- list.files("/Users/degn400/Downloads/Dataset33302.tar",
                    full.names = T) %>%
  .[grepl("Sum", .)]

# Parse through all data, pulling out all bins with true positives 
prebin <- do.call(rbind, lapply(files, function(file) {
  
  message(file)
  
  # Read data 
  data <- fread(file) %>% 
    data.frame() %>% 
    na.omit() %>%
    mutate(Bin = paste(Sample.name, "&", Peak.Index)) %>%
    select(-c(Sample.name, Peak.Index))
  
  return(data)
  
}))

# Get all true positives
true_counts <- prebin %>%
  filter(Truth.Annotation == "True.Positive") 

# Filter down to just bins
bin <- prebin %>% filter(Bin %in% true_counts$Bin)

# Run models to determine values 
model <- readRDS("Models/model.RDS")
red_model <- readRDS("Models/reduced_model.RDS")

# Run model
predictions <- predict(model, bin, type = "prob")
saveRDS(predictions, "~/Downloads/MQ_RF/predictions.RDS")
predictions <- readRDS("~/Downloads/MQ_RF/predictions.RDS")

# Run reduced model 
reduced_predictions <- predict(red_model, bin %>% select(Stein.Scott.Similarity.Nist, DFT.Correlation,
                                     DWT.Correlation, Stein.Scott.Similarity, 
                                     Weighted.Cosine.Correlation, Inner.Product.Distance),
                               type = "prob")
saveRDS(reduced_predictions, "~/Downloads/MQ_RF/reduced_predictions.RDS")
reduced_predictions <- readRDS("~/Downloads/MQ_RF/reduced_predictions.RDS")


# 1. Filter down pre-bins to bins with true positives
# 2. Make ranks. Ties should be rated with min, for example, 0.9, 0.5, 0.5, 0.1 should be 1, 2, 2, 4
Ranks <- bin %>%
  mutate(
    Full.Model = predictions$.pred_1,
    Reduced.Model = reduced_predictions$.pred_1
  ) %>%
  select(Bin, Full.Model, Reduced.Model, Stein.Scott.Similarity, Stein.Scott.Similarity.Nist,
         DFT.Correlation, DWT.Correlation, Weighted.Cosine.Correlation, Inner.Product.Distance, Truth.Annotation) %>%
  group_by(Bin) %>%
  mutate(
    Full.Model = rank(-1 * Full.Model, ties.method = "min"),
    Reduced.Model = rank(-1 * Reduced.Model, ties.method = "min"),
    Stein.Scott.Similarity = rank(-1 * Stein.Scott.Similarity, ties.method = "min"),
    Stein.Scott.Similarity.Nist = rank(-1 * Stein.Scott.Similarity.Nist, ties.method = "min"),
    DFT.Correlation = rank(-1 * DFT.Correlation, ties.method = "min"),
    DWT.Correlation = rank(-1 * DWT.Correlation, ties.method = "min"),
    Weighted.Cosine.Correlation = rank(-1 * Weighted.Cosine.Correlation, ties.method = "min"),
    Inner.Product.Distance = rank(Inner.Product.Distance, ties.method = "min")
  ) 

fwrite(Ranks %>% filter(Truth.Annotation == "True.Positive"), "Result_Data/TP_Ranks.txt", sep = "\t", quote = F, row.names = F)








