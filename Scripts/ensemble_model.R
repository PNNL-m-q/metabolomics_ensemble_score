library(tidyverse) # Load data maintenance packages
library(tidymodels) # Load statistical learning packages 
library(parallel) # Enable parallel computing in R
library(doParallel) # Enable parallel computing in R
library(randomForest) # Load the engine
library(data.table) # Fast file reader and write

seed_num <- 4656
set.seed(seed_num)

##########################
## EXPLORATORY ANALYSIS ##
##########################

# Load all data
Sum <- fread("~/Downloads/MQ_RF/All_Sum.txt")

# Pull sample types
SampleMetadata <- fread("~/Git_Repos/refined_ss_score/Metadata/Sample_Metadata.csv")

# Get counts per bin
BinCounts <- Sum %>%
  group_by(ID) %>%
  summarize(
    `TP Counts` = sum(Truth.Annotation == "True.Positive"),
    `TN Counts` = sum(Truth.Annotation == "True.Negative")
  ) %>%
  ungroup() %>%
  mutate(`Sample Name` = map_chr(ID, function(x) {strsplit(x, " & ") %>% unlist() %>% head(1)})) %>%
  left_join(SampleMetadata %>% select(`Sample Name`, `Sample Type`), by = "Sample Name")

# See if we lose any representation of a sample type, which we don't 
BinCounts %>%
  mutate(Fate = ifelse(`TN Counts` == 0, "Lost", "Retained")) %>%
  select(`Sample Type`, Fate) %>%
  ggplot(aes(x = `Sample Type`, fill = Fate)) + geom_bar() + ylab("Count") + theme_bw() 
 
# Remove all bins with no true negatives
Sum <- Sum %>%
  filter(ID %in% unlist(BinCounts[BinCounts$`TN Counts` >= 1, "ID"]))

###################
## INITIAL SPLIT ##
###################

# Split into 75:25 training:test 
toSplit <- unique(Sum$ID)
train_split <- sample(toSplit, round(0.75*length(toSplit)), replace = FALSE)
test_split <- toSplit[toSplit %in% train_split == FALSE]
fwrite(Sum %>% filter(ID %in% train_split), "~/Downloads/MQ_RF/Train_Sum.txt", quote = F, row.names = F, sep = "\t")
fwrite(Sum %>% filter(ID %in% test_split), "~/Downloads/MQ_RF/Test_Sum.txt", quote = F, row.names = F, sep = "\t")

####################################
## BUILD MODEL WITH TRAINING DATA ##
####################################

# Load training data 
train <- fread("~/Downloads/MQ_RF/Train_Sum.txt")

# Load score metadata
ScoreMetadata <- fread("~/Git_Repos/refined_ss_score/Metadata/Score_Metadata.txt")

# Pull true positives 
TPs <- train[train$Truth.Annotation == "True.Positive"]

# Randomly select one true negative per bin
TNs <- train %>% 
  filter(Truth.Annotation == "True.Negative") %>%
  group_by(ID) %>%
  slice_sample(n = 1)

# Create dataset
toModel <- rbind(TPs, TNs) %>%
  dplyr::select(Truth.Annotation, ScoreMetadata$Score) %>%
  dplyr::mutate(Truth.Annotation = factor(ifelse(Truth.Annotation == "True.Positive", "1", "0"), levels = c("0", "1"))) %>%
  data.frame()
toModel[,33] <- as.numeric(toModel[,33])
  
# Create a recipe 
my_recipe <- recipe(Truth.Annotation ~ ., data = toModel)

# Designate a resampling scheme. Let's do 4 fold validation 10 times 
folds <- vfold_cv(toModel,
                  v = 4,
                  repeats = 5,
                  strata = Truth.Annotation,
                  breaks = 4)

# Designate modeling engine. Let's allow for model tuning of trees and variables at each split. 
rf_engine <- rand_forest(trees = tune(), mtry = tune()) %>%
  set_engine("randomForest") %>%
  set_mode("classification")

# Designate test metrics 
test_metrics <- metric_set(accuracy)

# Specify number of trees based on a number of levels and m-try
rf_grid <- grid_regular(trees(), mtry(range = c(1, ncol(toModel)-1)), levels = 3) 

# Define the workflow
wf <- workflow_set(
  preproc = list("pre" = my_recipe),
  models = list("rf" = rf_engine)
)

# Pull the workflow ID
wfid <- wf$wflow_id

# Add the grid to the specific workflow id
train_wfs <- wf %>% option_add(id = wfid, grid = rf_grid)

# Designate half of the cores to run 
ncores <- floor(parallel::detectCores(logical = TRUE)/2)

# Set the number of cores
cl <- parallel::makePSOCKcluster(ncores)

# Register cores to use
doParallel::registerDoParallel(cl)

# Run tuning - will take 1-5 minutes
wf_complete <- train_wfs %>%
  workflowsets::workflow_map(resamples = folds,
                             verbose = TRUE, 
                             seed = seed_num,
                             metrics = test_metrics)

saveRDS(wf_complete, "~/Downloads/MQ_RF/wf_complete.RDS")
wf_complete <- readRDS("~/Downloads/MQ_RF/wf_complete.RDS")

# Stop parallel computing cores 
stopCluster(cl)

# Extract best fit
best_rf <- wf_complete %>%
  extract_workflow_set_result(id = "pre_rf") %>%
  select_best(metric = "accuracy")
best_rf

# Fit to training
train_fitted <- wf_complete %>%
  extract_workflow(id = "pre_rf") %>%
  finalize_workflow(best_rf) %>%
  fit(data = toModel)
train_fitted

saveRDS(train_fitted, "~/Downloads/MQ_RF/train_fitted.RDS")

# Pull testing example
test <- fread("~/Downloads/MQ_RF/Test_Sum.txt")  %>%
  dplyr::select(Truth.Annotation, ScoreMetadata$Score) %>%
  dplyr::mutate(Truth.Annotation = factor(ifelse(Truth.Annotation == "True.Positive", "1", "0"), levels = c("0", "1"))) %>%
  data.frame()
test[,33] <- as.numeric(test[,33])
table(test$Truth.Annotation) # 15059 true negatives, 2470 true positives 

# Get prediction matrix
test_pred <- predict(train_fitted, test) %>%
  bind_cols(predict(train_fitted, test, type = "prob")) %>%
  bind_cols(test %>% select(Truth.Annotation)) %>%
  rename(pred_class = .pred_class)
test_pred

saveRDS(test_pred, "~/Downloads/MQ_RF/test_pred.RDS")

# Pull model object
model <- train_fitted %>%
  extract_fit_parsnip() 

saveRDS(model, "~/Downloads/MQ_RF/model.RDS")
  
###############################
## REDUCED MODEL PERFORMANCE ##
###############################

# Set reduced variables 
reduced <- c("Truth.Annotation", "Stein.Scott.Similarity.Nist", "DFT.Correlation",
             "DWT.Correlation", "Stein.Scott.Similarity", "Weighted.Cosine.Correlation",
             "Inner.Product.Distance")
reducedModel <- toModel %>% select(all_of(reduced))

# Create a recipe 
red_recipe <- recipe(Truth.Annotation ~ ., data = reducedModel)

# Designate a resampling scheme. Let's do 4 fold validation 10 times 
red_folds <- vfold_cv(reducedModel,
                  v = 4,
                  repeats = 5,
                  strata = Truth.Annotation,
                  breaks = 4)

# Designate modeling engine. Let's allow for model tuning of trees and variables at each split. 
red_rf_engine <- rand_forest(trees = tune(), mtry = tune()) %>%
  set_engine("randomForest") %>%
  set_mode("classification")

# Designate test metrics 
test_metrics <- metric_set(accuracy)

# Specify number of trees based on a number of levels and m-try
red_rf_grid <- grid_regular(trees(), mtry(range = c(1, ncol(reducedModel)-1)), levels = 3) 

# Define the workflow
red_wf <- workflow_set(
  preproc = list("pre" = red_recipe),
  models = list("rf" = red_rf_engine)
)

# Pull the workflow ID
red_wfid <- red_wf$wflow_id

# Add the grid to the specific workflow id
red_train_wfs <- red_wf %>% option_add(id = red_wfid, grid = red_rf_grid)

# Run tuning - will take 1-5 minutes
red_complete <- red_train_wfs %>%
  workflowsets::workflow_map(resamples = red_folds,
                             verbose = TRUE, 
                             seed = seed_num,
                             metrics = test_metrics)

saveRDS(red_complete, "~/Downloads/MQ_RF/red_complete.RDS")
red_complete <- readRDS("~/Downloads/MQ_RF/red_complete.RDS")

# Extract best fit
red_best_rf <- red_complete %>%
  extract_workflow_set_result(id = "pre_rf") %>%
  select_best(metric = "accuracy")
red_best_rf

# Fit to training
red_train_fitted <- red_complete %>%
  extract_workflow(id = "pre_rf") %>%
  finalize_workflow(red_best_rf) %>%
  fit(data = reducedModel)
red_train_fitted

saveRDS(red_train_fitted, "~/Downloads/MQ_RF/red_train_fitted.RDS")

# Pull testing example
red_test <- fread("~/Downloads/MQ_RF/Test_Sum.txt")  %>%
  dplyr::select(Truth.Annotation, ScoreMetadata$Score) %>%
  dplyr::mutate(Truth.Annotation = factor(ifelse(Truth.Annotation == "True.Positive", "1", "0"), levels = c("0", "1"))) %>%
  data.frame() %>%
  select(reduced)
table(red_test$Truth.Annotation) # 15059 true negatives, 2470 true positives 

# Get prediction matrix
red_test_pred <- predict(red_train_fitted, red_test) %>%
  bind_cols(predict(red_train_fitted, red_test, type = "prob")) %>%
  bind_cols(red_test %>% select(Truth.Annotation)) %>%
  rename(pred_class = .pred_class)
red_test_pred

saveRDS(red_test_pred, "~/Downloads/MQ_RF/red_test_pred.RDS")

# Pull model object
red_model <- red_train_fitted %>%
  extract_fit_parsnip() 

saveRDS(red_model, "~/Downloads/MQ_RF/red_model.RDS")






