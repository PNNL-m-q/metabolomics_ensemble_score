library(tidyverse) # Load data maintenance packages
library(tidymodels) # Load statistical learning packages 
library(randomForest) # Load the engine
library(data.table) # Fast file reader and write
library(patchwork) # Arrange plots

# Set working directory
setwd("~/Git_Repos/refined_ss_score/")

# Set plot theme 
theme_set(theme_bw())

################################################################################
## Load results-----------------------------------------------------------------
################################################################################

# Load full workflow 
model <- readRDS("Models/model.RDS")

# Load reduced model
red_model <- readRDS("Models/reduced_model.RDS")

# Load prediction matrix
test_pred <- readRDS("Result_Data/test_pred.RDS")

# Load reduced prediction matrix
red_pred <- readRDS("Result_Data/reduced_test_pred.RDS")

# Load score metadata
ScoreMetadata <- fread("Metadata/Score_Metadata.txt")

# Merge Mean Decrease in Gini Index to ScoreMetadata and compare 
Imp <- importance(model$fit) %>% 
  data.frame() %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(
    Score = gsub(pattern = ".", replacement = " ", row.names(.), fixed = T),
    Score = ifelse(Score == "Avg  L1  L8  Distance", "Avg (L1, L8) Distance", Score),
    Score = ifelse(Score == "Squared chord Distance", "Squared-chord Distance", Score)
  ) %>%
  left_join(ScoreMetadata) %>%
  mutate(LogGini = log(MeanDecreaseGini))

################################################################################
## FIGURE 1: FULL MODEL PERFORMANCE--------------------------------------------- 
################################################################################

conf_mat(test_pred, truth = Truth.Annotation, estimate = pred_class) %>%
  autoplot(type = "heatmap")

# Confusion matrix 
Fig1a <- data.frame(
  Truth = factor(c("0", "1", "0", "1"), levels = c("1", "0")),
  Prediction = factor(c("1", "1", "0", "0"), levels = c("0", "1")),
  Count = c(558, 2443, 14501, 27),
  Labels = c("False Positives:\n558", "True Positives:\n2443", "True Negatives:\n14501", "False Negatives:\n27")
) %>%
  ggplot(aes(x = Truth, y = Prediction, fill = Count, label = Labels)) + 
  geom_tile() +
  geom_text() +
  ggtitle("Confusion Matrix") +
  scale_fill_gradient2(low = "white", high = "steelblue") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

TP <- 2443
FP <- 558
TN <- 14501
FN <- 27

# Metric table
data.frame(
  Metric = c("TPR", "TNR"),
  Value = c(TP / (TP + FN), TN / (TN + FP))
) %>% 
  pivot_wider(names_from = Metric, values_from = Value) %>%
  mutate(
    BA = (TPR + TNR) / 2,
    Recall = TP / (TP + FN),
    Precision = TP / (TP + FP),
    F1 = (2 * Precision * Recall) / (Precision + Recall),
    AUC = 0.996 # roc_auc(test_pred, truth = Truth.Annotation, .pred_0)
  ) %>%
  pivot_longer(cols = 1:6) %>%
  rename(Metric = name, Value = value) %>% knitr::kable()

Fig1b <- roc_curve(test_pred, truth = Truth.Annotation, .pred_0) %>%
  rename(Sensitivity = sensitivity, Specificity = specificity) %>%
  ggplot(aes(x = 1 - Specificity, y = Sensitivity)) + 
  geom_line() +
  ggtitle("ROC Curve, AUC = 0.996") +
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") +
  theme(plot.title = element_text(hjust = 0.5))

Fig1 <- Fig1a + Fig1b + plot_annotation(tag_levels = "A")
Fig1

################################################################################
## FIGURE 2: Scree plots and visualization of performance per metadata variable
################################################################################

## Scree plot ##

# Scree plot 
F2a <- data.table(
  `Gini Index` = Imp$MeanDecreaseGini,
  `Top Score` = factor(c(rep("Yes", 6), rep("No", 60)), levels = c("Yes", "No")),
  Score = 1:nrow(ScoreMetadata)
) %>%
  ggplot(aes(x = Score, y = `Gini Index`, color = `Top Score`)) + 
  geom_point(size = 2) +
  theme(plot.title = element_text(hjust = 0.5, size = 16), legend.text = element_text(size = 12),
        legend.position = "right",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_text(size = 14)) +
  scale_color_manual(values = c("red", "black")) +
  ggtitle("Scree Plot")
F2a

## Importances ##
boxplot_fun <- function(xvar, colorvar, title, fill_pallette = NULL) {
  plot <- data.frame(
    XVar = xvar,
    Colorvar = colorvar,
    LogGini = Imp$LogGini
  ) %>%
    ggplot(aes(x = XVar, y = LogGini)) + 
    geom_boxplot(aes(fill = XVar), outlier.shape = NA) +
    geom_jitter(aes(color = Colorvar), width = 0.1, height = 0, na.rm = T) +
    xlab("") +
    ggtitle(title) +
    ylab("Log Gini Index") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    scale_color_manual(values = c("red", "black"))
  if (is.null(fill_pallette)) {return(plot)} else {
    return(plot + scale_fill_brewer(palette = fill_pallette))
  }
}

# Add Top Score to importance
Imp <- Imp %>%
  mutate(TopScore = ifelse(Score %in% c("Stein Scott Similarity Nist", "Stein Scott Similarity",
                                        "DFT Correlation", "DWT Correlation", "Weighted Cosine Correlation", 
                                        "Inner Product Distance"), "Yes", "No"),
         TopScore = factor(TopScore, levels = c("Yes", "No")))
  
# Importances per metric type 
F2b <- boxplot_fun(factor(Imp$Type, levels = c("Similarity", "Distance")), Imp$TopScore, "Score Type", "Blues")
F2b

# Importances per metric family
familyorder <- Imp %>% 
  select(LogGini, Family) %>% 
  group_by(Family) %>% 
  summarise(Median = median(LogGini, na.rm = T)) %>% 
  arrange(-Median) %>% 
  select(Family) %>% 
  unlist()

F2c <- boxplot_fun(factor(Imp$Family, levels = familyorder), Imp$TopScore, "Family", "Greys") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Put plots together
(F2a + F2b) / F2c + plot_annotation(tag_levels = "A")


################################################################################
## FIGURE 3: TP Rankings--------------------------------------------------------
################################################################################

# Read in score ranks 
TPRanks <- fread("Result_Data/TP_Ranks.txt")

# Top 1 and Top 5
TPRanks %>%
  select(2:9) %>%
  pivot_longer(1:8) %>%
  rename(Model = name, Rank = value) %>%
  group_by(Model) %>%
  summarise(
    `Top 1 Proportion` = sum(Rank == 1) / nrow(TPRanks),
    `Top 5 Proportion` = sum(Rank <= 5) / nrow(TPRanks) 
  ) %>%
  arrange(-`Top 1 Proportion`) %>%
  mutate(Model = gsub(".", " ", Model, fixed = T),
         Model = factor(Model, levels = Model)) %>%
  pivot_longer(2:3) %>%
  rename(Rank = name, Proportion = value) %>%
  mutate(Proportion = round(Proportion, 3)) %>%
  ggplot(aes(x = Model, y = Proportion)) +
    geom_bar(stat = "identity") + 
    geom_text(aes(label = Proportion, vjust = -0.3)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ylab("Proportion of True Positive Bins") +
    xlab("") +
    ylim(c(0, 1.05)) +
    facet_wrap(.~Rank)

################################################################################
## Figure 4: FALSE POSITIVE AND FALSE NEGATIVE INVESTIGATION--------------------
################################################################################

F_Ranks <- left_join(fread("Result_Data/FP_FN_Ranks.txt"), fread("Result_Data/BinSizes.csv"), by = "Bin")

# Bin Size----------------------------------------------------------------------

F4a <- F_Ranks %>%
  select(Label, Model, Count) %>%
  mutate(Label = factor(Label, levels = c("True Positive", "False Positive", "True Negative", "False Negative"))) %>%
  rename(Truth = Label) %>%
  ggplot(aes(x = Truth, y = Count, fill = Truth)) +
    geom_boxplot() + 
    scale_fill_manual(values = c("steelblue4", "orchid4", "firebrick4", "orange4")) +  
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    ylim(c(0, 105)) + 
    xlab("") +
    ylab("Number of \nCandidate Matches") +
    facet_wrap(.~Model, scales = "free")

F_Ranks %>%
  select(Label, Model, Count) %>%
  mutate(Label = factor(Label, levels = c("True Positive", "False Positive", "True Negative", "False Negative"))) %>%
  rename(Truth = Label) %>%
  group_by(Truth, Model) %>%
  summarise(Median = median(Count))


# Metabolite counts-------------------------------------------------------------

CM <- fread("Metadata/Compound_Metadata.txt")


F4c <- F_Ranks %>%
  group_by(Model, Label, `Compound Name`) %>%
  summarise(Freq = n()) %>%
  ungroup() %>%
  group_by(Model, Label) %>%
  slice_max(order_by = Freq, n = 5) %>%
  mutate(
    Label = factor(Label, levels = c("True Positive", "False Positive", "True Negative", "False Negative")),
    `Compound Name` = gsub("\\[.+?\\]", "", `Compound Name`) %>% trimws(which = "both")
  ) %>%
  left_join(CM, by = "Compound Name") %>%
  arrange(-Freq) %>%
  mutate(`Compound Name` = factor(`Compound Name`, levels = unique(`Compound Name`))) %>%
  ggplot(aes(x = `Compound Name`, y = Freq, fill = `Compound Type`)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = Freq, vjust = -0.1)) +
    ylim(c(0, 400)) +
    ylab("Counts of Top 5 Compounds") + 
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_grid(rows = vars(Model), cols = vars(Label), scales = "free") 

# Sample type counts------------------------------------------------------------

SampleMetadata <- fread("Metadata/Sample_Metadata.csv") %>%
  select(`Sample Name`, `Sample Type`) 

F_Ranks <- F_Ranks %>%
  mutate(`Sample Name` = map_chr(Bin, function(x) {strsplit(x, " & ") %>% unlist() %>% head(1)})) %>%
  left_join(SampleMetadata, by = "Sample Name")

F4b <- F_Ranks %>%
  select(Model, Label, `Sample Type`) %>%
  mutate(`Sample Type` = ifelse(`Sample Type` %in% c("Standard Mixture", "Standard"), "Standard", "Non-Standard")) %>%
  group_by(Model, Label, `Sample Type`) %>%
  summarise(`Sample Count` = n()) %>%
  ungroup() %>%
  filter(Label %in% c("True Positive", "False Positive")) %>%
  pivot_wider(id_cols = c(Model, `Sample Type`), values_from = `Sample Count`, names_from = Label) %>%
  mutate(
    `# True Positives / # Positives` = round(`True Positive` / (`True Positive` + `False Positive`), 4)
  ) %>%
  arrange(-`# True Positives / # Positives`) %>%
  ggplot(aes(x = `Sample Type`, y = `# True Positives / # Positives`, fill = `Sample Type`)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = `# True Positives / # Positives`, vjust = -0.1)) +
    scale_fill_manual(values = c("steelblue", "darkorange")) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    xlab("") +
    ylim(c(0, 1.05)) +
    facet_wrap(.~Model)

(F4a + F4b) / F4c + plot_annotation(tag_levels = "A")








