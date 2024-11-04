library(data.table)
library(tidyverse)

# Set the working directory
# sewtd()

# List all files
files <- list.files(".")

all_sum <- do.call(rbind, lapply(files, function(file) {
  
  message(file)
  
  # Read data 
  data <- fread(file)
  
  # Make a data ID 
  data$ID <- paste0(data$`Sample name`, " & ", data$`Peak Index`)
  
  # Extract all true positive cases 
  TP_IDs <- data %>% filter(Truth.Annotation == "True.Positive") %>% select(ID) %>% unlist()
  
  # Check that there is no duplicates in the true positive IDs
  counts <- table(TP_IDs) %>% data.frame()
  if (any(counts$Freq > 1)) {message("Stop! There's a true positive count of greater than 1 for a bin.")}
  
  # Filter to true positive bins 
  filtered_data <- data %>% 
    filter(ID %in% TP_IDs & Truth.Annotation %in% c("True.Negative", "True.Positive"))
  
  return(filtered_data)
  
}))


