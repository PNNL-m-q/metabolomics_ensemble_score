library(tidyverse)
library(tidymodels)
library(randomForest)

#' Function to apply either the Full Ensemble Model (66 metrics) or the Reduced
#' Ensemble Model (6 metrics) to the output of CoreMS. Column names must have the
#' period in them (see read.csv and read.table).
#' 
#' @param corems_output The output csv/txt files from CoreMS read as a data.frame
#' @param model_type The input model as either "full" or "reduced". Default is "full"
#' 
#' @details The output will be a data.frame with the probability that each spectra
#'.   is either a match (1) or not a match (0)
#'    
#' @example{
#'    setwd("~/Git_Repos/metabolomics_ensemble_score/Models") 
#'    corems_example <- read.table("ExampleCoreMS.txt", sep = "\t", header = T)
#'    apply_model(corems_example)
#'    apply_model(corems_example, model_type = "reduced")
#' }

apply_model <- function(corems_output, model_type = "full", type = "response") {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Ensure corems_output is a data.frame
  if (!is.data.frame(corems_output)) {
    return("corems_output should be a data.frame")
  }
  
  # Model must be either "full" or "reduced"
  if (model_type != "full" & model_type != "reduced") {
    return("model must be either 'full' or 'reduced'")
  }
  
  ################
  ## PULL MODEL ##
  ################
  
  # Extract the model
  if (model_type == "full") {
    model <- readRDS("model.RDS")
  } else {
    model <- readRDS("reduced_model.RDS")
  }
  
  #############
  ## PREDICT ##
  #############
  
  return(predict(model, new_data = corems_example, type = "prob"))
  
}
