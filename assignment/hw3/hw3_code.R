
# Load necessary libraries
library(dplyr)             # For data manipulation
library(tidyverse)         # For data science workflows
library(ggplot2)           # For creating visualizations and plots
library(caret)             # For machine learning modeling and preprocessing
library(predtools)         # For predictive modeling utilities
library(magrittr)          # For improved piping operations
library(Metrics)           # For computing regression and classification metrics
library(pROC)              # For calculating and plotting AUC-ROC curves
library(PRROC)             # For computing precision-recall AUC (PR-AUC)
library(ROCR)              # For ROC curve performance evaluation
library(xgboost)           # For extreme gradient boosting (XGBoost) models
library(Matrix)            # For handling sparse matrices (used in XGBoost)
library(DataExplorer)      # For automated exploratory data analysis
library(CalibrationCurves) # For generating calibration curves
library(caTools)           # For splitting data and AUC computation
library(MASS)              # For stepwise selection
library(iml)               # For SHAP analysis
library(knitr)             # For generating tables in R Markdown
library(kableExtra)        # For enhanced table formatting options in LaTeX
library(SHAPforxgboost)    # For SHAP interpretation for XGBoost models
library(shapviz)           # For Visualization of SHAP values
library(gridExtra)         # For arranging multiple ggplots in a grid layout

################################################################################
# Database Connection Setup
################################################################################

# Setup for database
db_host <- "spinup-db001ec7.cluster-c9ukc6s0rmbg.us-east-1.rds.amazonaws.com"
db_user <- "introml568"
db_pass <- "m7bxMRtyMqPbcxyRRGML8"
db_name <- "urineculture"

# Connect to database
db_conn <- DBI::dbConnect(RPostgres::Postgres(),
                          host = db_host,
                          dbname = db_name,
                          user = db_user,
                          password = db_pass,
                          options="-c search_path=public")

# Clean up workspace (remove sensitive variables)
rm(db_name, db_host, db_pass, db_user)

################################################################################
# Data Collection and Processing
################################################################################

# Load data from database
raw_data <- tbl(db_conn, 'results') %>% dplyr::collect()

# Drop irrelevant columns
df <- raw_data %>% dplyr::select(-UCX_abnormal, -split, -alt_diag, -abxUTI)

# Inspect dataset
# str(df)
# summary(df)
# head(df)
# if (ncol(df) >= 216) {
#   head(df[, 210:216])
# }

# Observations:
# - NA values present
# - Inconsistent data types
# - Logical values stored as integers

# Check for redundant variables, i.e. only have NA values
# which(apply(df, 2, function(x) sum(is.na(x))) == nrow(df))

# Remove redundant variables
df <- df[, -which(apply(df, 2, function(x) sum(is.na(x))) == nrow(df))]

# Check for only one response
# which(apply(df, 2, function(x) length(unique(x))) == 1)

# Check for not on likert
non_likert <- apply(df, 2, function(x) length(unique(x))) > 6

# Check for numeric
num <- sapply(df, is.numeric)

# Check for binary 
binary <- apply(df, 2, function(x) sum(c('0', '1') %in%
                                         unique(x)) == length(unique(x)))

# Convert data types accordingly
df[, (non_likert & num)] <- sapply(df[, (non_likert & num)], as.numeric)
df[, !(non_likert & num)] <- sapply(df[, !(non_likert & num)], as.character)
df[, binary] <- sapply(df[, binary], as.integer)
df[, binary] <- sapply(df[, binary], as.logical)

# Drop identifier column
df$id <- NULL

# Clear temporary variables
rm(num, binary, non_likert, raw_data, db_conn)

################################################################################
# Functions created and used
################################################################################

replace_na_with_nr <- function(df) {
  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], 
                                         function(x) ifelse(is.na(x), 'not_reported', x))
  return(df)
}

evaluate_glm <- function(model, train, test) {
  
  # Set up training results
  train_pred <- predict(model, type = 'response')
  train_pred_binary <- ifelse(train_pred > 0.5, 1, 0)
  train_results <- data.frame(y = train$UTI_diag, 
                              pred = train_pred, 
                              pred_binary = train_pred_binary)
  
  # Set up test results
  test_pred <- predict(model, test, type = 'response')
  test_pred_binary <- ifelse(test_pred > 0.5, 1, 0)
  test_results <- data.frame(y = test$UTI_diag, 
                             pred = test_pred, 
                             pred_binary = test_pred_binary)
  
  # Calibration plots
  p1 <- calibration_plot(data = train_results, obs = "y", pred = "pred", 
                         title = "Calibration Plot - Training")
  p2 <- calibration_plot(data = test_results, obs = "y", pred = "pred", 
                         title = "Calibration Plot - Validation")
  
  # MSE and RMSE
  mse_train <- mse(train_results$y, train_results$pred)
  rmse_train <- rmse(train_results$y, train_results$pred)
  mse_test <- mse(test_results$y, test_results$pred)
  rmse_test <- rmse(test_results$y, test_results$pred)
  mses <- data.frame(MSE = c(mse_train, mse_test), RMSE = c(rmse_train, rmse_test))
  row.names(mses) <- c("Training", "Validation")
  
  # Cross table and confusion matrix
  cross_table = table(predicted = as.logical(test_results$pred_binary), 
                      actual = as.logical(test_results$y))
  confusion <- confusionMatrix(cross_table, positive = "TRUE")
  results <- list(train_results, test_results)
  
  # ROC
  test_roc <- roc(test_results$y ~ test_results$pred, plot = TRUE, print.auc = TRUE)
  
  return(list(train_results, test_results, p1, p2, mses, confusion, test_roc))
}


################################################################################
# Preprocessing (Data Split)
################################################################################

# Set seed for reproducibility
set.seed(123)

# Step 1: Split data into training (70%) and temporary set (30%)
split <- sample.split(df$UTI_diag, SplitRatio = 0.7)
# 70% training data
train_all <- subset(df, split == TRUE)  
# 30% temporary data (will be split into validation and testing)
temp <- subset(df, split == FALSE) 

# Step 2: Split the temporary set into validation (15%) and testing (15%)
# Split 30% into 15% validation and 15% testing
split_temp <- sample.split(temp$UTI_diag, SplitRatio = 0.5)  
validation_all <- subset(temp, split_temp == TRUE)  # 15% validation data
test_all <- subset(temp, split_temp == FALSE)   # 15% testing data

# Step 3: Verify the sizes of each dataset
# nrow(train_all)         # Should be ~70% of the original data
# nrow(validation_all)    # Should be ~15% of the original data
# nrow(test_all)          # Should be ~15% of the original data

# Step 4: Backup datasets
train_all_just_in_case <- train_all
validation_all_just_in_case <- validation_all
test_all_just_in_case <- test_all

# Initialize manual data selection
train <- train_all
validation <- validation_all
test <- test_all

################################################################################
############################# Preprocessing ####################################
################################################################################

# NOTE: This is only for the models tested with manually selected data.
# XGBoost will use automatically selected variables.

# Extract numeric variables for EDA
tmp_num <- train[, sapply(train, is.numeric)]

# Check missing values in numeric variables
# apply(tmp_num, 2, function(x) sum(is.na(x)))  # Not many missing values

# Check correlation between numeric variables
# cor(na.omit(tmp_num))  # Little collinearity observed

# Identify and remove outliers in selected numeric variables
t1 <- quantile(tmp_num$ua_spec_grav, c(0.01, .99)) # 1st and 99th percentile
t2 <- quantile(tmp_num$age, c(0.01, .99)) # 1st and 99th percentile

# Remove extreme outliers (values below 1st percentile or above 99th percentile)
train <- train[train$ua_spec_grav > t1[1] & train$ua_spec_grav < t1[2], ]
train <- train[train$age > t2[1] & train$age < t2[2], ]

# Convert categorical variable to character
train$ua_ph <- as.character(train$ua_ph)
test$ua_ph <- as.character(test$ua_ph)

# Cleanup temporary variables
rm(t1, t2, split_temp, tmp_num)

# Selecting mean values from grouped measurements
measurements <- c('Temp_', 'HR_', 'SBP_', 'DBP_', 'RR_', 'O2_Sat_', 'O2_Amou')

rm_measurements <- function(x, m) {
  for (i in 1:length(m)) {
    # Remove the first 4 columns of each measurement type (keeping Mean)
    x <- x[, -grep(m[i], colnames(x))[1:4]]
  }
  # Remove redundant columns for 'O2_Dependency' and 'GCS_'
  x <- x[, -grep('O2_Depe', colnames(x))[2]]
  x <- x[, -grep('GCS_', colnames(x))[2]]
  return(x)
}

# Apply the function to train, validation and test datasets
train <- rm_measurements(train, measurements)
validation <- rm_measurements(validation, measurements)
test <- rm_measurements(test, measurements)

# Re-check missing values after dropping redundant columns
# head(sort(apply(train, 2, function(x) sum(is.na(x))), decreasing = TRUE), 10)

# Since all categorical variables are either on a Likert scale or binary (TRUE/FALSE),
# Replace NA values with 'not_reported' to maintain interpretability.

# Handle missing values by replacing NA with 'not_reported' for categorical variables
na_to_nr <- function(x) {
  x[, sapply(x, is.character)][is.na(x[, sapply(x, is.character)])] <- 'not_reported'
  return(x)
}

# Check missing values in test dataset before applying transformation
# head(sort(apply(test, 2, function(x) sum(is.na(x))), decreasing = TRUE), 10)

# Apply missing value replacement to train, validation and test datasets
train <- na_to_nr(train)
validation <- na_to_nr(validation)
test <- na_to_nr(test)

# Summarize dataset after preprocessing
# plot_intro(train)  # Overview of cleaned training dataset
# plot_intro(validation)  # Overview of cleaned validation dataset
# plot_intro(test)   # Overview of cleaned test dataset


################################################################################
# Feature Selection (Manually)
################################################################################

# Initial logistic regression model using all available predictors
glm0 <- glm(UTI_diag ~ ., family = binomial, data = train)
# summary(glm0)

# Selecting variables with p-value < 0.1 for feature refinement
vals_glm <- c('patid', 'ua_bacteria', 'ua_bili', 'ua_blood', 'ua_clarity', 
              'ua_color', 'ua_epi', 'ua_glucose', 'ua_ketones', 'ua_leuk', 
              'ua_nitrite', 'ua_ph', 'ua_rbc', 'ua_urobili', 'ua_wbc', 
              'CVA_tenderness', 'abd_mass', 'abd_rebound', 'vag_discharge', 
              'abd_distended2', 'gen_neg', 'pelvic_pain', 'weakness', 
              'psychiatric_confusion', 'flank_pain', 'diff_urinating', 
              'dysuria', 'hematuria', 'polyuria', 'chief_complaint', 'age', 
              'gender', 'ethnicity', 'employStatus', 'insurance_status', 
              'UTI_diag')

train_glm <- train[, vals_glm]

# Fit logistic regression on the refined set of predictors
glm1 <- glm(UTI_diag ~ ., family = binomial, data = train_glm)
# summary(glm1)

# Further variable selection: Removing variables with less significance
# Excluded: 'abd_rebound', 'ua_ketones', 'pelvic_pain', 'diff_urinating', 'polyuria'
vals_glm <- c('patid', 'ua_bacteria', 'ua_bili', 'ua_blood', 'ua_clarity', 
              'ua_color', 'ua_epi', 'ua_glucose', 'ua_leuk', 
              'ua_nitrite', 'ua_ph', 'ua_rbc', 'ua_urobili', 'ua_wbc', 
              'CVA_tenderness', 'abd_mass', 'vag_discharge', 
              'abd_distended2', 'gen_neg', 'weakness', 
              'psychiatric_confusion', 'flank_pain', 
              'dysuria', 'hematuria', 'chief_complaint', 'age', 
              'gender', 'ethnicity', 'employStatus', 'insurance_status', 
              'UTI_diag')

train_glm <- train[, vals_glm]
test_glm <- test[, vals_glm]
validation_glm <- validation[, vals_glm]

# Fit logistic regression after further feature refinement
glm2 <- glm(UTI_diag ~ ., family = binomial, data = train_glm)
# summary(glm2)

# Function to convert categorical variables into binary indicators 
# 'insurance_status', 'ethnicity'
convert_categorical <- function(x) {
  x$insurance_status_Medicare <- x[, 'insurance_status'] == 'Medicare'
  x$ethnicity_Non_Hispanic <- x[, 'ethnicity'] == 'Non-Hispanic'
  x$insurance_status <- NULL
  x$ethnicity <- NULL
  return(x)
}

# Apply categorical conversion
train_glm <- convert_categorical(train_glm)
test_glm <- convert_categorical(test_glm)
validation_glm <- convert_categorical(validation_glm)

# Fit logistic regression with transformed categorical variables
glm3 <- glm(UTI_diag ~ ., family = binomial, data = train_glm)
# summary(glm3)

# Perform backward selection to optimize model
glm_back <- step(glm3, test = "LRT", trace = 0)  

# Selected variables after backward selection
val_back <- c('patid' , 'ua_bacteria', 'ua_clarity', 'ua_color', 'ua_epi', 
              'ua_leuk', 'ua_nitrite', 'ua_urobili', 'ua_wbc', 
              'CVA_tenderness', 'abd_distended2', 'gen_neg', 
              'psychiatric_confusion', 'flank_pain', 'dysuria', 'hematuria', 
              'chief_complaint', 'age', 'gender', 'ethnicity_Non_Hispanic', 
              'UTI_diag')

# Compare model fit using AIC and BIC
# AIC(glm_back)  # 43571.79
# BIC(glm_back)  # 43829.49

# AIC(glm3)   # 40798.76
# BIC(glm3)   # 41838.44

# Keep only the selected variables in train and test sets
train <- train_glm[, val_back]
test <- test_glm[, val_back]
validation <- validation_glm[, val_back]

# Fit logistic regression on the final refined set of predictors
glm4 <- glm(UTI_diag ~ ., family = binomial, data = train)
# summary(glm4)

# Model comparison using ANOVA
# anova(glm3, glm2)  # Compare glm3 (with transformed categorical variables) to glm2
# anova(glm4, glm3)  # Compare glm4 (refined model) to glm3

# Model evaluation
# evaluate_glm(glm4, train, validation)  # AUC=0.8402
# evaluate_glm(glm4, train, test)  # AUC=0.8438

# **Variable Reduction: Simplified Model**
# Selecting a subset of the most important variables for a more compact model
# Initial selection of 15 variables based on p-values from glm4
val_step <- c('patid', 'ua_bacteria', 'ua_clarity', 'ua_color', 'ua_epi', 
              'ua_leuk', 'ua_nitrite', 'ua_urobili', 'ua_wbc', 'CVA_tenderness', 
              'psychiatric_confusion', 'flank_pain', 'age', 'gender', 
              'ethnicity_Non_Hispanic', 'UTI_diag')

train <- train[, val_step]
test <- test[, val_step]
validation <- validation[, val_step]

# Fit logistic regression with the simplified variable set
glm5 <- glm(UTI_diag ~ ., family = binomial, data = train)
# summary(glm5)

# **Final Model Selection: Keeping 11 Variables**
# The 'patid' variable is retained for now, though its relevance is uncertain.
val_step <- c('patid', 'ua_bacteria', 'ua_clarity', 'ua_epi', 
              'ua_wbc', 'CVA_tenderness', 'psychiatric_confusion', 
              'flank_pain', 'age', 'gender', 
              'ethnicity_Non_Hispanic', 'UTI_diag')

train <- train[, val_step]
test <- test[, val_step]
validation <- validation[, val_step]

# Fit logistic regression with 11 selected variables
glm6 <- glm(UTI_diag~., family=binomial, data=train)
# summary(glm6)

# Perform another backward selection to ensure model simplicity
glm_back <- step(glm6, test = 'LRT', trace = 0)  # Automatic feature selection

# **Final Decision:** 
# Retain CVA_tenderness for interpretability, despite model suggesting removal.

# **Evaluation of Final Model**
# Achieved AUC of **0.8117** with a streamlined model—strong performance.
# evaluate_glm(glm6, train, test)

# Top by GLM:
# 'patid', 'ua_bacteria', 'ua_clarity', 'ua_epi', 
# 'ua_wbc', 'CVA_tenderness', 'psychiatric_confusion', 
# 'flank_pain', 'age', 'gender', 
# 'ethnicity_Non_Hispanic', 'UTI_diag'


################################################################################
# Calculate Odds Ratio (OR) and 95% CI
################################################################################

# Specify the selected variables
val_final <- c('patid', 'ua_bacteria', 'ua_clarity', 'ua_epi', 
               'ua_wbc', 'CVA_tenderness', 'psychiatric_confusion', 
               'flank_pain', 'age', 'gender', 
               'ethnicity_Non_Hispanic', 'UTI_diag')

train <- train[, val_final]
test <- test[, val_final]
validation <- validation[, val_final]

# Fit logistic regression with the selected variables
glm_final <- glm(UTI_diag~., family=binomial, data=train)
# summary(glm_final)

# Calculate odds ratios and 95% confidence intervals
fitOR <- exp(cbind(OR = coef(glm_final), confint(glm_final)))

# Save results
saveRDS(fitOR, 'OR.RDS')

# Print top 10 variables based on OR magnitude
fitOR_sorted <- fitOR[order(-fitOR[, 1]), ]  # Sort by OR in descending order
# print(head(fitOR_sorted, 10))  # Display top 10 variables

# Convert fitOR to a data frame to ensure compatibility with kable
fitOR <- data.frame(fitOR)

# Generate a LaTeX-formatted table for all variables
kable(
  fitOR, 
  booktabs = TRUE, align = "c",
  col.names = c("Variable", "Odds Ratio (OR)", "Lower 95% CI", "Upper 95% CI"),
  caption = "Table 1: Odds Ratios with 95% Confidence Intervals"
) %>%
  kable_styling(
    latex_options = c("striped", "scale_down", "HOLD_position"), 
    full_width = TRUE
  ) %>%
  column_spec(1:4, latex_column_spec = "c")  # Center all columns

# Select top 10 variables based on odds ratio magnitude
fitOR_sorted <- data.frame(head(fitOR_sorted, 10))

# Generate a LaTeX-formatted table for the top 10 variables
kable(
  fitOR_sorted, 
  booktabs = TRUE, align = "c",
  col.names = c("Variable", "Odds Ratio (OR)", "Lower 95% CI", "Upper 95% CI"),
  caption = "Table 2: Top 10 Variables by Odds Ratio"
) %>%
  kable_styling(
    latex_options = c("striped", "scale_down", "HOLD_position"), 
    full_width = TRUE
  ) %>%
  column_spec(1:4, latex_column_spec = "c")  # Center all columns


################################################################################
# Predict Probabilities for Validation and Test Data
################################################################################

train$predicted_prob <- predict(glm6, train, type = "response")
validation$predicted_prob <- predict(glm6, validation, type = "response")
test$predicted_prob <- predict(glm6, test, type = "response")

################################################################################
# Evaluate AUC-ROC Curve for Validation and Testing Data
################################################################################

# Compute ROC Curves
roc_valid <- roc(validation$UTI_diag, validation$predicted_prob)
roc_test <- roc(test$UTI_diag, test$predicted_prob)

# Plot ROC Curves
par(mfrow = c(1, 2))  # Side-by-side plots
plot(roc_valid, col = "blue", lwd = 2, main = "AUC-ROC Curve (Validation)")
legend("bottomright", legend = paste("AUC =", round(auc(roc_valid), 4)), 
       col = "blue", lwd = 2)

plot(roc_test, col = "red", lwd = 2, main = "AUC-ROC Curve (Testing)")
legend("bottomright", legend = paste("AUC =", round(auc(roc_test), 4)), 
       col = "red", lwd = 2)

dev.copy(png, "AUC_logistic.png", width = 1000, height = 500)
dev.off()

# Print AUC Values
# print(paste("Validation AUC-ROC:", round(auc(roc_valid), 4)))  # 0.8106
# print(paste("Testing AUC-ROC:", round(auc(roc_test), 4)))  # 0.8117

# The AUC-ROC values for both datasets are very similar and indicate that 
# the model is good at distinguishing between the classes, with values above 
# 0.80 suggesting excellent discriminative ability.


################################################################################
# Evaluate Precision-Recall AUC for Validation and Testing Data
################################################################################

# Convert target variable to numeric
y_valid_numeric <- as.numeric(validation$UTI_diag)
y_test_numeric <- as.numeric(test$UTI_diag)

# Compute PR Curves
pr_valid <- pr.curve(scores.class0 = validation$predicted_prob[y_valid_numeric == 0],  
                     scores.class1 = validation$predicted_prob[y_valid_numeric == 1],  
                     curve = TRUE)

pr_test <- pr.curve(scores.class0 = test$predicted_prob[y_test_numeric == 0],  
                    scores.class1 = test$predicted_prob[y_test_numeric == 1],  
                    curve = TRUE)

# Plot PR Curves
par(mfrow = c(1, 2))
plot(pr_valid, col = "blue", main = "Precision-Recall Curve (Validation)", lwd = 2)
legend("bottomright", legend = paste("AUC-PR =", round(pr_valid$auc.integral, 4)), 
       col = "blue", lwd = 2)

plot(pr_test, col = "red", main = "Precision-Recall Curve (Testing)", lwd = 2)
legend("bottomright", legend = paste("AUC-PR =", round(pr_test$auc.integral, 4)), 
       col = "red", lwd = 2)

dev.copy(png, "PR_AUC_logistic.png", width = 1000, height = 500)
dev.off()

# Print PR AUC Values
# print(paste("Validation AUC-PR:", round(pr_valid$auc.integral, 4)))  # 0.6349
# print(paste("Testing AUC-PR:", round(pr_test$auc.integral, 4)))  # 0.6352

# The PR AUC values are relatively high, especially considering the class 
# imbalance often present in healthcare datasets. This suggests that the model 
# is performing well, even for the less frequent class, by correctly identifying 
# positive instances.


################################################################################
# Evaluate Calibration Curve for Validation and Testing Data
################################################################################

# Function to compute and plot calibration curve
plot_calibration_curve <- function(pred_probs, true_labels, dataset_name) {
  # Compute calibration curve manually
  bin_cutoffs <- quantile(pred_probs, probs = seq(0, 1, length.out = 11))  # 10 bins
  bin_labels <- cut(pred_probs, breaks = bin_cutoffs, include.lowest = TRUE, labels = FALSE)
  
  # Aggregate observed proportions per bin
  cal_data <- data.frame(
    bin = bin_labels,
    predicted = pred_probs,
    observed = true_labels
  )
  
  cal_summary <- aggregate(observed ~ bin, data = cal_data, FUN = mean)
  cal_summary$predicted <- aggregate(predicted ~ bin, data = cal_data, FUN = mean)$predicted
  
  # Plot Calibration Curve
  ggplot(cal_summary, aes(x = predicted, y = observed)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = paste("Calibration Curve (", dataset_name, ")", sep = ""),
         x = "Predicted Probability",
         y = "Observed Proportion") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line()
    )
}

# Plot Calibration Curves for Validation and Testing
p1 <- plot_calibration_curve(validation$predicted_prob, y_valid_numeric, "Validation")
p2 <- plot_calibration_curve(test$predicted_prob, y_test_numeric, "Testing")
grid.arrange(p1, p2, ncol = 2)

dev.copy(png, "Calibration_logistic_1.png", width = 1000, height = 500)
dev.off()

# Plot Calibration Curves (val.prob) for Validation and Testing
par(mfrow = c(1,2))
cal_curve_valid <- val.prob(validation$predicted_prob, y_valid_numeric, pl = TRUE)
title(main = "Calibration Curve for Validation (val.prob)")

cal_curve_test <- val.prob(test$predicted_prob, y_test_numeric, pl = TRUE)
title(main = "Calibration Curve for Testing (val.prob)")
par(mfrow = c(1,1))

dev.copy(png, "Calibration_logistic_2.png", width = 1000, height = 500)
dev.off()


################################################################################
# XGBoost Implementation
################################################################################

# Preprocessing: Clean and prepare data
train_all <- train_all_just_in_case
validation_all <- validation_all_just_in_case
test_all <- test_all_just_in_case

# Handle missing values
train_all <- na_to_nr(train_all)
validation_all <- validation_all_just_in_case
test_all <- na_to_nr(test_all)

# Replace NA values in categorical variables
train_all$ua_ph[is.na(train_all$ua_ph)] <- 'not_reported'
validation_all$ua_ph[is.na(validation_all$ua_ph)] <- 'not_reported'
test_all$ua_ph[is.na(test_all$ua_ph)] <- 'not_reported'

# Convert data to dataframes
train_all <- as.data.frame(train_all)
validation_all <- as.data.frame(validation_all)
test_all <- as.data.frame(test_all)

# Prepare input features and target variables
X_train <- data.matrix(train_all[,-ncol(train_all)])
y_train <- train_all[,ncol(train_all)]

X_validation <- data.matrix(validation_all[,-ncol(validation_all)])
y_validation <- validation_all[,ncol(validation_all)]

X_test <- data.matrix(test_all[,-ncol(test_all)])   
y_test <- test_all[,ncol(test_all)]

# Convert the train, validation and test data into XGBoost matrix type
xgboost_train <- xgb.DMatrix(data=X_train, label=y_train)
xgboost_validation <- xgb.DMatrix(data=X_validation, label=y_validation)
xgboost_test <- xgb.DMatrix(data=X_test, label=y_test)

# Train XGBoost Model
model <- xgboost(data = xgboost_train,           # The data   
                 max.depth=3,                    # Max depth 
                 nrounds=50,
                 objective = "binary:logistic",  # Binary classification with logistic loss
                 verbose = 0)                    # Suppress training output
# summary(model)

# Model Evaluation: Confusion Matrix and Accuracy
# Predict on validation and test sets
pred_valid <- predict(model, xgboost_validation)
pred_test <- predict(model, xgboost_test)

# Convert probabilities to binary labels
pred_valid_binary <- as.factor(ifelse(pred_valid > 0.5, TRUE, FALSE))
pred_test_binary <- as.factor(ifelse(pred_test > 0.5, TRUE, FALSE))

# Create confusion matrices
conf_mat_valid <- confusionMatrix(as.factor(y_validation), pred_valid_binary)
conf_mat_test <- confusionMatrix(as.factor(y_test), pred_test_binary)

# Print confusion matrix results
# print(conf_mat_valid)
# print(conf_mat_test)

# Accuracy on Validation and Test
validation_accuracy <- sum(diag(conf_mat_valid$table)) / sum(conf_mat_valid$table)
test_accuracy <- sum(diag(conf_mat_test$table)) / sum(conf_mat_test$table)

# print(paste("Validation Accuracy: ", round(validation_accuracy, 4)))  # 0.889
# print(paste("Test Accuracy: ", round(test_accuracy, 4)))  # 0.896


################################################################################
# XGBoost Feature Importance Analysis
################################################################################

# Compute feature importance matrix
importance_matrix <- xgb.importance(colnames(xgboost_train), model = model)

# Print importance matrix
# print(importance_matrix)

# Adjust y-axis label position
par(mgp = c(3, 0, 0))

# Plot top 20 important features
xgb.plot.importance(importance_matrix[1:20,], 
                    main = "Top 20 Feature Importance - XGBoost",
                    xlab = "Importance Score",
                    ylab = "Features")

# Reset par to default after plotting
par(mgp = c(3, 1, 0))

dev.copy(png, "Feature_Importance_XGBoost.png", width = 1000, height = 500)
dev.off()

# Prepare the data for explanation
# topvars <- importance_matrix[1:20,]$Feature

# Plot the multi-tree structure of the model
# xgb.plot.multi.trees(model = model, use.names = FALSE, fill = TRUE)

# Top 10 XGBoost
# 'abx', 'patid', 'ua_wbc', 'ua_leuk', 'ua_nitrite', 
# 'ua_bacteria', 'dispo', 'antibiotics', 'age', 'chief_complaint'


################################################################################
# SHAP Analysis for XGBoost Model
################################################################################

# 1. SHAP Summary Plot (Top 20 Features)
# Compute SHAP values for the training dataset
shap_values <- shap.prep(xgb_model = model, X_train = X_train)

# Ensure SHAP values are a data.table
shap_long <- shap_values  # `shap.prep()` already returns a data.table

# Select top 20 features by mean absolute SHAP value
top_20_features <- shap_long[, .(mean_shap = mean(abs(value))), 
                             by = variable][order(-mean_shap)][1:20, variable]

# Filter dataset to include only the top 20 features
shap_top_20 <- shap_long[variable %in% top_20_features]

# SHAP Summary Plot (Top 20 Features)
ggplot(shap_top_20, aes(x = value, y = variable, color = value)) +
  geom_jitter(size = 0.8, alpha = 0.4) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  labs(title = "SHAP Summary Plot - Top 20 Features", 
       x = "SHAP Value (impact on model output)", y = "Feature") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

dev.copy(png, "SHAP Summary (Top 20 Features).png", width = 1000, height = 500)
dev.off()


################################################################################
# SHAP Feature Importance Plot - Top 20 Features
################################################################################

# 2. SHAP Feature Importance Plot (Top 20 Features)
# Compute feature importance as mean absolute SHAP values
shap_importance_df <- shap_long[, .(Importance = mean(abs(value))), by = variable]

# Select top 20 features by importance
top_20_importance <- shap_importance_df[order(-Importance)][1:20]

# SHAP Feature Importance Plot (Top 20 Features)
ggplot(top_20_importance, aes(x = reorder(variable, Importance), y = Importance)) +
  geom_col(fill = viridis::viridis(20, option = "cividis"), color = "black", width = 0.7) +
  coord_flip() +
  labs(title = "SHAP Feature Importance - Top 20 Features", 
       x = "Feature", y = "Mean |SHAP Value| (average impact on model output magnitude)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

dev.copy(png, "SHAP Feature Importance (Top 20 Features).png", width = 1000, height = 500)
dev.off()


################################################################################
# SHAP Dependence Plot for White Blood Cell Count (`ua_wbc`)
################################################################################

# 3. SHAP Dependence Plot for White Blood Cell Count (`ua_wbc`)
# Ensure 'ua_wbc' is in the dataset
if ("ua_wbc" %in% shap_long$variable) {
  
  # Extract SHAP values and corresponding feature values
  shap_wbc <- shap_long[variable == "ua_wbc", .(rfvalue, shap_value = value)]
  
  # Plot SHAP dependence
  ggplot(shap_wbc, aes(x = rfvalue, y = shap_value)) +
    geom_point(alpha = 0.5, size = 1, color = "blue") +  # Uniform point color
    geom_smooth(method = "loess", color = "red", linewidth = 1, se = FALSE) +
    scale_x_continuous(limits = c(1, 6), breaks = 1:6) +  # Set x-axis range
    labs(title = "SHAP Dependence Plot: White Blood Cell Count (ua_wbc)", 
         x = "White Blood Cell Count (ua_wbc)", 
         y = "SHAP Value for ua_wbc") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.position = "none",  # Remove legend since color gradient is removed
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
} else {
  print("Warning: 'ua_wbc' not found in SHAP dataset!")
}

dev.copy(png, "SHAP Dependence Plot (ua_wbc).png", width = 1000, height = 500)
dev.off()


################################################################################
# AUC-ROC Curve (XGBoost)
################################################################################

# Compute ROC Curves
roc_valid <- roc(y_validation, pred_valid)
roc_test <- roc(y_test, pred_test)

# Plot ROC Curves
par(mfrow = c(1, 2))  # Side-by-side plots
plot(roc_valid, col = "blue", lwd = 2, main = "AUC-ROC Curve (Validation)")
legend("bottomright", legend = paste("AUC =", round(auc(roc_valid), 4)), 
       col = "blue", lwd = 2)

plot(roc_test, col = "red", lwd = 2, main = "AUC-ROC Curve (Testing)")
legend("bottomright", legend = paste("AUC =", round(auc(roc_test), 4)), 
       col = "red", lwd = 2)

dev.copy(png, "AUC_XGBoost.png", width = 1000, height = 500)
dev.off()

# Print AUC Values
# print(paste("Validation AUC-ROC:", round(auc(roc_valid), 4)))  # 0.9348
# print(paste("Test AUC-ROC:", round(auc(roc_test), 4)))  # 0.9564

# The AUC-ROC values for both validation and test sets are very high (>0.93), 
# indicating that the model has excellent discriminatory power between the positive and
# negative classes. The slight improvement in the test set (0.9564) suggests that 
# the model generalizes well and maintains strong predictive performance on unseen data.


################################################################################
# Precision-Recall AUC (XGBoost)
################################################################################

# Compute PR AUC for validation and test
pr_valid <- pr.curve(scores.class0 = pred_valid[y_validation == 0],  
                     scores.class1 = pred_valid[y_validation == 1],  
                     curve = TRUE)

pr_test <- pr.curve(scores.class0 = pred_test[y_test == 0],  
                    scores.class1 = pred_test[y_test == 1],  
                    curve = TRUE)

# Plot PR Curves
par(mfrow = c(1, 2))
plot(pr_valid, col = "blue", main = "Precision-Recall Curve (Validation)", lwd = 2)
legend("bottomright", legend = paste("AUC-PR =", round(pr_valid$auc.integral, 4)), 
       col = "blue", lwd = 2)

plot(pr_test, col = "red", main = "Precision-Recall Curve (Testing)", lwd = 2)
legend("bottomright", legend = paste("AUC-PR =", round(pr_test$auc.integral, 4)), 
       col = "red", lwd = 2)

dev.copy(png, "PR_AUC_XGBoost.png", width = 1000, height = 500)
dev.off()

# Print PR AUC Values
# print(paste("Validation AUC-PR:", round(pr_valid$auc.integral, 4)))  # 0.5934
# print(paste("Test AUC-PR:", round(pr_test$auc.integral, 4)))  # 0.5878

# The AUC-PR values are relatively moderate (~0.59), which suggests that while the model 
# performs well in terms of overall classification (as indicated by AUC-ROC), there may be 
# some room for improvement in handling class imbalance. AUC-PR is particularly useful when 
# dealing with imbalanced datasets, as it focuses on the precision-recall tradeoff rather than
# the overall classification.


################################################################################
# Calibration Curve (XGBoost)
################################################################################

# Function to compute and plot calibration curve
plot_calibration_curve <- function(pred_probs, true_labels, dataset_name) {
  # Compute calibration curve manually
  bin_cutoffs <- quantile(pred_probs, probs = seq(0, 1, length.out = 11))  # 10 bins
  
  # Aggregate observed proportions per bin
  bin_labels <- cut(pred_probs, breaks = bin_cutoffs, include.lowest = TRUE, labels = FALSE)
  cal_data <- data.frame(
    bin = bin_labels,
    predicted = pred_probs,
    observed = true_labels
  )
  
  cal_summary <- aggregate(observed ~ bin, data = cal_data, FUN = mean)
  cal_summary$predicted <- aggregate(predicted ~ bin, data = cal_data, FUN = mean)$predicted
  
  # Plot Calibration Curve
  ggplot(cal_summary, aes(x = predicted, y = observed)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = paste("Calibration Curve (", dataset_name, ")", sep = ""),
         x = "Predicted Probability",
         y = "Observed Proportion") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line()
    )
}

# Plot Calibration Curves for Validation and Testing (XGBoost)
p5 <- plot_calibration_curve(pred_valid, y_validation, "Validation - XGBoost")
p6 <- plot_calibration_curve(pred_test, y_test, "Testing - XGBoost")
grid.arrange(p5, p6, ncol = 2)

dev.copy(png, "Calibration_XGBoost_1.png", width = 1000, height = 500)
dev.off()

# Plot Calibration Curves (val.prob) for Validation and Testing (XGBoost)
par(mfrow = c(1,2))
cal_curve_valid <- val.prob(pred_valid, y_validation, pl = TRUE)
title(main = "Calibration Curve for Validation (val.prob)")

cal_curve_test <- val.prob(pred_test, y_test, pl = TRUE)
title(main = "Calibration Curve for Testing (val.prob)")
par(mfrow = c(1,1))

dev.copy(png, "Calibration_XGBoost_2.png", width = 1000, height = 500)
dev.off()



################################################################################
# Log training parameters, metrics, and figures to MLFlow
################################################################################

# Load necessary libraries
library(reticulate)   # For interface between R and Python  
library(mlflow)       # For track experiments and manage ML models

# Activate the Conda environment
reticulate::use_condaenv("bis568ops", required=TRUE)

# Set the correct path to mlflow executable
Sys.setenv(PATH = paste("E:/miniconda3/envs/bis568ops/Scripts", 
                        Sys.getenv("PATH"), sep = ";"))

# Verify the path has been set correctly
Sys.getenv("PATH")

# Launch mlflow UI
# mlflow::mlflow_set_tracking_uri("http://127.0.0.1:5000")
mlflow::mlflow_ui()

################################################################################
# Logistic Regression Model
################################################################################

mlflow::mlflow_start_run()

# Log parameters (selected features)
selected_features <- paste(names(train)[-ncol(train)], collapse=", ")
mlflow::mlflow_log_param("model_type", "logistic_regression")
mlflow::mlflow_log_param("selected_features", selected_features)
# mlflow::mlflow_log_param("selected_features", "ua_bacteria, ua_clarity, ua_epi, 
#                          ua_wbc, CVA_tenderness, psychiatric_confusion, flank_pain, 
#                          age, gender, ethnicity_Non_Hispanic")

# Log metrics
mlflow::mlflow_log_metric("roc_valid", 0.8106)
mlflow::mlflow_log_metric("roc_test", 0.8117)
mlflow::mlflow_log_metric("pr_valid", 0.6349)
mlflow::mlflow_log_metric("pr_test", 0.6352)

# Log artifacts (figures)
mlflow::mlflow_log_artifact("AUC_logistic.png")
mlflow::mlflow_log_artifact("PR_AUC_logistic.png")
mlflow::mlflow_log_artifact("Calibration_logistic_1.png")
mlflow::mlflow_log_artifact("Calibration_logistic_2.png")

# Properly log the GLM model
mlflow::mlflow_log_model(
  model = carrier::crate(
    function(new_data) {
      predict(!!glm6, newdata = new_data, type = "response")
    },
    glm6 = glm6
  ),
  artifact_path = "logistic_model"
)

mlflow::mlflow_end_run()

################################################################################
# XGBoost
################################################################################

mlflow::mlflow_start_run()

# Log model parameters
mlflow::mlflow_log_param("model_type", "xgboost")
mlflow::mlflow_log_param("data", "xgboost_train")
mlflow::mlflow_log_param("max_depth", 3)
mlflow::mlflow_log_param("nrounds", model$niter)  # Log actual iterations (50)
mlflow::mlflow_log_param("objective", "binary:logistic")
mlflow::mlflow_log_param("missing_value_handling", "not_reported")

# Log performance metrics
# validation_accuracy <- sum(diag(conf_mat_valid$table)) / sum(conf_mat_valid$table)
# test_accuracy <- sum(diag(conf_mat_test$table)) / sum(conf_mat_test$table)
# mlflow::mlflow_log_metric("validation_accuracy", validation_accuracy)
# mlflow::mlflow_log_metric("test_accuracy", test_accuracy)
# mlflow::mlflow_log_metric("roc_valid", pROC::roc(y_validation, pred_valid)$auc)
# mlflow::mlflow_log_metric("roc_test", pROC::roc(y_test, pred_test)$auc)
mlflow::mlflow_log_metric("validation_accuracy", 0.889)
mlflow::mlflow_log_metric("test_accuracy", 0.896)
mlflow::mlflow_log_metric("roc_valid", 0.9348)
mlflow::mlflow_log_metric("roc_test", 0.9564)
mlflow::mlflow_log_metric("pr_valid", 0.5934)
mlflow::mlflow_log_metric("pr_test", 0.5878)

# Log confusion matrix metrics
mlflow::mlflow_log_metric("validation_sensitivity", conf_mat_valid$byClass["Sensitivity"])
mlflow::mlflow_log_metric("validation_specificity", conf_mat_valid$byClass["Specificity"])
mlflow::mlflow_log_metric("test_sensitivity", conf_mat_test$byClass["Sensitivity"])
mlflow::mlflow_log_metric("test_specificity", conf_mat_test$byClass["Specificity"])

# Log top features
top_features <- importance_matrix[1:10, Feature]
mlflow::mlflow_log_param("top_10_features", paste(top_features, collapse = ", "))

# Log all visualization artifacts
mlflow::mlflow_log_artifact("Feature_Importance_XGBoost.png")
mlflow::mlflow_log_artifact("SHAP Summary (Top 20 Features).png")
mlflow::mlflow_log_artifact("SHAP Feature Importance (Top 20 Features).png")
mlflow::mlflow_log_artifact("SHAP Dependence Plot (ua_wbc).png")
mlflow::mlflow_log_artifact("AUC_XGBoost.png")
mlflow::mlflow_log_artifact("PR_AUC_XGBoost.png")
mlflow::mlflow_log_artifact("Calibration_XGBoost_1.png")
mlflow::mlflow_log_artifact("Calibration_XGBoost_2.png")

# Properly log the XGBoost model
mlflow::mlflow_log_model(
  model = carrier::crate(
    function(new_data) {
      predict(!!model, as.matrix(new_data))
    },
    model = model
  ),
  artifact_path = "xgboost_model"
)

mlflow::mlflow_end_run()

