
# Load necessary libraries
library(dplyr)             # For data manipulation
library(tidyverse)         # For data science workflows
library(ggplot2)           # For creating visualizations and plots
library(caret)             # For machine learning modeling and preprocessing
library(randomForestSRC)   # For survival analysis and random forests
library(randomForest)      # For random forest models
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
raw_data <- tbl(db_conn, 'results') %>% collect()

# Drop irrelevant columns
df <- raw_data %>% dplyr::select(-UCX_abnormal, -split, -alt_diag, -abxUTI)

# Inspect dataset
str(df)
summary(df)
head(df)
if (ncol(df) >= 216) {
  head(df[, 210:216])
}

# Observations:
# - NA values present
# - Inconsistent data types
# - Logical values stored as integers

# Check for redundant variables, i.e. only have NA values
which(apply(df, 2, function(x) sum(is.na(x))) == nrow(df))

# Remove redundant variables
df <- df[, -which(apply(df, 2, function(x) sum(is.na(x))) == nrow(df))]

# Check for only one response
which(apply(df, 2, function(x) length(unique(x))) == 1)

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

# Evaluate random forest models
evaluate_rf <- function(model, train, test) {
  pred_train <- predict(model, type = "prob")
  pred_test <- predict(model, newdata = test, type = "prob")
  
  train_results <- data.frame(y = train$UTI_diag, pred = pred_train[,2], 
                              pred_binary = as.factor(ifelse(pred_train[,1] < 0.5, TRUE, FALSE)))
  test_results <- data.frame(y = test$UTI_diag, pred = pred_test[,2], 
                             pred_binary = as.factor(ifelse(pred_test[,1] < 0.5, TRUE, FALSE)))
  
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
# Exploratory Data Analysis (EDA)
################################################################################

# Summary Statistics
# summary(df)
# Checking structure of dataset
# str(df)

# ==================== Looking at the continuous variables vs response variable
# Continuous Variables - Violin and Box Plots
# Urine spec_grav
p1 <- ggplot(df, aes(x= UTI_diag, y=ua_spec_grav, fill = UTI_diag)) + 
  geom_violin() +
  labs(title="Urine spec_grav vs. response variable",
       subtitle = 'Response variable: Was UTI diagnosed?',
       x="UTI diagnosed?", y = "spec_grav") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# No spread of data, but some extreme outliers in FALSE
# The very long whisker extending upwards suggests that 
# some values are far outside the main distribution.

# pH value in urine
p2 <- ggplot(df, aes(x=as.factor(ua_ph), fill = UTI_diag)) +
  geom_bar() + 
  labs(title="pH value in urine vs. response variable",
       subtitle = 'Response variable: Was UTI diagnosed?') + 
  xlab("Ph Value") + ylab("Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Age
p3 <- ggplot(df, aes(x= UTI_diag, y=age, fill = UTI_diag)) + 
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title="Age vs. response variable",
       subtitle = 'Response variable: Was UTI diagnosed?',
       x="UTI diagnosed?", y = "Age") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# ==================== Looking at the categorical variables vs response variable
# Categorical Variables - Bar Charts
# Gender
p4 <- ggplot(df, aes(x=gender, fill = UTI_diag)) +
  geom_bar() +
  labs(title = "Gender Distribution", x="Gender", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Race
p5 <- ggplot(df, aes(x=race, fill = UTI_diag)) +
  geom_bar() +
  labs(title = "Race Distribution", x="Race", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Ethnicity
p6 <- ggplot(df, aes(x=ethnicity, fill = UTI_diag)) +
  geom_bar() +
  labs(title = "Ethnicity Distribution", x="Ethnicity", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Urine Clarity
p7 <- ggplot(df, aes(x=ua_clarity, fill = UTI_diag)) +
  geom_bar() +
  labs(title = "Urine Clarity Distribution", x="Clarity", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Urine Color
p8 <- ggplot(df, aes(x=ua_color, fill = UTI_diag)) +
  geom_bar() +
  labs(title = "Urine Color Distribution", x="Color", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Urine Bacteria
p9 <- ggplot(df, aes(x=ua_bacteria, fill = UTI_diag)) +
  geom_bar() +
  labs(title = "Urine Bacteria Distribution", x="Bacteria", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Load library
library(gridExtra)

# Arrange all plots in a single panel
grid.arrange(p1, p2, p3, ncol = 2)
grid.arrange(p4, p6, p7, ncol = 2)
grid.arrange(p8, p9, p5, 
             layout_matrix = rbind(c(1, 2), 
                                   c(3, 3)))


# ==================== Visualizing Descriptive Statistics ====================
# 1. Patient Age at Visit (in years)
# Age Distribution (histogram)
p1 <- ggplot(df, aes(x=age)) + 
  geom_histogram(bins = 30) + 
  labs(title="Age Distribution", x="Age at Visit (Years)", y="Count") + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Age Distribution by Group (histogram)
p1.1 <- ggplot(df, aes(x=age, fill=factor(gender))) + 
  geom_histogram(bins = 30, position="stack") +
  labs(title="Age Distribution by Gender", x="Age at Visit (Years)", y="Count", 
       fill="Gender") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Create a frequency table for age
age_count <- table(df$age)
# Find the maximum count
max_count <- max(age_count)
# Find the age(s) with the maximum count
age_with_max_count <- as.numeric(names(age_count)[age_count == max_count])
# print(age_with_max_count)

# 2. Gender Distribution (Bar Plot)
p2 <- ggplot(df, aes(x=gender, fill=factor(gender))) +
  geom_bar() +
  labs(title="Gender Distribution", x="Gender", y="Count", fill="Gender") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# 3. Clarity Distribution (Bar Plot)
p3 <- ggplot(df, aes(x=ua_clarity, fill=ua_clarity)) +
  geom_bar() +
  labs(title="Urine Clarity Distribution", x="Clarity", y="Count", fill="Clarity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Combine all bar plots/histograms into a single list
part1_plots <- list(p1, p1.1, p2, p3)

# Display the plots
grid.arrange(grobs = part1_plots, ncol = 2)


# ==================== Visualizing Descriptive Statistics ====================
# 4. Bacteria Distribution (Bar Plot)
p4 <- ggplot(df, aes(x=ua_bacteria, fill=ua_bacteria)) +
  geom_bar() +
  labs(title="Bacteria Distribution", x="Bacteria", y="Count", fill="Bacteria") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# 5. WBC Distribution (Bar Plot)
p5 <- ggplot(df, aes(x=ua_wbc, fill=ua_wbc)) +
  geom_bar() +
  labs(title="White Blood Cell (WBC) Distribution", x="WBC Count", y="Count", 
       fill="WBC") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Combine all bar plots/histograms into a single list
part2_plots <- list(p4, p5)

# Display the plots
grid.arrange(grobs = part2_plots, ncol = 2)

# ==================== Visualizing Descriptive Statistics ====================
# 6. Race Distribution (Bar Plot)
p6 <- ggplot(df, aes(x=race, fill=race)) +
  geom_bar() +
  labs(title="Race Distribution", x="Race", y="Count", fill="Race") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# 7. Ethnicity Distribution (Bar Plot)
p7 <- ggplot(df, aes(x=ethnicity, fill=ethnicity)) +
  geom_bar() +
  labs(title="Ethnicity Distribution", x="Ethnicity", y="Count", fill="Ethnicity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Combine all bar plots/histograms into a single list
part3_plots <- list(p6, p7)

# Display the plots
grid.arrange(grobs = part3_plots, nrow = 2)


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
nrow(train_all)     # Should be ~70% of the original data
nrow(validation_all)    # Should be ~15% of the original data
nrow(test_all)      # Should be ~15% of the original data

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
apply(tmp_num, 2, function(x) sum(is.na(x)))  # Not many missing values

# Check correlation between numeric variables
cor(na.omit(tmp_num)) # Little collinearity observed

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

################################################################################
############################# Categorical Data #################################
################################################################################

# Extract categorical variables from the training dataset
tmp_char <- train[, sapply(train, is.character)]

# Visualize missing values
missing_plot <- apply(tmp_char, 2, function(x) sum(is.na(x))) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  ggplot(aes(x = rowname, y = .)) +
  geom_bar(stat = 'identity', fill = 'tomato3') +
  labs(x = 'Variable', y = 'Number of Missing Values', 
       title = 'Number of Missing Values') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

print(missing_plot)

# Calculate percentage of missing values for each categorical variable
missing_values <- tmp_char %>%
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  summarise(total = n(),
            num.isna = sum(isna),
            pct = num.isna / total * 100) %>%
  arrange(desc(pct))

# Order categorical variables by percentage of missing values
levels <- missing_values$key

# Plot missing values percentage using a dot plot
percentage_plot <- missing_values %>%
  ggplot(aes(x = reorder(key, pct), y = pct)) +
  geom_point(color = "tomato3", size = 2) +
  geom_segment(aes(x = key, xend = key, y = 0, yend = pct), color = "steelblue") +
  scale_x_discrete(limits = levels) +
  coord_flip() +
  labs(title = "Percentage of Missing Values", x = "Variable", y = "% Missing") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

print(percentage_plot)

# Identify top 10 variables with the most missing values
head(sort(apply(tmp_char, 2, function(x) sum(is.na(x))), decreasing = TRUE), 10)

# Pause handling missing values for now, as they may contain useful information.
# Next, examine clusters of measurement variables that include First, Last, Max, Min, Mean.

# Measurement groups to analyze:
# 1. Temp_ (Temperature)
# 2. HR_ (Heart Rate)
# 3. SBP_ (Systolic Blood Pressure)
# 4. DBP_ (Diastolic Blood Pressure)
# 5. RR_ (Respiratory Rate)
# 6. O2_Sat_ (Oxygen Saturation)
# 7. O2_Amount_ (Oxygen Amount)
# 8. O2_Dependency_ (Check later, treating NA as 'not_reported')
# 9. GCS_ (Glasgow Coma Scale, check later, treating NA as 'not_reported')

# Extract measurements from the training dataset
# c1 <- train[, grep('Temp_', colnames(train))]   # Keep Temp_Mean [5]
# c2 <- train[, grep('HR_', colnames(train))]     # Keep HR_Mean [5]
# c3 <- train[, grep('SBP_', colnames(train))]    # Keep SBP_Mean [5]
# c4 <- train[, grep('DBP_', colnames(train))]    # Keep DBP_Mean [5]
# c5 <- train[, grep('RR_', colnames(train))]     # Keep RR_Mean [5]
# c6 <- train[, grep('O2_Sat_', colnames(train))] # Keep O2_Sat_Mean [5]
# c7 <- train[, grep('O2_Amou', colnames(train))] # Keep O2_Amount_Mean [5]

# Since most of these measurements are highly correlated, drop redundant columns.

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
head(sort(apply(train, 2, function(x) sum(is.na(x))), decreasing = TRUE), 10)

# Since all categorical variables are either on a Likert scale or binary (TRUE/FALSE),
# Replace NA values with 'not_reported' to maintain interpretability.

# Handle missing values by replacing NA with 'not_reported' for categorical variables
na_to_nr <- function(x) {
  x[, sapply(x, is.character)][is.na(x[, sapply(x, is.character)])] <- 'not_reported'
  return(x)
}

# Check missing values in test dataset before applying transformation
head(sort(apply(test, 2, function(x) sum(is.na(x))), decreasing = TRUE), 10)

# Apply missing value replacement to train, validation and test datasets
train <- na_to_nr(train)
validation <- na_to_nr(validation)
test <- na_to_nr(test)

# Summarize dataset after preprocessing
plot_intro(train)  # Overview of cleaned training dataset
plot_intro(validation)  # Overview of cleaned validation dataset
plot_intro(test)   # Overview of cleaned test dataset

# Write out data used to train LR/tree-based models
write.csv(train, 'dat.csv')


################################################################################
# Feature Selection (Manually)
################################################################################

# Initial logistic regression model using all available predictors
glm0 <- glm(UTI_diag ~ ., family = binomial, data = train)
summary(glm0)

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
summary(glm1)

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
summary(glm2)

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
summary(glm3)

# Perform backward selection to optimize model
glm_back <- step(glm3, test = "LRT")  

# Selected variables after backward selection
val_back <- c('patid' , 'ua_bacteria', 'ua_clarity', 'ua_color', 'ua_epi', 
              'ua_leuk', 'ua_nitrite', 'ua_urobili', 'ua_wbc', 
              'CVA_tenderness', 'abd_distended2', 'gen_neg', 
              'psychiatric_confusion', 'flank_pain', 'dysuria', 'hematuria', 
              'chief_complaint', 'age', 'gender', 'ethnicity_Non_Hispanic', 
              'UTI_diag')

# Compare model fit using AIC and BIC
AIC(glm_back)  # 43571.79
BIC(glm_back)  # 43829.49

AIC(glm3)   # 40798.76
BIC(glm3)   # 41838.44

# Keep only the selected variables in train and test sets
train <- train_glm[, val_back]
test <- test_glm[, val_back]
validation <- validation_glm[, val_back]

# Fit logistic regression on the final refined set of predictors
glm4 <- glm(UTI_diag ~ ., family = binomial, data = train)
summary(glm4)

# Model comparison using ANOVA
anova(glm3, glm2)  # Compare glm3 (with transformed categorical variables) to glm2
anova(glm4, glm3)  # Compare glm4 (refined model) to glm3

# Model evaluation
evaluate_glm(glm4, train, validation)  # AUC=0.8402
evaluate_glm(glm4, train, test)  # AUC=0.8438

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
summary(glm5)

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
summary(glm6)

# Perform another backward selection to ensure model simplicity
glm_back <- step(glm6, test = 'LRT')  # Automatic feature selection

# **Final Decision:** 
# Retain CVA_tenderness for interpretability, despite model suggesting removal.

# **Evaluation of Final Model**
# Achieved AUC of **0.8117** with a streamlined model—strong performance.
evaluate_glm(glm6, train, test)

# Top by GLM:
# 'patid', 'ua_bacteria', 'ua_clarity', 'ua_epi', 
# 'ua_wbc', 'CVA_tenderness', 'psychiatric_confusion', 
# 'flank_pain', 'age', 'gender', 
# 'ethnicity_Non_Hispanic', 'UTI_diag'


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

# Print AUC Values
print(paste("Validation AUC-ROC:", round(auc(roc_valid), 4)))  # 0.8106
print(paste("Testing AUC-ROC:", round(auc(roc_test), 4)))  # 0.8117

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

# Print PR AUC Values
print(paste("Validation AUC-PR:", round(pr_valid$auc.integral, 4)))  # 0.6349
print(paste("Testing AUC-PR:", round(pr_test$auc.integral, 4)))  # 0.6352

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
cal_curve_valid <- val.prob(validation$predicted_prob, y_valid_numeric, pl = TRUE)
title(main = "Calibration Curve for Validation (val.prob)")
plot_calibration_curve(validation$predicted_prob, y_valid_numeric, "Validation")

cal_curve_test <- val.prob(test$predicted_prob, y_test_numeric, pl = TRUE)
title(main = "Calibration Curve for Testing (val.prob)")
plot_calibration_curve(test$predicted_prob, y_test_numeric, "Testing")


################################################################################
# Random Forest Model Implementation
################################################################################

# Using 12 Key Variables Identified from Logistic Regression
# (instead of the 11 variables selected in the final glm model)
val_step <- c('patid', 'ua_bacteria', 'ua_clarity', 'ua_epi', 'ua_leuk',
              'ua_nitrite', 'ua_wbc', 'CVA_tenderness', 'psychiatric_confusion',
              'flank_pain', 'age', 'gender', 'UTI_diag')

# Extract the selected features for training, validation, and testing datasets
train <- train_all[, val_step]
validation <- validation_all[, val_step]   # Add validation dataset
test <- test_all[, val_step]

# Convert UTI diagnosis into a factor for classification
train$UTI_diag <- as.factor(as.character(train$UTI_diag))
validation$UTI_diag <- as.factor(as.character(validation$UTI_diag))
test$UTI_diag <- as.factor(as.character(test$UTI_diag))

# Function to handle missing values in categorical variables
handle_missing_values <- function(x) {
  char_vars <- sapply(x, is.character)  # Identify character columns
  x[, char_vars][is.na(x[, char_vars])] <- 'not_reported'  # Replace NAs with 'not_reported'
  return(x)
}

# Apply missing value handling
train <- handle_missing_values(train)
validation <- handle_missing_values(validation)  # Handle missing values in validation set
test <- handle_missing_values(test)

# Training the Random Forest Model
rf <- randomForest(UTI_diag ~ ., data = train, ntree = 500, norm.votes = FALSE, 
                   do.trace = 10, importance = TRUE)

# Display model summary
print(rf)

# Extract and visualize variable importance
# Evaluating Feature Importance 
importance(rf)

# Generate the variable importance plot
varImpPlot(rf, 
           main = "Variable Importance Plot - UTI Diagnosis",   # Title of the plot
           col.main = "blue",                                   # Title color
           font.main = 2,                                       # Make title bold
           cex.main = 1.5,                                      # Title size
           pch = 16,                                            # Point style
           col = "darkred",                                     # Point color
           cex.axis = 1.2,                                      # Axis label size
           cex.lab = 1.0,                                       # Axis labels size
           las = 2,                                             # Axis label orientation
           mar = c(5, 10, 4, 2)                                 # Adjust margins to fit labels
)

## Based on variable importance, we might consider removing 'CVA_tenderness' and/or 'psychiatric_confusion'.

# Top 10 - Random Forrest:
# 'patid', 'ua_bacteria', 'ua_clarity', 'ua_leuk', 'ua_nitrite', 'ua_wbc', 'age',
# 'ua_epi', 'CVA_tenderness', 'flank_pain'

## Evaluate the model performance on the validation set
validation_predictions <- predict(rf, validation)
validation_confusion_matrix <- table(validation$UTI_diag, validation_predictions)
validation_accuracy <- sum(diag(validation_confusion_matrix)) / sum(validation_confusion_matrix)
print(paste("Validation Accuracy: ", round(validation_accuracy, 4)))  # 0.8273

## Evaluate the model performance on the test set
test_predictions <- predict(rf, test)
test_confusion_matrix <- table(test$UTI_diag, test_predictions)
test_accuracy <- sum(diag(test_confusion_matrix)) / sum(test_confusion_matrix)
print(paste("Test Accuracy: ", round(test_accuracy, 4)))   # 0.83

################################################################################
# AUC-ROC Curve for Validation and Testing Data (Random Forest Model)
################################################################################

# Predict probabilities for validation and test data
validation$predicted_prob <- predict(rf, validation, type = "prob")[, 2]
test$predicted_prob <- predict(rf, test, type = "prob")[, 2]

# Compute ROC Curves
roc_valid_rf <- roc(validation$UTI_diag, validation$predicted_prob)
roc_test_rf <- roc(test$UTI_diag, test$predicted_prob)

# Plot ROC Curves
par(mfrow = c(1, 2))  # Side-by-side plots
plot(roc_valid_rf, col = "blue", lwd = 2, main = "AUC-ROC Curve (Validation)")
legend("bottomright", legend = paste("AUC =", round(auc(roc_valid_rf), 4)), 
       col = "blue", lwd = 2)

plot(roc_test_rf, col = "red", lwd = 2, main = "AUC-ROC Curve (Testing)")
legend("bottomright", legend = paste("AUC =", round(auc(roc_test_rf), 4)), 
       col = "red", lwd = 2)

# Print AUC Values
print(paste("Validation AUC-ROC:", round(auc(roc_valid_rf), 4))) # 0.8755
print(paste("Testing AUC-ROC:", round(auc(roc_test_rf), 4)))  # 0.8738

# The AUC-ROC values indicate that the Random Forest model has good discriminatory ability, 
# with both validation and test scores above 0.87. The minimal difference between validation 
# and test AUC-ROC suggests that the model generalizes well without significant overfitting. 
# However, compared to the XGBoost model (which had AUC-ROC values above 0.93), the Random 
# Forest model appears to be slightly less effective in classification.

################################################################################
# Precision-Recall AUC for Validation and Testing Data (Random Forest Model)
################################################################################

# Convert target variable to numeric
y_valid_numeric_rf <- as.numeric(validation$UTI_diag) - 1
y_test_numeric_rf <- as.numeric(test$UTI_diag) - 1

# Compute PR Curves
pr_valid_rf <- pr.curve(scores.class0 = validation$predicted_prob[y_valid_numeric_rf == 0],  
                        scores.class1 = validation$predicted_prob[y_valid_numeric_rf == 1],  
                        curve = TRUE)

pr_test_rf <- pr.curve(scores.class0 = test$predicted_prob[y_test_numeric_rf == 0],  
                       scores.class1 = test$predicted_prob[y_test_numeric_rf == 1],  
                       curve = TRUE)

# Plot PR Curves
par(mfrow = c(1, 2))
plot(pr_valid_rf, col = "blue", main = "Precision-Recall Curve (Validation)", lwd = 2)
legend("bottomright", legend = paste("AUC-PR =", round(pr_valid_rf$auc.integral, 4)), 
       col = "blue", lwd = 2)

plot(pr_test_rf, col = "red", main = "Precision-Recall Curve (Testing)", lwd = 2)
legend("bottomright", legend = paste("AUC-PR =", round(pr_test_rf$auc.integral, 4)), 
       col = "red", lwd = 2)

# Print PR AUC Values
print(paste("Validation AUC-PR:", round(pr_valid_rf$auc.integral, 4)))  # 0.6134
print(paste("Testing AUC-PR:", round(pr_test_rf$auc.integral, 4)))  # 0.6136

# The AUC-PR values are quite similar between validation and test sets, which indicates stable
# model performance across different data splits. The AUC-PR score of ~0.61 suggests that the 
# model maintains a reasonable balance between precision and recall. Interestingly, the Random
# Forest model's AUC-PR is slightly better than that of XGBoost (~0.59 in the results below), 
# which may indicate that it performs better in terms of handling class imbalance or capturing
# relevant patterns in the data.

################################################################################
# Calibration Curve for Validation and Testing Data (Random Forest Model)
################################################################################

# Function to compute and plot calibration curve with handling for duplicate bin cutoffs
plot_calibration_curve <- function(pred_probs, true_labels, dataset_name) {
  # Compute calibration curve manually
  bin_cutoffs <- quantile(pred_probs, probs = seq(0, 1, length.out = 11))  # 10 bins
  
  # Remove duplicates from bin_cutoffs
  bin_cutoffs <- unique(bin_cutoffs)
  
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

# Plot Calibration Curves for Validation and Testing (Random Forest)
cal_curve_valid <- val.prob(validation$predicted_prob, y_valid_numeric_rf, pl = TRUE)
title(main = "Calibration Curve for Validation (val.prob)")
plot_calibration_curve(validation$predicted_prob, y_valid_numeric_rf, "Validation - RF")

cal_curve_test <- val.prob(test$predicted_prob, y_test_numeric_rf, pl = TRUE)
title(main = "Calibration Curve for Testing (val.prob)")
plot_calibration_curve(test$predicted_prob, y_test_numeric_rf, "Testing - RF")


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

# # Convert the train, validation and test data into XGBoost matrix type
xgboost_train <- xgb.DMatrix(data=X_train, label=y_train)
xgboost_validation <- xgb.DMatrix(data=X_validation, label=y_validation)
xgboost_test <- xgb.DMatrix(data=X_test, label=y_test)

# Train XGBoost Model
model <- xgboost(data = xgboost_train,           # The data   
                 max.depth=3,                    # Max depth 
                 nrounds=50,
                 objective = "binary:logistic")  # Binary classification with logistic loss
summary(model)

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
print(conf_mat_valid)
print(conf_mat_test)

# Accuracy on Validation and Test
validation_accuracy <- sum(diag(conf_mat_valid$table)) / sum(conf_mat_valid$table)
test_accuracy <- sum(diag(conf_mat_test$table)) / sum(conf_mat_test$table)

print(paste("Validation Accuracy: ", round(validation_accuracy, 4)))  # 0.889
print(paste("Test Accuracy: ", round(test_accuracy, 4)))  # 0.896

################################################################################
# XGBoost Feature Importance Analysis
################################################################################

# Compute feature importance matrix
importance_matrix <- xgb.importance(colnames(xgboost_train), model = model)

# Print importance matrix
print(importance_matrix)

# Adjust y-axis label position
par(mgp = c(3, 0, 0))

# Plot top 15 important features
xgb.plot.importance(importance_matrix[1:15,], 
                    main = "Top 15 Feature Importance - XGBoost",
                    xlab = "Importance Score",
                    ylab = "Features")

# Reset par to default after plotting
par(mgp = c(3, 1, 0))

# Plot the multi-tree structure of the model
xgb.plot.multi.trees(model = model, use.names = FALSE, fill = TRUE)

# Top 10 XGBoost
# 'abx', 'patid', 'ua_wbc', 'ua_leuk', 'ua_nitrite', 
# 'ua_bacteria', 'dispo', 'antibiotics', 'age', 'chief_complaint'

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

# Print AUC Values
print(paste("Validation AUC-ROC:", round(auc(roc_valid), 4)))  # 0.9348
print(paste("Test AUC-ROC:", round(auc(roc_test), 4)))  # 0.9564

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

# Print PR AUC Values
print(paste("Validation AUC-PR:", round(pr_valid$auc.integral, 4)))  # 0.5934
print(paste("Test AUC-PR:", round(pr_test$auc.integral, 4)))  # 0.5878

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

# Plot Calibration Curves for Validation and Testing
cal_curve_valid <- val.prob(pred_valid, y_validation, pl = TRUE)
title(main = "Calibration Curve for Validation (val.prob)")
plot_calibration_curve(pred_valid, y_validation, "Validation - XGBoost")

cal_curve_test <- val.prob(pred_test, y_test, pl = TRUE)
title(main = "Calibration Curve for Testing (val.prob)")
plot_calibration_curve(pred_test, y_test, "Testing - XGBoost")

################################################################################

# Most common 10 variables:
# 'patid', 'ua_bacteria', 'ua_clarity', 'ua_wbc', 'age',
# 'dispo', 'antibiotics', 'ua_leuk', 'abx', 'ua_nitrite'

# Write out post feature-selection data
vars <- c('patid', 'ua_bacteria', 'ua_clarity', 'ua_wbc', 'age',
          'dispo', 'antibiotics', 'ua_leuk', 'abx', 'ua_nitrite')
write.csv(train_all_just_in_case[, vars], 'finaldat.csv')

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# Exploratory data analysis (Cont'd)

# Boxplot of Age by Gender
p1 <- ggplot(df, aes(x = gender, y = age, fill = gender)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Age vs. Gender", x = "Gender", y = "Age", fill="Gender") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Boxplot of Age by Clarity
p2 <- ggplot(df, aes(x = ua_clarity, y = age, fill = ua_clarity)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Age vs. Urine Clarity", x = "Clarity", y = "Age", fill="Clarity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Boxplot of Age by Ethnicity
p3 <- ggplot(df, aes(x = ethnicity, y = age, fill = ethnicity)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Age vs. Ethnicity", x = "Ethnicity", y = "Age", fill="Ethnicity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Boxplot of Age by WBC
p4 <- ggplot(df, aes(x = as.factor(ua_wbc), y = age, fill = ua_wbc)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(title = "Age vs. WBC Count", x = "WBC Count", y = "Age", fill="WBC") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Arrange all plots in a single panel
grid.arrange(p1, p2, p3, p4, ncol = 2)


# Bar plot of WBC by Gender
p1 <- ggplot(df, aes(x = as.factor(ua_wbc), fill = gender)) +
  geom_bar(position = "dodge") +
  labs(title = "WBC vs. Gender", x = "WBC", y = "Count", fill="Gender") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Bar plot of WBC by Clarity
p2 <- ggplot(df, aes(x = as.factor(ua_wbc), fill = ua_clarity)) +
  geom_bar(position = "dodge") +
  labs(title = "WBC vs. Clarity", x = "WBC", y = "Count", fill="Clarity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Bar plot of WBC by Ethnicity
p3 <- ggplot(df, aes(x = as.factor(ua_wbc), fill = ethnicity)) +
  geom_bar(position = "dodge") +
  labs(title = "WBC vs. Ethnicity", x = "WBC", y = "Count", fill="Ethnicity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Arrange plots in a single panel
grid.arrange(p1, p2, p3, ncol = 2)


# Bar plot of Bacteria by Gender
p1 <- ggplot(df, aes(x = ua_bacteria, fill = gender)) +
  geom_bar(position = "dodge") +
  labs(title = "Bacteria vs. Gender", x = "Bacteria", y = "Count", fill="Gender") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Bar plot of Bacteria by Clarity
p2 <- ggplot(df, aes(x = ua_bacteria, fill = ua_clarity)) +
  geom_bar(position = "dodge") +
  labs(title = "Bacteria vs. Clarity", x = "Bacteria", y = "Count", fill="Clarity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Bar plot of Bacteria by Ethnicity
p3 <- ggplot(df, aes(x = ua_bacteria, fill = ethnicity)) +
  geom_bar(position = "dodge") +
  labs(title = "Bacteria vs. Ethnicity", x = "Bacteria", y = "Count", fill="Ethnicity") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line()
  )

# Arrange plots in a single panel
grid.arrange(p1, p2, p3, ncol = 2)

