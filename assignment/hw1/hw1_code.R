# Load necessary libraries
library(dplyr)
library(omopr)
library(tidyverse)

################################################################################
# Database Connection Setup
################################################################################

# Setup for database
db_host <- rstudioapi::askForPassword("Hostname: ")
db_user <- rstudioapi::askForPassword("Username: ")
db_pass <- rstudioapi::askForPassword("Database password")
db_name <- rstudioapi::askForPassword("Database name: ")

# Connect to database
db_conn <- DBI::dbConnect(RPostgres::Postgres(),
                          host = db_host,
                          dbname = db_name,
                          user = db_user,
                          password = db_pass,
                          options="-c search_path=omop")

# Initiate omopr/dplyr connection
table_refs = omopr_init(db_conn)

# Clean up workspace
rm(db_host, db_pass, db_user, db_name)

################################################################################
# Extract Relevant Tables
################################################################################

# Visit Occurrence: Contains information about each patient visit.
# Alternative: Could use `visit_detail` table for more granular details, but it may increase complexity.
visit_data <- table_refs[["visit_occurrence"]] %>% concept_names()

# Person: Contains demographic information about patients.
# Alternative: Could use `observation` table for additional demographic details, but it may introduce redundancy.
person_data <- table_refs[["person"]] %>% concept_names()

# Condition Occurrence: Contains information about patient conditions.
# Alternative: Could use `observation` or `procedure` tables, but they may not capture conditions as comprehensively.
condition_data <- table_refs[["condition_occurrence"]] %>% concept_names()

# Measurement: Contains measurements (e.g., lab results, vitals) for patients.
# Alternative: Could use `observation` table, but it may not include structured measurement data.
measurement_data <- table_refs[["measurement"]] %>% concept_names('measurement_concept_id')

# Death: Contains death records for patients.
# Alternative: Could use `observation` table for death records, but it may not be as standardized.
death_data <- table_refs[["death"]] %>% select(person_id, death_date)

# Function to remove concept_id columns and retain concept names
get_concept_names <- function(data) {
  # Identify columns containing 'concept_id' (excluding gender-related columns)
  concept_cols <- which(grepl('concept_id', colnames(data)) & !grepl('gender', colnames(data)))
  # Remove concept_id columns to simplify the dataset
  new_data <- data %>% select(-colnames(data)[concept_cols])
  return(new_data)
}

# Apply the function to remove concept_id columns
visit_data <- get_concept_names(visit_data)
person_data <- get_concept_names(person_data)
condition_data <- get_concept_names(condition_data)

################################################################################
# Select Relevant Columns
################################################################################

# Columns from the visit_occurrence table
visit_cols <- c("visit_occurrence_id", "person_id", "visit_start_datetime", 
                "visit_end_datetime", "provider_id", "visit_concept_name", 
                "visit_type_concept_name",  "visit_source_concept_name", 
                "admitted_from_concept_name", "discharged_to_concept_name", 
                "preceding_visit_occurrence_id")

# Columns from the person table
person_cols <- c("person_id", "gender_concept_id", "birth_datetime", 
                 "gender_source_value", "ethnicity_source_value", "race_source_value")

# Columns from the condition_occurrence table
condition_cols <- c("condition_occurrence_id", "person_id", "provider_id",
                    "condition_start_datetime", "condition_end_datetime",
                    "visit_occurrence_id", "condition_concept_name")

# Columns from the measurement table
measurement_cols <- c("measurement_id", "person_id", "value_as_number", 
                      "provider_id", "visit_detail_id", "value_source_value",
                      "unit_source_value", "visit_occurrence_id",
                      "measurement_concept_id", "measurement_concept_name")

# Select only the relevant columns from each table
visit_data <- visit_data %>% select(all_of(visit_cols))
person_data <- person_data %>% select(all_of(person_cols))
condition_data <- condition_data %>% select(all_of(condition_cols))
measurement_data <- measurement_data %>% select(all_of(measurement_cols))

################################################################################
# Merge and Collect Data
################################################################################

# Join person and death data on person_id to include death information
person_death_data <- person_data %>% inner_join(death_data, 'person_id')

# Join visit data with person and death data on person_id
joined_data <- visit_data %>% inner_join(person_death_data, 'person_id')

# Join condition and measurement data on person_id, visit_occurrence_id, and provider_id
final_data <- reduce(list(joined_data, condition_data, measurement_data), 
                     inner_join, 
                     by = c('person_id', 'visit_occurrence_id', 'provider_id'))

# Clean up workspace
rm(visit_data, person_data, condition_data, measurement_data, death_data, 
   visit_cols, person_cols, condition_cols, measurement_cols, joined_data, 
   get_concept_names)

# Columns selected for the requirement in HW1
selected_cols <- c('person_id', 'visit_occurrence_id', 
                   'visit_start_datetime', 'visit_end_datetime', 'birth_datetime', 
                   'gender_concept_id', 'race_source_value', 'ethnicity_source_value',
                   'value_as_number', 'condition_concept_name',
                   'measurement_concept_name', 'preceding_visit_occurrence_id',
                   "admitted_from_concept_name", "discharged_to_concept_name")

# Collect the data into memory for further processing
collected_data <- final_data %>% select(all_of(selected_cols)) %>% collect()

################################################################################
# Apply Inclusion and Exclusion Criteria
################################################################################

# The following code applies inclusion and exclusion criteria to the dataset.

processed_data <- collected_data %>% 
  # Create new columns
  mutate(age_at_visit_days = difftime(visit_start_datetime, birth_datetime, 
                                      units = 'days'),
         visit_length_days = difftime(visit_end_datetime, visit_start_datetime, 
                                      units = 'days'),
         has_diabetes = grepl('diabetes',
                              tolower(condition_concept_name)),
         has_hypertension = grepl('hypertension',
                                  tolower(condition_concept_name))) %>%
  # Convert age_at_visit_days to years and visit_length_days to numeric
  mutate(age_at_visit_years = as.numeric(gsub(' .*', '', 
                                              age_at_visit_days))/365,
         visit_length_days = as.numeric(gsub(' .*', '', 
                                             visit_length_days))) %>% 
  # Arrange data by person_id and visit_start_datetime to ensure chronological order
  arrange(person_id, visit_start_datetime) %>%
  # Create a column to identify new visits for each patient
  mutate(is_new_visit = !(duplicated(person_id) & 
                            duplicated(visit_occurrence_id))) %>% 
  group_by(person_id) %>% 
  # Create a cumulative visit number for each patient
  mutate(visit_number = cumsum(is_new_visit)) %>% 
  # Calculate the total number of visits for each patient
  mutate(total_visits = max(visit_number)) %>% 
  # Apply inclusion/exclusion criteria:
  # 1. Include visits with visit_number between 1 and 37
  # 2. Include patients aged 52 to 56 years
  # 3. Exclude patients with 3-17 total visits and systolic blood pressure > 140
  filter(visit_number > 0,
         visit_number < 38,
         age_at_visit_years >= 52,
         age_at_visit_years <= 56,
         !((total_visits %in% 3:17) & 
             !(tolower(measurement_concept_name) == 'systolic blood pressure' & 
                 value_as_number > 140))) %>% 
  # Remove duplicate measurements for the same visit
  distinct(person_id, visit_number, measurement_concept_name, .keep_all = TRUE) %>% 
  # Select and rename columns for the final dataset
  select(visit_date = 'visit_start_datetime',
         person_id = 'person_id',
         gender_concept_id = 'gender_concept_id',
         age_at_visit = 'age_at_visit_years',
         race = 'race_source_value',
         ethnicity = 'ethnicity_source_value',
         visit_length = 'visit_length_days',
         has_hypertension,
         has_diabetes,
         total_visits,
         admitted_from_concept_name,
         discharged_to_concept_name,
         measurement_name = 'measurement_concept_name',
         measurement_value = 'value_as_number') %>%
  # Pivot measurements into wide format for easier analysis
  pivot_wider(names_from = measurement_name, values_from = measurement_value, 
              values_fill = NA)

################################################################################
# Final Processing and Export
################################################################################

# Remove rows without any associated measurements
num_cols <- ncol(processed_data)
rows_without_measurements <- apply(processed_data[, 13:num_cols], 1, 
                                   function(x) sum(is.na(x))) == length(13:num_cols)
processed_data <- processed_data[!rows_without_measurements, ]

# Convert visit_date to character format for consistency
processed_data$visit_date <- as.character(as.POSIXct(processed_data$visit_date))

# Select only the required columns for Part 1 of the homework
final_part1_data <- processed_data %>% 
  select(-total_visits, -admitted_from_concept_name, -discharged_to_concept_name)

# Save the processed data
saveRDS(processed_data, 'data.RDS')
write.csv(final_part1_data, 'dat.csv')

################################################################################
################################################################################
################################################################################
################################################################################

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

################################################################################
# Bar Plots/Histograms for Descriptive Statistics
################################################################################

# 1. Patient Age at Visit (in years)
# Bar plot with automatic age binning (width = 1 year)
p1_age <- ggplot(processed_data, aes(x = cut_width(
  age_at_visit, width = 1, boundary = floor(min(processed_data$age_at_visit))))) +
  geom_bar(fill = "#0072B2", alpha = 0.8) +
  labs(title = 'Age at Visit for Cohort',
       subtitle = 'Age is grouped into 1-year intervals',
       x = 'Age at Visit (Years)',
       y = 'Count') +
  scale_x_discrete(labels = function(x) {
    gsub("\\((\\d+)\\.(\\d+),\\s*(\\d+)\\.(\\d+)\\]", "(\\1,\\3]", x)
  }) +  # Format labels to remove decimals
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

# 2. Gender Distribution (8532 = Female, 8507 = Male)
p2_gender <- ggplot(processed_data, aes(x = factor(gender_concept_id, 
                                                   levels = c(8532, 8507),
                                                   labels = c('Female', 'Male')))) +
  geom_bar(fill = "#0072B2", alpha = 0.8) +
  labs(title = 'Gender Distribution for Cohort',
       x = 'Gender',
       y = 'Count') +
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

# 3. Race Distribution
p3_race <- ggplot(processed_data, aes(x = factor(race))) +
  geom_bar(fill = "#0072B2", alpha = 0.8) +
  labs(title = 'Race Distribution for Cohort',
       x = 'Race',
       y = 'Count') +
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

# 4. Visit Length (in days)
p4_visit_length <- ggplot(processed_data, aes(x = visit_length)) +
  geom_bar(fill = "#0072B2", alpha = 0.8) +
  labs(title = 'Visit Length for Cohort',
       subtitle = 'Visit length is split by days',
       x = 'Visit Length (Days)',
       y = 'Count') +
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
part1_plots <- list(p1_age, p2_gender, p3_race, p4_visit_length)

# Display the plots
grid.arrange(grobs = part1_plots, ncol = 2)

# Description of what we learn from the bar plots/histograms:
# - The age distribution shows that the age at visit is roughly uniformly distributed.
# - The gender distribution shows a higher proportion of males, indicating a gender imbalance in the cohort.
# - The race distribution highlights that the cohort is predominantly white, with limited representation from other racial groups.
# - Visit lengths are typically very short, with most visits lasting 0 days, though some extend up to 20 days.

################################################################################
# Box Plots for Age at Visit by Different Variables
################################################################################

# 5. Age at Visit vs. Total Number of Visits
p5_age_visits <- ggplot(processed_data, aes(y = age_at_visit)) +
  geom_boxplot(aes(fill = factor(total_visits)), alpha = 0.7) +
  facet_grid(. ~ total_visits) +
  ylab('Age at Visit (Years)') +
  xlab('Visit Number') +
  labs(title = 'Age at Visit by Total Number of Visits',
       subtitle = 'Faceted by Total Number of Visits') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", color = "black"),   # Customizing the text style
    strip.background = element_rect(color = "gray", fill = "lightgray", size = 1, linetype = "solid"),  # Shadowed box for facet labels
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(),
    panel.grid = element_blank()
  )

# 6. Age at Visit vs. Gender
p6_age_gender <- ggplot(processed_data, aes(y = age_at_visit)) +
  geom_boxplot(aes(fill = factor(gender_concept_id, 
                                 levels = c(8532, 8507),
                                 labels = c('Female', 'Male'))), alpha = 0.7) +
  facet_grid(. ~ factor(gender_concept_id, 
                        levels = c(8532, 8507),
                        labels = c('Female', 'Male'))) +
  scale_fill_brewer(palette = "Set1") +
  ylab('Age at Visit (Years)') +
  xlab('Gender') +
  labs(title = 'Age at Visit by Gender',
       subtitle = 'Faceted by Gender') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", color = "black"),   # Customizing the text style
    strip.background = element_rect(color = "gray", fill = "lightgray", size = 1, linetype = "solid"),  # Shadowed box for facet labels
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(),
    panel.grid = element_blank()
  )

# 7. Age at Visit vs. Admitted From
p7_age_admitted <- ggplot(processed_data, aes(y = age_at_visit)) +
  geom_boxplot(aes(fill = factor(admitted_from_concept_name, levels = unique(admitted_from_concept_name))),
               alpha = 0.7) +
  facet_grid(. ~ admitted_from_concept_name) +
  scale_fill_manual(values = rep("gray40", length(unique(processed_data$admitted_from_concept_name)))) +  # Use gray for all boxes
  ylab('Age at Visit (Years)') +
  xlab('Admitted Location') +
  labs(title = 'Age at Visit by Admitted From',
       subtitle = 'Faceted by Admitted Location') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", color = "black"),   # Customizing the text style
    strip.background = element_rect(color = "gray", fill = "lightgray", size = 1, linetype = "solid"),  # Shadowed box for facet labels
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(),
    panel.grid = element_blank()
  )

# 8. Age at Visit vs. Discharged To
p8_age_discharged <- ggplot(processed_data, aes(y = age_at_visit)) +
  geom_boxplot(aes(fill = factor(discharged_to_concept_name, levels = unique(discharged_to_concept_name))),
               alpha = 0.7) +
  facet_grid(. ~ discharged_to_concept_name) +
  scale_fill_manual(values = rep("gray40", length(unique(processed_data$discharged_to_concept_name)))) +  # Use gray for all boxes
  ylab('Age at Visit (Years)') +
  xlab('Discharged Location') +
  labs(title = 'Age at Visit by Discharged To',
       subtitle = 'Faceted by Discharged Location') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", color = "black"),   # Customizing the text style
    strip.background = element_rect(color = "gray", fill = "lightgray", size = 1, linetype = "solid"),  # Shadowed box for facet labels
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(),
    panel.grid = element_blank()
  )

# 9. Age at Visit vs. Race
p9_age_race <- ggplot(processed_data, aes(y = age_at_visit)) +
  geom_boxplot(aes(fill = race), alpha = 0.7) +
  facet_grid(. ~ race) +
  scale_fill_brewer(palette = "Set2") +
  ylab('Age at Visit (Years)') +
  xlab('Race') +
  labs(title = 'Age at Visit by Race',
       subtitle = 'Faceted by Race') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", color = "black"),   # Customizing the text style
    strip.background = element_rect(color = "gray", fill = "lightgray", size = 1, linetype = "solid"),  # Shadowed box for facet labels
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(),
    panel.grid = element_blank()
  )

# 10. Age at Visit vs. Ethnicity
p10_age_ethnicity <- ggplot(processed_data, aes(y = age_at_visit)) +
  geom_boxplot(aes(fill = ethnicity), alpha = 0.7) +
  facet_grid(. ~ ethnicity) +
  scale_fill_brewer(palette = "PiYG") +
  ylab('Age at Visit (Years)') +
  xlab('Ethnicity') +
  labs(title = 'Age at Visit by Ethnicity',
       subtitle = 'Faceted by Ethnicity') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", color = "black"),   # Customizing the text style
    strip.background = element_rect(color = "gray", fill = "lightgray", size = 1, linetype = "solid"),  # Shadowed box for facet labels
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(),
    panel.grid = element_blank()
  )

# Combine all box plots into a single list
part2_plots <- list(p5_age_visits, p6_age_gender, p7_age_admitted, 
                    p8_age_discharged, p9_age_race, p10_age_ethnicity)

# Display the plots in a 3x2 grid
grid.arrange(grobs = part2_plots, ncol = 2, nrow = 3)

# Description of what we learn from the box plots:
# - Hispanic and male patients are generally younger, while female and non-Hispanic patients are older.
# - Black patients are notably older than both Asian and White patients.
# - The middle two plots provide no useful information about the locations where patients are admitted or discharged.

################################################################################
# Save Plots
################################################################################

# Save the plots for future use
saveRDS(part1_plots, 'Plots1.RDS')
saveRDS(part2_plots, 'Plots2.RDS')

