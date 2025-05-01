library(EpiControl)
library(VGAM)
library(parallel)
library(pbapply)
library(zoo)  # For rolling sum operations
library(ggplot2)
library(EpiEstim)

library(readr)
library(dplyr)
library(tidyr)

# Read data
ebola_data <- read_csv("code_to_test_functionality/rstb20160308_si_001.csv")

# Extract and convert date column
date_reported <- as.Date(ebola_data$DateReport)  # Ensure it's in Date format
date_inferred <- as.Date(ebola_data$DateOnsetInferred)

# Create a sequence of all dates from min to max
all_dates <- data.frame(date_reported = seq(min(date_reported, na.rm = TRUE),
                                            max(date_reported, na.rm = TRUE),
                                            by = "day"))

# Count occurrences of each reported date
date_counts <- as.data.frame(table(date_reported, useNA = "no"))
date_counts_i <- as.data.frame(table(date_inferred, useNA = "no"))

colnames(date_counts) <- c("date_reported", "count")
colnames(date_counts_i) <- c("date_reported", "count")
date_counts$date_reported <- as.Date(date_counts$date_reported)
date_counts_i$date_reported <- as.Date(date_counts_i$date_reported)# Convert back to Date format

# Merge full date range with counts, filling missing dates with 0
date_counts_filled <- all_dates %>%
  left_join(date_counts, by = "date_reported") %>%
  replace_na(list(count = 0))

date_counts_i_filled <- all_dates %>%
  left_join(date_counts_i, by = "date_reported") %>%
  replace_na(list(count = 0))

# Print the complete frequency table
print(date_counts_filled)
print(date_counts_i_filled)

date_counts_filled$count_i <- date_counts_i_filled$count

ggplot(date_counts_filled, aes(x = date_reported, y = count)) +
  geom_line(color = "blue") +  # Line plot for cases over time
  geom_line(color = "red", aes(x = date_reported, y = count_i)) +
  labs(title = "Ebola Cases Over Time",
       x = "Date",
       y = "Number of Cases") +
  theme_minimal()

cases_v <- date_counts_filled$count_i

Ygen <- dgamma(1:ndays, 15.0/2.1, 1/2.1)
Ygen <- Ygen/sum(Ygen)
Ygen <- c(0, Ygen)

res <- estimate_R(incid = cases_v,
                  method = "non_parametric_si",
                  config = make_config(si_distr = Ygen))


plot(res)
