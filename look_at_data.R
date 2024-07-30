library(haven)
library(ggplot2)
library(readxl)
library(patchwork)
library(plotly)

# from here: https://www.thelancet.com/journals/lanpub/article/PIIS2468-2667(22)00060-3/fulltext

data <- read_dta('covid_data/Data_from_OxCGRT.dta')

# Filter rows where country is "United Kingdom"
uk_data <- subset(data, country == "United Kingdom")

pop_data <- read_dta('covid_data/World_population_world_bank.dta')

uk_pop <- pop_data[251, 2]

# Display the filtered data
print(uk_data)

# Convert date column to Date type if necessary
uk_data$date_mdy <- as.Date(uk_data$date_mdy, format="%m/%d/%Y")

# Initialize episimdata with relevant columns
uk_data$daily_cases_100k[is.na(uk_data$daily_cases_100k)] <- 0
episimdata <- data.frame(date_mdy = uk_data$date_mdy, C = (uk_data$daily_cases_100k), Lambda_C = rep(1, nrow(uk_data)), R_coeff = rep(1, nrow(uk_data)), Rest = rep(1, nrow(uk_data)), R0est = rep(1, nrow(uk_data)))

# Set initial variables
gen_time <- 6.5
gen_time_var <- 2.1

R_est_wind <- 5  # Define your window size
Ygen <- dgamma(1:nrow(uk_data), gen_time/gen_time_var, 1/gen_time_var)
Ygen <- Ygen/sum(Ygen)

# Estimate R
for (ii in 2:nrow(episimdata)) {
  if (ii-1 < R_est_wind) {
    episimdata[ii, 'Rest'] <- mean(episimdata[1:(ii-1), 'C']) / mean(episimdata[1:(ii-1), 'Lambda_C'])
    R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'R_coeff']) / sum(Ygen[1:(ii-1)])
  } else {
    if ( mean(episimdata[(ii-R_est_wind):(ii-1), 'Lambda_C']) == 0){
      episimdata[ii, 'Rest'] <- 0
      R_coeff_tmp <- 1
    } else {
      episimdata[ii, 'Rest'] <- mean(episimdata[(ii-R_est_wind):(ii-1), 'C']) / mean(episimdata[(ii-R_est_wind):(ii-1), 'Lambda_C'])
      R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'R_coeff']) / sum(Ygen[1:(ii-1)])
    }
  }
  episimdata[ii, 'R0est'] <- episimdata[ii, 'Rest'] / R_coeff_tmp
  episimdata[ii, 'Lambda_C'] <- sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])
  episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
}

uk_data <- uk_data[uk_data$date_mdy >= as.Date("2020-02-21") & uk_data$date_mdy <= as.Date("2021-10-01"), ]

uk_data$daily_cases_100k <- round((uk_data$daily_cases_100k)*67215293/1e5, digits =  0)

mob_data <- read_excel("covid_data/googlemobilitydataset290922.xlsx", sheet = "Mobility data Index 7 day MA")
mob_data$Date <- as.Date(mob_data$Date)
mob_data <- mob_data[mob_data$Date >= as.Date("2020-02-21") & mob_data$Date <= as.Date("2021-10-01"), ]

episimdata <- episimdata[episimdata$date_mdy >= as.Date("2020-02-21") & episimdata$date_mdy <= as.Date("2021-10-01"), ]

ggplot(data = mob_data, aes(x = Date)) +
  geom_line(aes(y = `Retail and recreation`, color = "Retail and recreation")) +
  geom_line(aes(y = `Grocery and pharmacy`, color = "Grocery and pharmacy")) +
  geom_line(aes(y = Parks, color = "Parks")) +
  geom_line(aes(y = `Transit stations`, color = "Transit stations")) +
  geom_line(aes(y = Workplaces, color = "Workplaces")) +
  geom_line(aes(y = Residential, color = "Residential")) +
  labs(title = "Mobility Trends Over Time",
       x = "Date",
       y = "Value",
       color = "Legend") +
  theme_minimal()


mobility_plot <- ggplot(data = mob_data, aes(x = Date)) +
  geom_line(aes(y = `Retail and recreation`, color = "Retail and recreation")) +
  geom_line(aes(y = `Grocery and pharmacy`, color = "Grocery and pharmacy")) +
  geom_line(aes(y = `Transit stations`, color = "Transit stations")) +
  geom_line(aes(y = Workplaces, color = "Workplaces")) +
  labs(title = "Mobility Trends Over Time",
       x = "Date",
       y = "Value",
       color = "Legend") +
  theme_minimal()

# Plot for daily cases
cases_plot <- ggplot(uk_data, aes(x = date_mdy, y = daily_cases_100k)) +
  geom_line(color = "blue") +
  labs(title = "Daily COVID-19 Cases in the United Kingdom",
       x = "Date",
       y = "Daily Cases UK") +
  theme_minimal()

# Plot for stringency index
stringency_plot <- ggplot(uk_data, aes(x = date_mdy, y = stringencyindex)) +
  geom_line(color = "red") +
  labs(title = "Stringency Index in the United Kingdom",
       x = "Date",
       y = "Stringency Index") +
  theme_minimal()

# Combine the plots using patchwork
combined_plot <- mobility_plot / cases_plot / stringency_plot
print(combined_plot)
# Plot date_mdy and stringencyindex

ggplot(uk_data, aes(x = date_mdy)) +
  geom_line(aes(y = daily_cases_100k/1000, color = "Daily Cases per 100k")) +
  geom_line(aes(y = stringencyindex, color = "Stringency Index")) +
  scale_y_continuous(
    name = "Daily Cases*1000",
    sec.axis = sec_axis(~., name = "Stringency Index")
  ) +
  labs(title = "Daily COVID-19 Cases and Stringency Index in the United Kingdom",
       x = "Date",
       color = "Legend") +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "red"),
    axis.title.y.left = element_text(color = "blue"))

# Display the resulting episimdata
print(episimdata)
Rplot <- ggplot(episimdata, aes(x = date_mdy, y = R0est)) +
  geom_line(color = "green") +
  labs(title = "Estimated R Values in the United Kingdom",
       x = "Date",
       y = "Estimated R") +
  theme_minimal()

combined_plot <- mobility_plot / stringency_plot / Rplot

# Print the combined plot
print(combined_plot)

### redo with plotly

mobility_plot <- plot_ly(mob_data, x = ~Date) %>%
  add_lines(y = ~`Retail and recreation`, name = "Retail and recreation", line = list(color = 'blue')) %>%
  add_lines(y = ~`Grocery and pharmacy`, name = "Grocery and pharmacy", line = list(color = 'red')) %>%
  add_lines(y = ~`Transit stations`, name = "Transit stations", line = list(color = 'purple')) %>%
  add_lines(y = ~Workplaces, name = "Workplaces", line = list(color = 'orange')) %>%
  layout(title = "Mobility Trends Over Time",
         xaxis = list(title = "Date"),
         yaxis = list(title = "Value"),
         legend = list(title = list(text = "Legend")))

# Plot for daily cases
Rplot <- plot_ly(episimdata, x = ~date_mdy, y = ~R0est, type = 'scatter', mode = 'lines',
                      line = list(color = 'blue'), name = 'Rest') %>%
  layout(title = "Daily COVID-19 Cases in the United Kingdom",
         xaxis = list(title = "Date"),
         yaxis = list(title = "Daily Cases per 100k"))

# Plot for stringency index
stringency_plot <- plot_ly(uk_data, x = ~date_mdy, y = ~stringencyindex, type = 'scatter', mode = 'lines',
                           line = list(color = 'red'), name = 'Stringency Index') %>%
  layout(title = "Stringency Index in the United Kingdom",
         xaxis = list(title = "Date"),
         yaxis = list(title = "Stringency Index"))

# Combine the plots using subplot from plotly
combined_plot1 <- subplot(mobility_plot, Rplot, stringency_plot, nrows = 3, shareX = TRUE)

# Display the combined plot
combined_plot1

# Merge episimdata and uk_data on date_mdy if needed
merged_data <- merge(episimdata, uk_data, by = "date_mdy")

# Create the plot
scatter_plot <- plot_ly(merged_data, x = ~stringencyindex, y = ~R0est, type = 'scatter', mode = 'markers',
                        marker = list(color = 'blue')) %>%
  layout(title = "Rest vs Stringency Index",
         xaxis = list(title = "Stringency Index"),
         yaxis = list(title = "Rest"))

# Display the plot
scatter_plot

merged_data2 <- merge(mob_data, merged_data, by.x = "Date", by.y = "date_mdy")

retail_plot <- plot_ly(merged_data2, x = ~R0est, y = ~`Retail and recreation`, type = 'scatter', mode = 'markers',
                       marker = list(color = 'blue')) %>%
  layout(title = "Retail and Recreation vs Rest",
         xaxis = list(title = "Rest"),
         yaxis = list(title = "Retail and Recreation"))

# Plot Grocery and pharmacy vs Stringency Index
grocery_plot <- plot_ly(merged_data2, x = ~R0est, y = ~`Grocery and pharmacy`, type = 'scatter', mode = 'markers',
                        marker = list(color = 'red')) %>%
  layout(title = "Grocery and Pharmacy vs Stringency Index",
         xaxis = list(title = "Rest"),
         yaxis = list(title = "Grocery and Pharmacy"))

# Plot Transit stations vs Stringency Index
transit_plot <- plot_ly(merged_data2, x = ~R0est, y = ~`Transit stations`, type = 'scatter', mode = 'markers',
                        marker = list(color = 'purple')) %>%
  layout(title = "Transit Stations vs Stringency Index",
         xaxis = list(title = "Rest"),
         yaxis = list(title = "Transit Stations"))

# Plot Workplaces vs Stringency Index
workplaces_plot <- plot_ly(merged_data2, x = ~R0est, y = ~Workplaces, type = 'scatter', mode = 'markers',
                           marker = list(color = 'orange')) %>%
  layout(title = "Workplaces vs Stringency Index",
         xaxis = list(title = "Rest"),
         yaxis = list(title = "Workplaces"))

# Combine the plots using subplot from plotly
combined_plot2 <- subplot(retail_plot, grocery_plot, transit_plot, workplaces_plot, nrows = 2, shareX = TRUE, shareY = FALSE)

# Display the combined plot
combined_plot2

retail_plot3 <- plot_ly(merged_data2, x = ~stringencyindex, y = ~`Retail and recreation`, type = 'scatter', mode = 'markers',
                       marker = list(color = 'blue')) %>%
  layout(title = "Retail and Recreation vs Rest",
         xaxis = list(title = "Stringency Index"),
         yaxis = list(title = "Retail and Recreation"))

# Plot Grocery and pharmacy vs Stringency Index
grocery_plot3 <- plot_ly(merged_data2, x = ~stringencyindex, y = ~`Grocery and pharmacy`, type = 'scatter', mode = 'markers',
                        marker = list(color = 'red')) %>%
  layout(title = "Grocery and Pharmacy vs Stringency Index",
         xaxis = list(title = "Stringency Index"),
         yaxis = list(title = "Grocery and Pharmacy"))

# Plot Transit stations vs Stringency Index
transit_plot3 <- plot_ly(merged_data2, x = ~stringencyindex, y = ~`Transit stations`, type = 'scatter', mode = 'markers',
                        marker = list(color = 'purple')) %>%
  layout(title = "Transit Stations vs Stringency Index",
         xaxis = list(title = "Stringency Index"),
         yaxis = list(title = "Transit Stations"))

# Plot Workplaces vs Stringency Index
workplaces_plot3 <- plot_ly(merged_data2, x = ~stringencyindex, y = ~Workplaces, type = 'scatter', mode = 'markers',
                           marker = list(color = 'orange')) %>%
  layout(title = "Workplaces vs Stringency Index",
         xaxis = list(title = "Stringency Index"),
         yaxis = list(title = "Workplaces"))

combined_plot3 <- subplot(retail_plot3, grocery_plot3, transit_plot3, workplaces_plot3, nrows = 2, shareX = TRUE, shareY = FALSE)
combined_plot3

