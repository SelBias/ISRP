
# This file aims to reproduce results of "A Poisson Autoregressive Model to Understand COVID-19 Contagion Dynamics", written by Agosto and Giudici in 2020. 

library(dplyr)
library(ggplot2)

# Step 1. Load data
# Data are downloaded from https://covid19.who.int/data at 2023.08.25 11PM KST
global_countries = read.csv("./covid_WHO/WHO-COVID-19-global-data.csv")
global_countries = global_countries[, c("Date_reported", "Country", "New_cases", "Cumulative_cases", "New_deaths", "Cumulative_deaths")]
global_countries$Date_reported = as.Date(global_countries$Date_reported)

selected_countries = global_countries %>% filter(Country %in% c("Republic of Korea", "China", "Italy", "Iran (Islamic Republic of)"))

# south_korea = global_countries %>% filter(Country == "Republic of Korea")
# china = global_countries %>% filter(Country == "China")
# italy = global_countries %>% filter(Country == "Italy")
# iran = global_countries %>% filter(Country == "Iran (Islamic Republic of)")

# Step 2. Reproduction of Figure 1 
# It seems that there is a huge difference on data for China during Feb 10 ~ 15 2020. 
# In my guess, these data might be updated after the publication of the paper. 
fig_1 = selected_countries %>% 
    filter(Date_reported < "2020-04-02") %>% 
    filter(Date_reported > "2020-01-20") %>% 
    ggplot() + 
    theme_bw() + 
    geom_line(mapping = aes(x = Date_reported, y = New_cases, group = Country, color = Country))

ggsave(fig_1, filename = "./fig_COVID_Poisson_AR/fig_1.png", width = 10, height = 6)






# us_counties = read.csv("./covid-19-data-master/us-counties.csv")
# santa_cruz = us_counties %>% filter(county == "Santa Cruz") %>% filter(state == "California")
# santa_cruz = santa_cruz[, -c(2,3,4)]

# cumulative_santa_cruz = data.frame(ts(santa_cruz))
# daily_santa_cruz = data.frame(cbind(
#     cumulative_santa_cruz$date[-1], 
#     lapply(cumulative_santa_cruz, diff, lag = 1)$cases, 
#     lapply(cumulative_santa_cruz, diff, lag = 1)$deaths
# ))

# head(daily_santa_cruz)
# colnames(daily_santa_cruz) = c("date", "cases", "deaths")

# ggplot(daily_santa_cruz) + 
#     theme_bw() + 
#     geom_line(mapping = aes(x = date, y = cases), color = "black") + 
#     geom_line(mapping = aes(x = date, y = deaths), color = "red")
