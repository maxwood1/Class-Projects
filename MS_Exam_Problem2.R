#Applied Stat MS Exam Problem 2
#Max Woodbury

library(readxl)
library(dplyr)
library(tibble)
library(reshape2)


## read in data

setwd("C:/Users/maxwo/OneDrive/Desktop/STAT582")

lifedata <- read_excel("LifeData.xlsx", sheet=1)
num_sheets <- length(excel_sheets("LifeData.xlsx"))

for(i in 2:num_sheets) {
  temp <- read_excel("LifeData.xlsx", sheet=i)
  lifedata <- lifedata %>% bind_rows(temp)
}

colnames(lifedata)[1:2] <- c("Country_Code","Country_Name")


## look at data

View(lifedata %>% filter(Indicator %in% c("Female life expectancy at Birth", "Male life Expectancy at Birth",
                                          "Life Expectancy at Birth (Total)")))


## data transformation: flip rows and columns

lifedata <- lifedata[,-1]
lifedata <- melt(lifedata)
colnames(lifedata)[3] <- "Year"
lifedata$Year <- as.numeric(as.character(lifedata$Year))



## Life Expectancy at Birth (Total)


## year vs life expectancy (total) for all countries

year_life <- lifedata %>% filter(Indicator == "Life Expectancy at Birth (Total)")
countries <- unique(lifedata$Country_Name)

colors <- rainbow(16)
count=1
for(country in countries) {
  country_life <- year_life %>% filter(Country_Name==country)
  
  missing <- which(is.na(country_life$value))
  if (length(missing) > 0) {
    country_life <- country_life[-which(is.na(country_life$value)),]
  }
  
  if(country=="Argentina") {
    plot(country_life$Year, country_life$value, type='l', xlab="Year", ylab="Life Expectancy at Birth (Total)", 
         lwd=2, xlim=c(1946,2012), ylim=c(20,100))
  }
  else {
    lines(country_life$Year, country_life$value, type='l', lwd=2, col=colors[count])
  }
  count=count+1
}


## linear regression for overall year vs total life expectancy

year_life.lm <- lm(value ~ Year, data=year_life)
summary(year_life.lm)
#note linear reg intercept is at year 0 which doesn't make sense to interpret

plot(fitted(year_life.lm), resid(year_life.lm), xlab="Fitted Values", ylab="Residuals")
#assumptions of linear model violated because of non-constant variance

#weighted regression to even out residuals (standardized)
yl <- year_life[-which(is.na(year_life$value)),]
wt <- 1 / lm(abs(year_life.lm$residuals) ~ year_life.lm$fitted.values)$fitted.values^2
wt.lm <- lm(value ~ Year, data = yl, weights=wt)
summary(wt.lm)

#plot standardize residuals to check constant variance assumption
stdres <- rstandard(wt.lm)
plot(wt.lm$fitted.values, stdres)

#main idea: increasing trend in life expectancy over time



## which countries have the largest increase in life expectancy?

#get slopes for each individual countries' regression line
slopes <- c()

for(country in countries) {
  country_life <- year_life %>% filter(Country_Name==country)
  
  country.lm <- summary(lm(value ~ Year, data=country_life))
  slopes <- slopes %>% bind_rows(data.frame(country, country.lm$coefficients[2,1]))
}
colnames(slopes) <- c("Country","Slope")

View(slopes)


#two different trends
#separate countries into two groups by slopes

#low life exp countries
low <- c("China","India","Mexico","Thailand","Turkey", "Philippines")
#high life exp countries
high <- setdiff(countries, low)

#low life exp countries
colors <- rainbow(length(low))
count=1
for(country in low) {
  country_life <- year_life %>% filter(Country_Name==country)
  
  if(country==low[1]) {
    plot(country_life$Year, country_life$value, type='l', xlab="Year", ylab="Life Expectancy at Birth (Total)", 
         lwd=2, xlim=c(1946,2012), ylim=c(20,100))
  }
  else {
    lines(country_life$Year, country_life$value, type='l', lwd=2, col=colors[count])
  }
  count=count+1
}

#plot low increase countries
colors <- rainbow(length(high))
count=1
for(country in high) {
  country_life <- year_life %>% filter(Country_Name==country)
  
  if(country==high[1]) {
    plot(country_life$Year, country_life$value, type='l', xlab="Year", ylab="Life Expectancy at Birth (Total)", 
         lwd=2, xlim=c(1946,2012), ylim=c(20,100))
  }
  else {
    lines(country_life$Year, country_life$value, type='l', lwd=2, col=colors[count])
  }
  count=count+1
}

#two regression lines

yl_low <- year_life %>% filter(Country_Name %in% low)
low.lm <- lm(value ~ Year, data=yl_low)
summary(low.lm)
plot(fitted(low.lm), resid(low.lm))

yl_high <- year_life %>% filter(Country_Name %in% high)
high.lm <- lm(value ~ Year, data=yl_high)
summary(high.lm)
plot(fitted(high.lm), resid(high.lm))

#show these on new plot with color grouping

for(country in countries) {
  country_life <- year_life %>% filter(Country_Name==country)
  missing <- which(is.na(country_life$value))
  if (length(missing) > 0) {
    country_life <- country_life[-which(is.na(country_life$value)),]
  }
  
  if(country %in% high) {
    color="black"
  }
  else {
    color="red"
  }
  
  if(country=="Argentina") {
    plot(country_life$Year, country_life$value, type='l', xlab="Year", ylab="Life Expectancy at Birth (Total)", 
         lwd=1.5, xlim=c(1946,2012), ylim=c(20,100), col=color)
  }
  else {
    lines(country_life$Year, country_life$value, type='l', lwd=1.5, col=color)
  }
}

abline(low.lm, lwd=4,col="red")
abline(high.lm, lwd=4)



## ANCOVA for between country comparison of total life expectancy

year_life$Country_Name <- as.factor(year_life$Country_Name)

#added year as covariate to explain the variability over time and increase country significance
year_life.aov <- aov(value ~ Year + Country_Name, data=year_life)
summary(year_life.aov)

#there is a significant difference between the average life expectancy of the countries



#look at other variables


#take averages over all observed years for each country and plot these points for each country to look for trends

plot_life_avg <- function(lifedata, var) {
  test1 <- c()
  
  test <- lifedata %>% filter(Indicator %in% c(var, "Life Expectancy at Birth (Total)"))
  for(i in seq(1, nrow(test), by=2)){
    if(!is.na(test$value[i]) & !is.na(test$value[i+1])) {
      test1 <- bind_rows(test1,test[c(i:(i+1)),])
    }
    i=i+2;
  }
  
  life_subset <- test1 %>% filter(Indicator == "Life Expectancy at Birth (Total)")
  var_subset <- test1 %>% filter(Indicator == var)
  
  life_countries <- unique(life_subset$Country_Name)
  var_countries <- unique(var_subset$Country_Name)
  countries <- intersect(life_countries, var_countries)
  
  life_subset <- life_subset %>% filter(Country_Name %in% countries)
  var_subset <- var_subset %>% filter(Country_Name %in% countries)
  
  avg_data <- c()
  
  for(country in countries) {
    country_life <- life_subset %>% filter(Country_Name==country)
    country_var <- var_subset %>% filter(Country_Name==country)
    
    life_years <- unique(country_life$Year)
    var_years <- unique(country_var$Year)
    years <- intersect(life_years, var_years)
    
    country_life <- country_life %>% filter(Year %in% years)
    country_var <- country_var %>% filter(Year %in% years)
    
    country_data <- as.data.frame(cbind(country_var$value, country_life$value))
    colnames(country_data) <- c("x", "life")
    
    avgx <- mean(country_data$x)
    avglife <- mean(country_data$life)
    
    avg_data <- bind_rows(avg_data,as.data.frame(cbind(country, avgx, avglife)))
  }
  
  avg_data$avgx <- as.numeric(as.character(avg_data$avgx))
  avg_data$avglife <- as.numeric(as.character(avg_data$avglife))
  
  maxi <- max(avg_data$avgx,na.rm=T)
  mini <- min(avg_data$avgx,na.rm=T)
  
  #make box plot
  for(i in 1:nrow(avg_data)){
    if(avg_data$country[i] %in% high) {
      color="black"
    } else {
      color="red"
    }
    if(i==1) {
      plot(avg_data$avgx[i], avg_data$avglife[i], xlab=var, ylab="Life Expectancy at Birth (Total)", lwd=2, 
           xlim=c(mini,maxi), ylim=c(40,100), col=color)
    } else {
      points(avg_data$avgx[i], avg_data$avglife[i], lwd=2, col=color)
    }
  }

}


#plot all variables
for(var in unique(lifedata$Indicator)) {
  plot_life_avg(lifedata, var)
}

#plot interesting variables 
variables <- c("Average Years of Education","Gender Equality Years of Education", "Educational Inequality Gini Coefficient",
               "Working week in manufacturing", "Global Extreme Poverty Dollar a Day", "Composite Measure of Wellbeing", 
               "GDP per Capita")
for(var in variables) {
  plot_life_avg(lifedata, var)
}

#interesting clustering in addition to linear trends

#these plots above show how factors impact life expectancy between countries



#can do some analysis over time too
#like the countries that increase the most in life exp
#also increase the most in ... over time

allvars <- unique(lifedata$Indicator)

similar <- c("Educational Inequality Gini Coefficient","Gender Equality Years of Education","Global Extreme Poverty Dollar a Day")

countries <- unique(lifedata$Country_Name)

for(var in similar) {
  year_var <- lifedata %>% filter(Indicator == var)
  
  for(country in countries) {
    country_life <- year_var %>% filter(Country_Name==country)
    missing <- which(is.na(country_life$value))
    if (length(missing) > 0) {
      country_life <- country_life[-which(is.na(country_life$value)),]
    }
    
    maxi <- max(year_var$value,na.rm=T)
    mini <- min(year_var$value,na.rm=T)
    
    if(country %in% high) {
      color="black"
    } else {
      color="red"
    }
    
    if(country=="Argentina") {
      plot(country_life$Year, country_life$value, type='l', xlab="Year", ylab=var, 
           lwd=1.5, xlim=c(1946,2012), ylim=c(mini,maxi), col=color)
    }
    else {
      lines(country_life$Year, country_life$value, type='l', lwd=1.5, col=color)
    }
  }
}

#Similar trend of decreased variance over time for these factors like the trend of life expectancy
