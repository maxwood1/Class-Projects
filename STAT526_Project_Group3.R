# Spring 2021
# STAT 526: Advanced Statistical Methodology
# Group 3 Final Project
# Predictive Modeling for Home Values in Buffalo, NY
# Used Linear Regression and Linear Mixed Models with Random Effects

packages <- c("readr","dplyr","janitor","lubridate","ggplot2","ggcorrplot","generalhoslem","nlme")
lapply(packages, require, character.only = TRUE)


### Get data ready for analysis ###

buffalo <- read_csv("Buffalo_Data.csv") %>% clean_names()

#convert deed dates to R date format
buffalo$deed_date <- as.POSIXct(buffalo$deed_date, format="%m/%d/%Y")

#subset by single-family dwellings with a deed date in the second half of 2018
#include only important variables
buffalo.sub <- buffalo %>% 
          filter(sale_price>1, deed_date >= as.Date("2018-07-01"), property_class_description=="ONE FAMILY DWELLING") %>%
          dplyr::select(sale_price, total_value, front, depth, year_built, total_living_area, overall_condition, number_of_beds, number_of_baths, zip_code_5_digit)

#Factorize zip code
buffalo.sub$zip_code_5_digit <- as.factor(buffalo.sub$zip_code_5_digit)

nrow(buffalo.sub)


#### Which house attributes most significantly affect sale price and tax appraisal respectively? ####

#First do correlation plot for numeric variables to see which ones are interesting or problematic
cordata <- buffalo.sub %>% 
  dplyr::select(sale_price, total_value, front, depth, year_built, total_living_area, overall_condition, number_of_beds, number_of_baths)

cm <- cor(cordata, method="pearson")
ggcorrplot(cm)
# Number of baths - living area corr .69 but still keeping both


#box plot sale prices per zip code to see what they look like
plot(buffalo.sub$zip_code_5_digit, buffalo.sub$sale_price)


ggplot(buffalo.sub, aes(x=zip_code_5_digit,y=sale_price))+
  geom_boxplot(fill="#ff7f7f", lwd=1, fatten=1, outlier.size=2.75) +
  theme(axis.title=element_text(size=14))

#box plot of sale price and tax appraisal per zip code

which(buffalo.sub$sale_price == max(buffalo.sub$sale_price))

#View(buffalo.sub[411,])
#View(buffalo[which(buffalo$sale_price == 3350000),])

#drop max observation
buffalo.sub <- buffalo.sub[-411,]



################# MODEL ONE #################

#Linear Regression: sale price
linsalemod <-lm(sale_price ~  front + depth + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths + zip_code_5_digit,data = buffalo.sub)
summary(linsalemod)

#take out insignificant predictors with stepwise (depth)
linsalestep <-step(linsalemod, trace = F)
summary(linsalestep)
#sig: front, year_built, total_living_area, overall_condition, number_of_beds, number_of_baths, zip_code_5_digit

#Diagnostics

#Box Cox
boxcox(linsalestep)

#Residuals
#cone shaped, but this can be somewhat explained
plot(fitted(linsalestep),resid(linsalestep),xlab="Fitted",ylab="Residuals",main="Linear Model Sale Price Residuals")

#Q-Q plot
plot(linsalestep, which=2)
qqnorm(resid(linsalestep),main="Linear Model Sale Price Q-Q Plot")
qqline(resid(linsalestep),col=2)

#Histogram of Resid with normal curve
linsaledf <- as.data.frame(resid(linsalestep))
names(linsaledf) <- c("x")

ggplot(linsaledf, aes(x)) +
  geom_histogram(aes(y = ..density..), bins=20) +
  stat_function(fun = dnorm,
                args = list(mean = mean(linsaledf$x), sd = sd(linsaledf$x)),
                col = "#1b98e0",
                size=1.5) +
  xlab("Residuals") +
  ggtitle("Sale Price Histogram of Residuals")

#We don't have predictors to adequately show the effect of expensive house prices (trim and stuff)
#R-squared good



################# MODEL TWO #################

#Linear mixed model sale price

#not interested in any specific zip code effect, but still want to control for it

#linear mixed model for sale price
mixmod <- lme(sale_price ~ front + depth + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixmod)

#get rid of insignificant predictors (depth)
mixmod1 <- lme(sale_price ~ front + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixmod1)
#same predictors as linear model

#use LRT to test if reduced model is better
pchisq(2*as.numeric(logLik(mixmod)-logLik(mixmod1)),3,lower=F)
#fail to reject, so smaller model is better
#parametric bootstrapping not needed here since p-values here tend to be small

#Use LRT to test signficance of zip code random effect
mixmod2 <- lm(sale_price ~ front + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, data=buffalo.sub)
pchisq(2*as.numeric(logLik(mixmod1)-logLik(mixmod2)),1,lower=F)
#parametric bootstrapping not needed here either
#p-values tend to be large, so a significant p-value is actually significant

#Diagnostics

#Residuals
qqnorm(resid(mixmod1),main="Mixed Linear Model Sale Price Q-Q Plot")
qqline(resid(mixmod1), col = 2)
plot(fitted(mixmod1),resid(mixmod1),xlab="Fitted",ylab="Residuals",main="Mixed Linear Model Sale Price Residuals")

#Random effects
qqnorm(ranef(mixmod1)[[1]],main="Sale Price Zip Code Random Effects")
qqline(ranef(mixmod1)[[1]],col=2)
#14202 causing issues

#Log-Transform sale price mixed

mixmodlog <- lme(log(sale_price) ~ front + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, 
                 random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixmodlog)
mixmodlog1 <- lme(log(sale_price) ~ front + year_built + total_living_area + overall_condition + number_of_baths, 
                  random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixmodlog1)
mixmodlog2 <- lme(log(sale_price) ~ front + year_built + total_living_area + overall_condition, 
                  random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixmodlog2)

#Log-transformed Residuals
qqnorm(resid(mixmodlog2),main="Log-Transformed Mixed Model Sale Price Q-Q Plot")
qqline(resid(mixmodlog2), col = 2)
plot(fitted(mixmodlog2),resid(mixmodlog2),xlab="Fitted",ylab="Residuals",main="Log-Transformed Mixed Model Sale Price Residuals")

#Log-transformed Random Effects
qqnorm(ranef(mixmodlog2)[[1]],main="Log-Transformed Sale Price Zip Code Random Effects")
qqline(ranef(mixmodlog2)[[1]],col=2)



################# MODEL THREE #################

#Linear Regression: Tax Appraisal
linappmod <-lm(total_value ~  front + depth + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths + zip_code_5_digit,data = buffalo.sub)
summary(linappmod)

#take out insignificant predictors with stepwise (kept depth this time)
linappstep <-step(linappmod, trace = F)
summary(linappstep)

#Diagnostics

#Box Cox
boxcox(linappstep)

#Residuals
#cone shaped, but this can be somewhat explained
plot(fitted(linappstep),resid(linappstep),xlab="Fitted",ylab="Residuals",main="Linear Model Tax Appraisal Residuals")

#Q-Q plot
plot(linappstep, which=2)
qqnorm(resid(linappstep),main="Linear Model Tax Appraisal Q-Q Plot")
qqline(resid(linappstep),col=2)

#Histogram of Resid with normal curve
linappdf <- as.data.frame(resid(linappstep))
names(linappdf) <- c("x")

ggplot(linappdf, aes(x)) +
  geom_histogram(aes(y = ..density..), bins=20) +
  stat_function(fun = dnorm,
                args = list(mean = mean(linappdf$x), sd = sd(linappdf$x)),
                col = "#1b98e0",
                size=1.5)



################# MODEL FOUR #################

#linear mixed model for tax appraisal
mixappmod <- lme(total_value ~ front + depth + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixappmod)

#get rid of insignificant predictors (depth)
mixappmod1 <- lme(total_value ~ front + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixappmod1)
#took out depth this time (different)

#use LRT to test if reduced model is better
pchisq(2*as.numeric(logLik(mixappmod)-logLik(mixappmod1)),3,lower=F)
#fail to reject, so smaller model is better
#parametric bootstrapping not needed here since p-values here tend to be small

#Use LRT to test signficance of zip code random effect
mixappmod2 <- lm(total_value ~ front + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, data=buffalo.sub)
pchisq(2*as.numeric(logLik(mixappmod1)-logLik(mixappmod2)),1,lower=F)
#parametric bootstrapping not needed here either
#p-values tend to be large, so a significant p-value is actually significant

#Diagnostics

#Residuals
qqnorm(resid(mixappmod1), main="Mixed Linear Model Tax Appraisal Q-Q Plot")
qqline(resid(mixappmod1), col = 2)
plot(fitted(mixappmod1),resid(mixappmod1),xlab="Fitted",ylab="Residuals",main="Mixed Linear Model Tax Appraisal Residuals")

#Random effects
qqnorm(ranef(mixappmod1)[[1]],main="Tax Appraisal Zip Code Random Effects")
qqline(ranef(mixappmod1)[[1]],col=2)
#14202 causing issues

#Log-Transform Tax Appraisal Mixed

mixapplog <- lme(log(total_value) ~ front + year_built + total_living_area + overall_condition + number_of_beds + number_of_baths, 
                 random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixapplog)
mixapplog1 <- lme(log(total_value) ~ front + year_built + total_living_area + overall_condition + number_of_baths, 
                  random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixapplog1)
mixapplog2 <- lme(log(total_value) ~ front + year_built + total_living_area + overall_condition, 
                  random=~1|zip_code_5_digit, method="ML", data=buffalo.sub)
summary(mixapplog2)

#Residuals
qqnorm(resid(mixapplog2),main="Log-Transformed Mixed Model Tax Appraisal Q-Q Plot")
qqline(resid(mixapplog2), col = 2)
plot(fitted(mixapplog2),resid(mixapplog2),xlab="Fitted",ylab="Residuals",main="Log-Transformed Mixed Model Tax Appraisal Residuals")

#Random Effects
qqnorm(ranef(mixapplog2)[[1]],main="Log-Transformed Tax Appraisal Zip Code Random Effects")
qqline(ranef(mixapplog2)[[1]],col=2)




