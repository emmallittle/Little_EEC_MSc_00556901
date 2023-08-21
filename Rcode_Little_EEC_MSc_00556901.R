#clean up space and set directory
rm(list=ls())
setwd("C:/Users/emmalittle/Desktop/Imperial Masters/Research Project/Submission")
list.files()
tempdir()
dir.create(tempdir())
options("install.lock"=FALSE)
R.version.string
citation()
require(dplyr)         # for general data wrangling and plotting
require(ggplot2)
require(tidyr)
require(viridis)
require(ggtext)
require(psych)         # for pairs function testing collinearity
require(usdm)          # for vif function testing collinearity
require(MASS)          # for glm.nb
require(lme4)          # for glmer
require(DHARMa)        # for testing model performance
require(performance)   # for model performance
require(ggeffects)     # for ggpredict to investigate model predictions
require(lattice)       # for heatmap of interaction
require(metR)          # for heatmap of interaction
require(sjPlot)        # for effect size charts
require(stargazer)     # for table of results


#################################################################################
#################################################################################
#######################   EMPIRICAL STUDY  ######################################
#################################################################################
#################################################################################



## RESPONSE VARIABLE

# Create column of bat passes, bat feeding buzzes and bat activity
require(dplyr)
bat_data <- read.csv("bat_data.csv", header=TRUE)
bat_data <- mutate(bat_data, bat_passes = bat_data$com_pip_passes + 
                     bat_data$sop_pip_passes + bat_data$nat_pip_passes +
                     bat_data$noid_pip_passes + bat_data$nls_passes +
                     bat_data$myotis_passes + bat_data$unknown_passes) 
bat_data <- mutate(bat_data, bat_fb = bat_data$com_pip_fb + bat_data$sop_pip_fb +
                     bat_data$noid_pip_fb)
bat_data <- mutate(bat_data, bat_activity = if_else(bat_data$bat_passes==0,0,1))

# sum bat passes and feeding buzzes by species by sample point by date
daily_bat_data <- bat_data %>% 
  group_by(Sample_point, Start_date, .add = TRUE) %>%
  dplyr::select(bat_passes, bat_fb, bat_activity, com_pip_passes, sop_pip_passes, 
                nat_pip_passes, noid_pip_passes, myotis_passes, nls_passes, unknown_passes, 
                com_pip_fb, sop_pip_fb, noid_pip_fb, social_calls) %>% 
  summarize_all(sum)
daily_bat_data
daily_bat_data <- as.data.frame(daily_bat_data)

# extract number of minutes recorded per night
# minutes of data
results <- data.frame(Sample_point=character(), Start_date=character(), mins_rec=numeric())
results
for(z in unique(bat_data$Sample_point)){
  for(i in unique(bat_data$Start_date)) {
    temp_file <- bat_data %>% subset(bat_data$Start_date==i & bat_data$Sample_point==z)
    count <- length(unique(temp_file$Time))
    temp <- data.frame(Start_date=i, mins_rec=count, Sample_point=z)
    results <- rbind(results, temp)
  }
}

daily_bat_data <- merge(daily_bat_data, results, by.x=c('Sample_point', 'Start_date'), 
                   by.y=c('Sample_point', 'Start_date'))

# Exclude two rows of data that failed during night,  per min rate would overstate
# as bats most active at dusk
empirical_data <- daily_bat_data
empirical_data <- filter(empirical_data, mins_rec == 102)


## EXPLANATORY VARIABLES
# 1/ LUX
lux_data <- read.csv("light_measurement_data.csv", header=TRUE)
# filter rows only vertical hold
lux_data <- filter(lux_data, Lux_positioning == "Vertical")
# create average lux from two recordings
lux_data <- mutate (lux_data, lux = rowMeans(cbind(lux_data$Lux_Night_1,lux_data$Lux_Night_2)))
# select columns 
lux_data <- dplyr::select(lux_data, Sample_point, lux)
# merge with dataset
empirical_data <- merge(empirical_data, lux_data, 'Sample_point', 'Sample_point', all.x=TRUE)

# 2/ LANDSCAPE DATA
landscape_data <- read.csv("sample_points_landscape_data.csv", header=TRUE)
# select columns 
landscape_data <- dplyr::select(landscape_data, -Date, -Site, -Surveyor, -Tree_ID, -Long, -Lat, -Accuracy, -Comments)
# merge with dataset
empirical_data <- merge(empirical_data, landscape_data, 'Sample_point', 'Sample_point', all.x=TRUE)

# 3/ TREE DENSITY
tree_density_data <- read.csv("tree_density_data.csv", header=TRUE)
# convert cm to m
tree_density_data <- mutate(tree_density_data, dist_nearest_tree_m = tree_density_data$Distance_nearest_tree/100)
# extract mean distance to four quarters' trees
results <- data.frame(Sample_point=character(), mean_distance=numeric())
results
for(i in unique(tree_density_data$Sample_point)) {
  temp_file <- tree_density_data %>% subset(tree_density_data$Sample_point==i)
  average <- mean((temp_file$dist_nearest_tree_m))
  temp <- data.frame(Sample_point=i, mean_distance=average)
  results <- rbind(results, temp)
}
# calculate density at sample point of trees per hectare
results <- mutate (results, tree_density = 10000*(1/(results$mean_distance)^2))
# select wanted columns
tree_density_results <- dplyr::select(results, -mean_distance)
# merge with dataset
empirical_data <- merge(empirical_data, tree_density_results, 'Sample_point', 'Sample_point', all.x=TRUE)

# 4/ CANOPY COVER
cc_data <- read.csv("canopy_cover_data.csv", header=TRUE)
# calculate two gap values and the average
cc_data <- mutate (cc_data, gap_para = cc_data$Para_255_value / (cc_data$Para_0_value + cc_data$Para_255_value))
cc_data <- mutate (cc_data, gap_perp = cc_data$Perp_255_value / (cc_data$Perp_0_value + cc_data$Perp_255_value))
cc_data <- mutate (cc_data, gap_avg = rowMeans(cbind(cc_data$gap_para,cc_data$gap_perp)))
# extract the mean for each sample point
results <- data.frame(Sample_point=character(), mean_gap=numeric())
results
for(i in unique(cc_data$Sample_point)) {
  temp_file <- cc_data %>% subset(cc_data$Sample_point==i)
  average <- mean((temp_file$gap_avg))
  temp <- data.frame(Sample_point=i, mean_gap=average)
  results <- rbind(results, temp)
}
gap_data <- results
# merge with dataset
empirical_data <- merge(empirical_data, gap_data, 'Sample_point', 'Sample_point', all.x=TRUE)
# convert to canopy cover
empirical_data <- mutate(empirical_data, canopy_cover = 1-empirical_data$mean_gap)

## 5/ INVERTEBRATE COUNT DATA
invert_data <- read.csv("invertebrate_data.csv", header=TRUE)
# select wanted columns
invert_data <- dplyr::select(invert_data, Sample_point, Invert_count)
# merge with dataset
empirical_data <- merge(empirical_data, invert_data, 'Sample_point', 'Sample_point', all.x=TRUE)

## 6/ DAILY ENVIRONMENTAL DATA
environ_data <- read.csv("daily_environmental_data.csv", header=TRUE)
# filter rows for dates required
environ_data <- filter(environ_data, Date == "02/05/2023" | Date == "05/05/2023" |Date == "06/05/2023"  |Date == "07/05/2023"  )
# select columns for night data
environ_data <- dplyr::select(environ_data, -contains("day"), -Comment)
# merge with dataset
site_data <- dplyr::select(bat_data, Site, Sample_point)
site_data <- site_data[!duplicated(site_data), ] 
empirical_data <- merge(empirical_data, site_data, 'Sample_point', 'Sample_point', all.x=TRUE)
empirical_data <- merge(empirical_data, environ_data, by.x=c('Site', 'Start_date'), 
                        by.y=c('Site', 'Date'))

## 7/ AMPLITUDE
dB_data <- read.csv("dB_data.csv", header=TRUE)
# exclude repeated data when had interference issue
dB_data <- filter(dB_data, Comment != "constant noise at low frequency, all files")
dB_data <- filter(dB_data, Comment != "constant noise at low frequency, all files for three days")

# extract the minimum dB for each sample point
results <- data.frame(Sample_point=character(), min_dB=numeric())
results
for(i in unique(dB_data$Sample_point)) {
  temp_file <- dB_data %>% subset(dB_data$Sample_point==i)
  mindB <- min((temp_file$Rel_dB))
  temp <- data.frame(Sample_point=i, min_dB=mindB)
  results <- rbind(results, temp)
}
# merge with dB
dB_data <- merge(dB_data, results, 'Sample_point', 'Sample_point', all.x=TRUE)
# calculate amplitude
dB_data <- mutate(dB_data, amp = 10^((Rel_dB-min_dB)/20))
# calculate variance in amplitude by sample point
results <- data.frame(Sample_point=character(), mean_dB=numeric())
results
for(i in unique(dB_data$Sample_point)) {
  temp_file <- dB_data %>% subset(dB_data$Sample_point==i)
  average <- mean((temp_file$Rel_dB))
  temp <- data.frame(Sample_point=i, mean_dB=average)
  results <- rbind(results, temp)
}
sound_data <- results

# calculate variance in amplitude by sample point
results <- data.frame(Sample_point=character(), var_amp=numeric())
results
for(i in unique(dB_data$Sample_point)) {
  temp_file <- dB_data %>% subset(dB_data$Sample_point==i)
  varamp <- var((temp_file$amp))
  temp <- data.frame(Sample_point=i, var_amp=varamp)
  results <- rbind(results, temp)
}

sound_data <- merge(sound_data, results,'Sample_point', 'Sample_point', all.x=TRUE)

# calculate amplitude by sample point
min(sound_data$mean_dB)
sound_data <- mutate(sound_data, amp = 10^((mean_dB-min(sound_data$mean_dB))/20))

# merge with dataset
emp_data <- merge(empirical_data, sound_data, 'Sample_point', 'Sample_point', all.x=TRUE)



## Alter Site, Date and Sample Point to factors

str(emp_data)
emp_data$Sample_point <- as.factor(emp_data$Sample_point)
emp_data$Site <- as.factor(emp_data$Site)
emp_data$Start_date <- as.factor(emp_data$Start_date)
class(emp_data$Sample_point)

emp_data$total_bat_passes <- emp_data$bat_passes
colnames(emp_data)[colnames(emp_data) == 'bat_passes'] <- 'total_bat_passes'
colnames(emp_data)[colnames(emp_data) == 'bat_fb'] <- 'total_bat_fb'
colnames(emp_data)[colnames(emp_data) == 'bat_activity'] <- 'total_bat_activity'

# inspect the data
length(unique(emp_data$Sample_point))
sum(emp_data$mins_rec)
sum(emp_data$total_bat_activity)
sum(emp_data$total_bat_passes)
sum(emp_data$total_bat_fb)
sum(emp_data$social_calls)

sum(emp_data$com_pip_passes)
sum(emp_data$sop_pip_passes)
sum(emp_data$nls_passes)
sum(emp_data$myotis_passes)
sum(emp_data$noid_pip_passes)
sum(emp_data$nat_pip_passes)
sum(emp_data$unknown_passes)

sum(emp_data$com_pip_fb)
sum(emp_data$sop_pip_fb)
sum(emp_data$noid_pip_fb)

### MODELLING

## Explore collinearity in continuous explanatory variables
## Using pairs panels
# test: lux, amp, temperature, tree density, canopy cover, invert count
require(psych)
pairs.panels(emp_data[,c(19, 46, 41, 20, 34, 35)])
# all <=0.55

## Repeat using VIF
require(usdm)
vif(emp_data[,c(19, 46, 41, 20, 34, 35)])
# all below 2

## Check for zeroes
colSums(emp_data == 0, na.rm = T)/length(emp_data$total_bat_passes)
# total bat feeding buzzes 51% zeroes
# nls and myotis passes 66-69% zeroes
# com pip and sop pip feeding buzzes 62-66% zeroes
# social calls 29% zeroes
# no streetlamps within 100m still 55% zeroes


### INVESTIGATING TOTAL BAT PASSES
# Use GLM with poisson
M1 <- glm(total_bat_passes~scale(lux) + scale(amp) + 
            scale(lux)*scale(amp) +
            scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
            scale(Mean_temp_night), 
          data = emp_data, family = "poisson")
summary(M1)


# check the fit
sum(cooks.distance(M1)>1)
sum(cooks.distance(M1)>1)/length(emp_data$total_bat_passes)
# 14 outliers, 8% of the observations
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
7084.9/171
# The model is highly over-dispersed - 41x higher conditional variance to 
# conditional mean. 
# Possible reasons: the outliers, missing covariates / interactions (with canopy cover)
# zero-inflation: only two zeros in 134 observations

## Alternative - use neg binomial distribution of errors
require(MASS)
M2 <- glm.nb(total_bat_passes~scale(lux) + scale(amp) + 
               scale(lux)*scale(amp) +
               scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
               scale(Mean_temp_night), 
             data = emp_data)
summary(M2)
# pseudoR2 = 1 - (residual dev / null dev) - 28% of the variance in bat passes can be explained 
1-(199.73/278.02)

sum(cooks.distance(M2)>1)
# no outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
199.73/171
# The model is only very slightly overdispersed, 1.17x higher conditional variance to 
# conditional mean

# Use DHARMa package to test performance
require(DHARMa)
testResiduals(M2)
par(mfrow=c(2,2))
testDispersion(M2)
testQuantiles(M2)
testUniformity(M2)
testZeroInflation(M2)
# no dispersion issues, residuals in quantile deviations ok - model looks ok

## Run stripping out species as random effect, need to pivot longer the dataframe
### Adapting dataset for long data for by species glmm for total passes
# Create column of bat passes, bat feeding buzzes and bat activity
bat_data <- read.csv("bat_data.csv", header=TRUE)
bat_data <- mutate(bat_data, bat_passes = bat_data$com_pip_passes + 
                     bat_data$sop_pip_passes + bat_data$nat_pip_passes +
                     bat_data$noid_pip_passes + bat_data$nls_passes +
                     bat_data$myotis_passes + bat_data$unknown_passes) 
bat_data <- mutate(bat_data, bat_fb = bat_data$com_pip_fb + bat_data$sop_pip_fb +
                     bat_data$noid_pip_fb)
bat_data <- mutate(bat_data, bat_activity = if_else(bat_data$bat_passes==0,0,1))

# sum the passes by species by sample point and date
sum_data <- bat_data %>% 
  group_by(Sample_point, Start_date, .add = TRUE) %>%
  dplyr::select(com_pip_passes, sop_pip_passes, nat_pip_passes, noid_pip_passes, myotis_passes, nls_passes, 
                unknown_passes) %>% 
  dplyr::summarize_all(sum)
sum_data
sum_data <- as.data.frame(sum_data)


# extract number of minutes recorded per night
# minutes of data
results <- data.frame(Sample_point=character(), Start_date=character(), mins_rec=numeric())
results
for(z in unique(bat_data$Sample_point)){
  for(i in unique(bat_data$Start_date)) {
    temp_file <- bat_data %>% subset(bat_data$Start_date==i & bat_data$Sample_point==z)
    count <- length(unique(temp_file$Time))
    temp <- data.frame(Start_date=i, mins_rec=count, Sample_point=z)
    results <- rbind(results, temp)
  }
}

sum_data <- merge(sum_data, results, by.x=c('Sample_point', 'Start_date'), 
                  by.y=c('Sample_point', 'Start_date'))

sum_data_long <- filter(sum_data, mins_rec == 102)
sum_data_long <- dplyr::select(sum_data_long, -mins_rec)
require(tidyr)
sum_data_long <- pivot_longer(sum_data_long, names_to = "Species", values_to = "Passes", cols = 3:9)

# Add in explanatory variables from earlier
sum_data_long <- merge(sum_data_long, lux_data, 'Sample_point', 'Sample_point', all.x=TRUE)
sum_data_long <- merge(sum_data_long, landscape_data, 'Sample_point', 'Sample_point', all.x=TRUE)
sum_data_long <- merge(sum_data_long, tree_density_results, 'Sample_point', 'Sample_point', all.x=TRUE)
sum_data_long <- merge(sum_data_long, gap_data, 'Sample_point', 'Sample_point', all.x=TRUE)
# convert to canopy cover
sum_data_long <- mutate(sum_data_long, canopy_cover = 1-sum_data_long$mean_gap)
sum_data_long <- merge(sum_data_long, invert_data, 'Sample_point', 'Sample_point', all.x=TRUE)
sum_data_long <- merge(sum_data_long, site_data, 'Sample_point', 'Sample_point', all.x=TRUE)
sum_data_long <- merge(sum_data_long, environ_data, by.x=c('Site', 'Start_date'), 
                       by.y=c('Site', 'Date'))
sum_data_long <- merge(sum_data_long, sound_data, 'Sample_point', 'Sample_point', all.x=TRUE)

# create factors
str(sum_data_long)
sum_data_long$Sample_point <- as.factor(sum_data_long$Sample_point)
sum_data_long$Site <- as.factor(sum_data_long$Site)
sum_data_long$Start_date <- as.factor(sum_data_long$Start_date)
sum_data_long$Species <- as.factor(sum_data_long$Species)
class(sum_data_long$Species)

require(lme4)
M1glmm <- glmer.nb(Passes~scale(lux) + scale(amp) + 
                     scale(lux)*scale(amp) +
                     scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                     scale(Mean_temp_night)+ 
                     (1|Species), 
                   data = sum_data_long)
summary(M1glmm)
require(performance)
model_performance(M1glmm)
# 0.78 R2 conditional 

# check fit of model
testDispersion(M1glmm)
testQuantiles(M1glmm)
testUniformity(M1glmm)
testZeroInflation(M1glmm)



### REPEAT FOR COMMON PIP 
# Use GLM with poisson
M1cp <- glm(com_pip_passes~scale(lux) + scale(amp) + 
              scale(lux)*scale(amp) +
              scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
              scale(Mean_temp_night), 
            data = emp_data, family = "poisson")
summary(M1cp)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(4585.3/5642.3)
# 19% of the variance in com pip passes can be explained 
sum(cooks.distance(M1cp)>1)
sum(cooks.distance(M1cp)>1)/length(emp_data$com_pip_passes)
# 8 outliers (4.5% of data)
# Check for over-dispersion: Dispersion parameter = resid dev / DF for resid dev 
4585.3/171
# The model is highly over-dispersed - 27x higher conditional variance to 
# conditional mean. 

## Try negative binomial distrib of errors
M2cp <- glm.nb(com_pip_passes~scale(lux) + scale(amp) + 
                 scale(lux)*scale(amp) +
                 scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                 scale(Mean_temp_night), 
               data = emp_data)
summary(M2cp)

### check the fit
# pseudoR2 = 1 - (residual dev / null dev) # 30% of the variance in bat passes can be explained 
1-(204.74/261.33)
# 22% of the variance in common pip passes can be explained 
sum(cooks.distance(M2cp)>1)
sum(cooks.distance(M2cp)>1)/length(emp_data$com_pip_passes)
#  no outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
204.74/171
# The model is slightly over-dispersed - 1.2x higher conditional variance to 
# conditional mean. 

# Use DHARMa package to test performance
testResiduals(M2cp)
testDispersion(M2cp)
testQuantiles(M2cp)
testUniformity(M2cp)
testZeroInflation(M2cp)
# no significant deviation of residuals, some quantile deviations but looks ok




### REPEAT FOR SOPRANO PIP 
# Use GLM with poisson
M1sp <- glm(sop_pip_passes~scale(lux) + scale(amp) + 
              scale(lux)*scale(amp) +
              scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
              scale(Mean_temp_night), 
            data = emp_data, family = "poisson")
summary(M1sp)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(4239.4/5759.7)
# 26% of the variance in sop pip passes can be explained

sum(cooks.distance(M1sp)>1)
sum(cooks.distance(M1sp)>1)/length(emp_data$sop_pip_passes)
# 7 outliers (4% of data)
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
4239.4/171
# The model is highly over-dispersed - 25x higher conditional variance to 
# conditional mean. 


## Try negative binomial distrib of errors
M2sp <- glm.nb(sop_pip_passes~scale(lux) + scale(amp) + 
                 scale(lux)*scale(amp) +
                 scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                 scale(Mean_temp_night), 
               data = emp_data)
summary(M2sp)

### check the fit
# pseudoR2 = 1 - (residual dev / null dev) # 30% of the variance in bat passes can be explained 
1-(208.04/273.59)
# 24% of the variance in sop pip passes can be explained
sum(cooks.distance(M2sp)>1)
#  no outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
208.04/171
# The model is slightly over-dispersed - 1.2x higher conditional variance to 
# conditional mean. 

# Use DHARMa package to test performance
testResiduals(M2sp)
testDispersion(M2sp)
testQuantiles(M2sp)
testUniformity(M2sp)
testZeroInflation(M2sp)
# some variance in residual vs predicted quantiles, but model looks ok




### REPEAT FOR Myotis
M1myot <- glm(myotis_passes~scale(lux) + scale(amp) + 
                scale(lux)*scale(amp) +
                scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                scale(Mean_temp_night), 
              data = emp_data, family = "poisson")
summary(M1myot)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(567.36/745.81)
# 24% of the variance in myotis passes can be explained 
sum(cooks.distance(M1myot)>1)
sum(cooks.distance(M1myot)>1)/length(emp_data$myotis_passes)
# no outliers 
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
567.36/171
# The model is over-dispersed - 3.3x higher conditional variance to conditional mean. 


## Try negative binomial distrib of errors
M2myot <- glm.nb(myotis_passes~scale(lux) + scale(amp) + 
                   scale(lux)*scale(amp) +
                   scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                   scale(Mean_temp_night), 
                 data = emp_data)
summary(M2myot)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev) # 30% of the variance in bat passes can be explained 
1-(138.24/184.35)
# 25% of the variance in myotis passes can be explained 
sum(cooks.distance(M2myot)>1)
sum(cooks.distance(M2myot)>1)/length(emp_data$myotis_passes)
# no outliers 
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
138.24/171
# The model is under-dispersed - 0.81x lower conditional variance to conditional mean. 

# Use DHARMa package to test performance
testResiduals(M2myot)
testDispersion(M2myot)
testQuantiles(M2myot)
testUniformity(M2myot)
testZeroInflation(M2myot)
# model looks ok



### REPEAT FOR Nyctalus and serotine grouping
M1nls <- glm(nls_passes~scale(lux) + scale(amp) + 
               scale(lux)*scale(amp) +
               scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
               scale(Mean_temp_night), 
             data = emp_data, family = "poisson")
summary(M1nls)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(326.03/356.62)
# only 8.5% of the variance in nls passes can be explained 
sum(cooks.distance(M1nls)>1)
sum(cooks.distance(M1nls)>1)/length(emp_data$nls_passes)
# one outlier
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
326.03/171
# The model is over-dispersed - 1.9x higher conditional variance to conditional mean. 


## Try negative binomial distrib of errors
M2nls <- glm.nb(nls_passes~scale(lux) + scale(amp) + 
                  scale(lux)*scale(amp) +
                  scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                  scale(Mean_temp_night), 
                data = emp_data)
summary(M2nls)
# pseudoR2 = 1 - (residual dev / null dev) only 9% of the variance in nls passes can be explained 
1-(156.74/171.62)
sum(cooks.distance(M2nls)>1)
sum(cooks.distance(M2nls)>1)/length(emp_data$nls_passes)
# one outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
156.74/171
# The model is slightly under-dispersed - 0.92x lower conditional variance to conditional mean. 




### INVESTIGATING TOTAL FEEDING RATE
# Use GLM with poisson and offset of bat passes (logged +1 to account for zeroes)
M1fb <- glm(total_bat_fb ~scale(lux) + scale(amp) + 
              scale(lux)*scale(amp) +
              scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
              scale(Mean_temp_night) + offset(log(total_bat_passes+1)), 
            data = emp_data, family = "poisson")
summary(M1fb)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(616.76/1497.57)
# 59% of the variance in bat feeding rate can be explained
sum(cooks.distance(M1fb)>1)
sum(cooks.distance(M1fb)>1)/length(emp_data$total_bat_fb)
# two outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
616.76/171
# The model is over-dispersed - 3.6x higher conditional variance to conditional mean. 

## Try using negative binomial to deal with overdispersion
M2fb <- glm.nb(total_bat_fb ~scale(lux) + scale(amp) + 
                 scale(lux)*scale(amp) +
                 scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                 scale(Mean_temp_night) + offset(log(total_bat_passes+1)), 
               data = emp_data)
summary(M2fb)

# pseudoR2 = 1 - (residual dev / null dev) 28% of the variance in bat feeding rate can be explained
1-(162.17/225.92)
sum(cooks.distance(M2fb)>1)
sum(cooks.distance(M2fb)>1)/length(emp_data$total_bat_fb)
# no outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
162.17/171
# The model is slightly under-dispersed - 0.95x lower conditional variance to conditional mean. 

# Use DHARMa package to test performance
testResiduals(M2fb)
testDispersion(M2fb)
testQuantiles(M2fb)
testUniformity(M2fb)
testZeroInflation(M2fb)
# predicted v observed residuals show significant deviations, but no dispersion or zero-inflation errors


## Repeat for common pip
# Use GLM with poisson and offset of bat passes (logged +1 to account for zeroes)
M1fbcp <- glm(com_pip_fb ~scale(lux) + scale(amp) + 
                scale(lux)*scale(amp) +
                scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                scale(Mean_temp_night) + offset(log(com_pip_passes+1)), 
              data = emp_data, family = "poisson")
summary(M1fbcp)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(358.48/855.48)
# 58% of the variance in bat feeding rate can be explained
sum(cooks.distance(M1fbcp)>1)
sum(cooks.distance(M1fbcp)>1)/length(emp_data$com_pip_fb)
# two outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
358.48/171
# The model is over-dispersed - 2.1x higher conditional variance to conditional mean. 

## Try using negative binomial to deal with overdispersion
M2fbcp <- glm.nb(com_pip_fb ~scale(lux) + scale(amp) + 
                   scale(lux)*scale(amp) +
                   scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                   scale(Mean_temp_night) + offset(log(com_pip_passes+1)), 
                 data = emp_data)
summary(M2fbcp)

### check the fit
# pseudoR2 = 1 - (residual dev / null dev) 27% of the variance in com pip feeding rate can be explained
1-(131.01/180.23)
sum(cooks.distance(M2fbcp)>1)
sum(cooks.distance(M2fbcp)>1)/length(emp_data$com_pip_fb)
# no outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
131.01/171
# The model is under-dispersed - 0.77x lower conditional variance to conditional mean. 

# Use DHARMa package to test performance
testResiduals(M2fbcp)
testDispersion(M2fbcp)
testQuantiles(M2fbcp)
testUniformity(M2fbcp)
testZeroInflation(M2fbcp)
# model looks ok, no significant issues detected



## Repeat for soprano pip
# Use GLM with poisson and offset of bat passes (logged +1 to account for zeroes)
M1fbsp <- glm(sop_pip_fb ~scale(lux) + scale(amp) + 
                scale(lux)*scale(amp) +
                scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                scale(Mean_temp_night) + offset(log(sop_pip_passes+1)), 
              data = emp_data, family = "poisson")
summary(M1fbsp)

### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(436.2/868.67)
# 49.7% of the variance in soprano pip feeding rate can be explained
sum(cooks.distance(M1fbsp)>1)
sum(cooks.distance(M1fbsp)>1)/length(emp_data$sop_pip_fb)
# one outlier
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
436.2/171
# The model is over-dispersed - 2.6x higher conditional variance to conditional mean. 

## Try using negative binomial to deal with overdispersion
M2fbsp <- glm.nb(sop_pip_fb ~scale(lux) + scale(amp) + 
                   scale(lux)*scale(amp) +
                   scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                   scale(Mean_temp_night) + offset(log(sop_pip_passes+1)), 
                 data = emp_data)
summary(M2fbsp)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev) 26% of the variance in sop pip feeding rate can be explained
1-(135.2/181.9)
sum(cooks.distance(M2fbsp)>1)
sum(cooks.distance(M2fbsp)>1)/length(emp_data$sop_pip_fb)
# no outliers
# Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
135.2/171
# The model is under-dispersed - 0.79x lower conditional variance to conditional mean. 
# Use DHARMa package to test performance
require(DHARMa)
testResiduals(M2fbsp)
par(mfrow=c(2,2))
testDispersion(M2fbsp)
testQuantiles(M2fbsp)
testUniformity(M2fbsp)
testZeroInflation(M2fbsp)
# model looks ok

exp(1.7614*2-.3884*2-1.0048*2)/exp(1.7614-.3884-1.0048)
exp(1.7614*2-.3884*2)/exp(1.7614-.3884)



### INVESTIGATING COMMUNICATION
# Use GLM with poisson for social calls
M1soc <- glm(social_calls ~scale(lux) + scale(amp) + 
               scale(lux)*scale(amp) +
               scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
               scale(Mean_temp_night) + offset(log(total_bat_passes+1)), 
             data = emp_data, family = "poisson")
summary(M1soc)


### check the fit
# pseudoR2 = 1 - (residual dev / null dev)
1-(6785.4/7524)
# 10% of the variance in bat social calls can be explained 
sum(cooks.distance(M1soc)>1)
sum(cooks.distance(M1soc)>1)/length(emp_data$social_calls)
# 15 outliers (8%)
## Check for overdispersion Dispersion parameter = resid dev / DF for resid dev 
6785.4/171
# The model is highly over-dispersed - 39.7x higher conditional variance to 
# conditional mean. 


## Try using negative binomial to deal with overdispersion
M2soc <- glm.nb(social_calls ~scale(lux) + scale(amp) + 
                  scale(lux)*scale(amp) +
                  scale(Distance_from_water_body) + scale(canopy_cover) + scale(Invert_count) + 
                  scale(Mean_temp_night) +offset(log(total_bat_passes+1)), 
                data = emp_data)
summary(M2soc)

### check the fit
# pseudoR2 = 1 - (residual dev / null dev) - 7.6% of the variance in bat social calls can be explained
1-(205.12/222.03)
testResiduals(M2soc)
testDispersion(M2soc)
testQuantiles(M2soc)
testUniformity(M2soc)
testZeroInflation(M2soc)
# deals with the dispersion issue, residuals v predicted signif diff in quantiles

sum(cooks.distance(M2soc)>1)
# no outliers

## Check for overdispersion: Dispersion parameter = resid dev / DF for resid dev 
205.12/171
# The model is slightly over-dispersed - 1.2x higher conditional variance to 
# conditional mean. 



# investigate relationship between feeding and social calls, to use as offset instead
lm_soc <- lm(social_calls ~ total_bat_fb, data = emp_data)
summary(lm_soc)
anova(lm_soc)

scatter_soc_feeding <- ggplot(data = emp_data, 
                          aes(x=total_bat_fb, y=social_calls)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
print(scatter_soc_feeding)
# conclude that no relationship so continue to use bat passes as offset

######################################################################################
######################################################################################
################################### plot interaction ##################################   
################################Feeding 
require(ggeffects)
require(ggplot2)
interaction_fb_1 <- plot(ggpredict(M2fbsp, terms = c("lux", "amp[2.0,3.2,3.9]"), ci = FALSE))+
  theme_classic() +
  xlab("Illuminance (lux)") +
  ylab("Predicted feeding attempts") +
  geom_line(linewidth = 2)+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 30), 
        legend.position = 'bottom', legend.title = element_text(size=30),
        legend.text=element_text(size=30),)
interaction_fb_1

interaction_fb_2<-plot(ggpredict(M2fbsp, terms = c("amp", "lux[0,0.07,0.34]"), ci = FALSE))+
  theme_classic() +
  xlab("Amplitude") +
  ylab("Predicted feeding attempts") +
  geom_line(linewidth = 2)+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 30), 
        legend.position = 'bottom', legend.title = element_text(size=30),
        legend.text=element_text(size=30),)
interaction_fb_2

summary(emp_data$lux)
summary(emp_data$amp)




########################################################################################
###################################### 2D heatmap plots #################################
# Code adapted from Flavia Bellotto-Trigo, June 2023


# prep data for heatmap
pred <- expand.grid(lux = seq(min(scale(emp_data$lux)), max(scale(emp_data$lux)), length.out = 100),
                    amp = seq(min(scale(emp_data$amp)), max(scale(emp_data$amp)), length.out = 100),
                    sop_pip_passes = seq(min(scale(emp_data$sop_pip_passes)), max(scale(emp_data$sop_pip_passes)), length.out = 100))

pred <- mutate(pred, Distance_from_water_body = 0, canopy_cover = 0,  Invert_count = 0, 
               Mean_temp_night = 0)


# use data from model
pred_fb <- mutate(pred, prediction = predict(M2fbsp, pred, type="response"))



# plot heatmap with contours
require(lattice)
require(metR)
require(viridis)
heatmap_fb <- ggplot(pred_fb,
                     aes(x = lux, y = amp, z= log(prediction))) +
  geom_tile(aes(fill = log(prediction))) + 
  geom_contour( colour = 'black') +
  scale_fill_viridis() +
  theme_classic() +
  xlab("Illuminance (lux; scaled)") +
  ylab("Amplitude (scaled)") +
  labs(fill='log(Feeding)') +
  scale_x_continuous(limits = c(0,4.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1,3), expand = c(0, 0)) +
  geom_text_contour(skip = 0, nudge_y = 0, nudge_x = 0, size = 7,
                    label.placer = label_placer_flattest(ref_angle = 65)) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=20),
        legend.text=element_text(size=17),)
heatmap_fb




################################## create table of stat outputs
require(stargazer)
stargazer(M1glmm, M2sp, M2cp, M2myot, M2nls, M2fbsp, M2fbcp, M2soc, 
          dep.var.labels = c("Overall passes", "Pippyg passes", "Pippip passes", "Myotis passes", "NycEpt passes", "Pippyg feeding", "Pippip feeding", "Social calls"),
          covariate.labels=c("Illuminance (lux)", "Amplitude","Distance from water", "Canopy Cover", "Invertebrate count", "Temperature", "Lux*Amplitude", "Intercept"),
          keep.stat = "n",
          star.cutoffs = c(0.05, 0.01, 0.001),
          out = "Model_stats.htm")


################################# Plot effect sizes ##############################
require(sjPlot)
require(ggtext)

M1glmm_plot <- plot_model(M1glmm, show.intercept = FALSE, show.values = FALSE, value.size=4, 
                    title = "", axis.title = "IRR for total bat passes")+
          theme_classic() +
          geom_hline(yintercept=1, linetype = "dotted")+
          theme(axis.title.x = ggtext::element_markdown())
  
M2cp_plot <- plot_model(M2cp, show.intercept = FALSE, show.values = FALSE, value.size=4, 
                        title = "", axis.title = "IRR for *P. pipistrellus* passes")+
          theme_classic() +
          geom_hline(yintercept=1, linetype = "dotted")+
          theme(axis.title.x = ggtext::element_markdown())

M2sp_plot <- plot_model(M2sp, show.intercept = FALSE, show.values = FALSE, value.size=4, 
                        title = "", axis.title = "IRR for *P. pygmaeus*  passes")+
           theme_classic() +
            geom_hline(yintercept=1, linetype = "dotted")+
          theme(axis.title.x = ggtext::element_markdown())

M2myot_plot <- plot_model(M2myot, show.intercept = FALSE, show.values = FALSE,  value.size=4, 
                          title = "", axis.title = "IRR for *Myotis* passes")+
           theme_classic() +
          geom_hline(yintercept=1, linetype = "dotted")+
          theme(axis.title.x = ggtext::element_markdown())

M2nls_plot <- plot_model(M2nls, show.intercept = FALSE, show.values = FALSE,  value.size=4, 
                         title = "", axis.title = "IRR for *Nyctalus*/*Eptesicus* passes")+
           theme_classic() +
           geom_hline(yintercept=1, linetype = "dotted")+
            theme(axis.title.x = ggtext::element_markdown())

M2fbsp_plot <- plot_model(M2fbsp, show.intercept = FALSE, show.values = FALSE, value.size=4, 
                          title = "", axis.title = "IRR for *P. pygmaeus* feeding") +
           theme_classic() +
           geom_hline(yintercept=1, linetype = "dotted")+
           theme(axis.title.x = ggtext::element_markdown())
M2fbsp_plot


M2fbcp_plot <- plot_model(M2fbcp, show.intercept = FALSE, show.values = FALSE, value.size=4, 
                          title = "", axis.title = "IRR for *P. pipistrellus* feeding")+
           theme_classic() +
           geom_hline(yintercept=1, linetype = "dotted")+
           theme(axis.title.x = ggtext::element_markdown())

M2soc_plot <- plot_model(M2soc, show.intercept = FALSE, show.values = FALSE, value.size=4, 
                         title = "", axis.title = "IRR for social calls")+
          theme_classic() +
          geom_hline(yintercept=1, linetype = "dotted")+
          theme(axis.title.x = ggtext::element_markdown())

passes_plot <- plot_models(M1glmm, M2cp, M2sp, M2myot, M2nls, 
                             axis.labels = c(
                             "Illuminance*Amplitude", "Temperature", "Invertebrate Count", "Canopy Cover", "Distance from water", "Amplitude", "Illuminance"),
                           m.labels = c("Overall passes", "Pippip passes", "Pippyg passes", "Myotis passes", "NLS passes"),
                           show.intercept = FALSE, show.values = FALSE, value.size=4, 
                           title = "", p.shape = TRUE)+
  theme_classic() +
  theme(plot.title=element_text(size=12, hjust=0.4))+
  ggtitle("Incidence Rate Ratios for fixed effects of bat activity models")+
  geom_hline(yintercept=1, linetype = "dotted", color = "red")
passes_plot

fb_plot <- plot_models(M2fbcp, M2fbsp, 
                           axis.labels = c(
                             "Illuminance*Amplitude", "Temperature", "Invertebrate Count", "Canopy Cover", "Distance from water", "Amplitude", "Illuminance"),
                           m.labels = c("Pippip feeding rate", "Pippyg feeding rate"),
                           show.intercept = FALSE, show.values = FALSE, value.size=4, 
                           title = "", p.shape = TRUE)+
  theme_classic() +
  theme(plot.title=element_text(size=12, hjust=0.4))+
  ggtitle("Incidence Rate Ratios for fixed effects of bat feeding models")+
  geom_hline(yintercept=1, linetype = "dotted", color = "red")
fb_plot

soc_plot <- plot_models(M2soc, 
                       axis.labels = c(
                         "Illuminance*Amplitude", "Temperature", "Invertebrate Count", "Canopy Cover", "Distance from water", "Amplitude", "Illuminance"),
                       m.labels = "Social call rate",
                       show.intercept = FALSE, show.values = FALSE, value.size=4, 
                       title = "", p.shape = TRUE)+
  theme_classic() +
  theme(plot.title=element_text(size=12, hjust=0.4))+
  ggtitle("Incidence Rate Ratios for fixed effects of bat communication model")+
  geom_hline(yintercept=1, linetype = "dotted", color = "red")
soc_plot







#######################################################################################
#######################################################################################
###########################    EXPERIMENTAL STUDY   ##################################
#######################################################################################
#######################################################################################




## RESPONSE VARIABLE

# Create column of bat passes, bat feeding buzzes and bat activity
bat_data <- read.csv("bat_data_manipulation_study.csv", header=TRUE)
bat_data <- mutate(bat_data, bat_passes = bat_data$com_pip_passes + 
                     bat_data$sop_pip_passes + bat_data$nat_pip_passes +
                     bat_data$noid_pip_passes + bat_data$nls_passes +
                     bat_data$myotis_passes + bat_data$unknown_passes) 
bat_data <- mutate(bat_data, bat_fb = bat_data$com_pip_fb + bat_data$sop_pip_fb +
                     bat_data$noid_pip_fb)
bat_data <- mutate(bat_data, bat_activity = if_else(bat_data$bat_passes==0,0,1))


bat_data$Transect <- as.factor(bat_data$Transect)
bat_data$Distance_from_treatment <- as.factor(bat_data$Distance_from_treatment)
bat_data$Treatment <- as.factor(bat_data$Treatment)

# extract start date by transect and by treatment
sum_data <- bat_data %>% 
  group_by(Distance_from_treatment, Transect, Treatment, .add = TRUE) %>%
  dplyr::select(com_pip_passes, sop_pip_passes, nat_pip_passes, noid_pip_passes, myotis_passes, nls_passes, 
         unknown_passes, bat_passes, com_pip_fb, sop_pip_fb, noid_pip_fb, unknown_fb, bat_fb, bat_activity, social_calls) %>% 
  summarize_all(sum)
sum_data
sum_data <- as.data.frame(sum_data)

start_data <- bat_data %>% 
  group_by(Distance_from_treatment, Transect, Treatment, .add = TRUE) %>%
  dplyr::select(Start_date) %>% 
  summarize_all(unique)
start_data
start_data <- as.data.frame(start_data)

location_data <- bat_data %>% 
  group_by(Distance_from_treatment, Transect, Treatment, .add = TRUE) %>%
  dplyr::select(Location) %>% 
  summarize_all(unique)
location_data
location_data <- as.data.frame(location_data)

mins_data <- bat_data %>% 
  group_by(Distance_from_treatment, Transect, Treatment, .add = TRUE) %>%
  dplyr::select(Time) %>% 
  summarize_all(length)
mins_data
mins_data <- as.data.frame(mins_data)

sum_data <- merge(sum_data, mins_data, by.x=c('Distance_from_treatment', 'Transect', 'Treatment'),
                   by.y=c('Distance_from_treatment', 'Transect', 'Treatment'))
colnames(sum_data)[colnames(sum_data) == 'Time'] <- 'Mins_recorded'

sum_data <- merge(sum_data, start_data, by.x=c('Distance_from_treatment', 'Transect', 'Treatment'),
                  by.y=c('Distance_from_treatment', 'Transect', 'Treatment'))
sum_data <- merge(sum_data, location_data, by.x=c('Distance_from_treatment', 'Transect', 'Treatment'),
                  by.y=c('Distance_from_treatment', 'Transect', 'Treatment'))


exp_data <- sum_data

## EXPLANATORY VARIABLES
# 1/ LUX and SPL
lux_data <- read.csv("light_sound_measurement_data.csv", header=TRUE)
#filter rows for point 2
lux_data <- filter(lux_data, Distance_from_treatment == 5 | Distance_from_treatment == 100 | Distance_from_treatment == 250)
# create average lux from two recordings for vetical hold
lux_data <- mutate (lux_data, lux = rowMeans(cbind(lux_data$Lux_vert1,lux_data$Lux_vert2)))
# select columns 
lux_data <- dplyr::select(lux_data, Distance_from_treatment, Transect, Treatment, lux, spl)
# merge with dataset
exp_data <- merge(exp_data, lux_data, by.x=c('Distance_from_treatment', 'Transect', 'Treatment'),
                  by.y=c('Distance_from_treatment', 'Transect', 'Treatment'))

## 6/ DAILY ENVIRONMENTAL DATA
environ_data <- read.csv("3hrly_environmental_data_manipulation_study.csv", header=TRUE)
# filter rows for night data
environ_data <- filter(environ_data, day_night == "night")

# extract mean temp
results <- data.frame(Start_date=character(), Air_temp=numeric())
results
for(i in unique(environ_data$Start_date)) {
  temp_file <- environ_data %>% subset(environ_data$Start_date==i)
  mean <- mean(temp_file$Air_temp)
  temp <- data.frame(Start_date=i, Air_temp=mean)
  results <- rbind(results, temp)
}
daily_data <- results

# extract mean wind
results <- data.frame(Start_date=character(), Wind_speed=numeric())
results
for(i in unique(environ_data$Start_date)) {
  temp_file <- environ_data %>% subset(environ_data$Start_date==i)
  mean <- mean(temp_file$Wind_speed)
  temp <- data.frame(Start_date=i, Wind_speed=mean)
  results <- rbind(results, temp)
}
daily_data <- merge(daily_data, results,'Start_date', 'Start_date', all.x=TRUE)

# select relevent dates
daily_data <- filter(daily_data, Start_date == "01/06/2023" | Start_date == "02/06/2023"  | 
                       Start_date == "07/06/2023"  | Start_date == "18/05/2023" | 
                       Start_date == "19/05/2023"  | Start_date == "21/05/2023"  | 
                       Start_date == "22/05/2023"  | Start_date == "23/05/2023"  | 
                       Start_date == "24/05/2023"  | Start_date == "25/05/2023"  | 
                       Start_date == "30/05/2023"  | Start_date == "31/05/2023"   )
# merge with dataset
exp_data <- merge(daily_data, exp_data,'Start_date', 'Start_date', all.x=TRUE)

str(exp_data)
exp_data$Start_date <- as.factor(exp_data$Start_date)
exp_data$Location <- as.factor(exp_data$Location)

# inspect the data
sum(exp_data$Mins_rec)
sum(exp_data$bat_activity)
sum(exp_data$bat_passes)
sum(exp_data$bat_fb)
sum(exp_data$social_calls)
sum(exp_data$com_pip_passes)
sum(exp_data$sop_pip_passes)
sum(exp_data$nls_passes)
sum(exp_data$myotis_passes)

sum(exp_data$nat_pip_passes+exp_data$noid_pip_passes+exp_data$unknown_passes)


##################################### light and noise data ######################
# LUX and SPL
lux_data <- read.csv("light_sound_measurement_data.csv", header=TRUE)
#filter rows for light and sound treatment
#lux_data <- filter(lux_data, Distance_from_treatment == 250  | Distance_from_treatment == 5)
# create average lux from two recordings for vetical hold
lux_data <- mutate (lux_data, lux = rowMeans(cbind(lux_data$Lux_vert1,lux_data$Lux_vert2)))
# select columns 
lux_data <- dplyr::select(lux_data, Distance_from_treatment, Location, Transect, Treatment, lux, spl)

lux_spl_data <- lux_data %>% 
  group_by(Distance_from_treatment, Treatment, .add = TRUE) %>%
  dplyr::select(lux, spl) %>% 
  summarize_all(sum)
lux_spl_data
lux_spl_data <- as.data.frame(lux_spl_data)
lux_spl_data <- mutate(lux_spl_data, mean_lux = lux_spl_data$lux/3)
lux_spl_data <- mutate(lux_spl_data, mean_spl = lux_spl_data$spl/3)


lux_spl_var_data <- lux_data %>% 
  group_by(Distance_from_treatment, Treatment, .add = TRUE) %>%
  dplyr::select(lux, spl) %>% 
  summarize_all(var)
lux_spl_var_data
lux_spl_var_data <- as.data.frame(lux_spl_var_data)
colnames(lux_spl_var_data)[colnames(lux_spl_var_data) == 'lux'] <- 'var_lux'
colnames(lux_spl_var_data)[colnames(lux_spl_var_data) == 'spl'] <- 'var_spl'
lux_spl_data <- merge(lux_spl_data, lux_spl_var_data, by.x=c('Distance_from_treatment', 'Treatment'),
                      by.y=c('Distance_from_treatment', 'Treatment'))
lux_spl_data <- mutate(lux_spl_data, se_lux = sqrt(lux_spl_data$var_lux/length(unique(bat_data$Transect))))
lux_spl_data <- mutate(lux_spl_data, se_spl = sqrt(lux_spl_data$var_spl/length(unique(bat_data$Transect))))

lux_spl_data$Distance_from_treatment <- as.factor(lux_spl_data$Distance_from_treatment)

##################################### Plot data ###################################
require(ggplot2)
require(viridis)

box_lux <- ggplot(lux_spl_data, aes(fill=Distance_from_treatment, y=mean_lux, x=Treatment)) +
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, alpha=.75)   +
  geom_errorbar(aes(ymin = mean_lux - se_lux, ymax = mean_lux +se_lux), position = "dodge") +
  xlab("Treatment") +
  ylab("Mean illuminance") +
  labs(fill='Distance from treatment site (m)') +
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 30), legend.position="none")
box_lux

box_spl <- ggplot(lux_spl_data, aes(fill=Distance_from_treatment, y=mean_spl, x=Treatment)) +
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, alpha=.75)   +
  geom_errorbar(aes(ymin = mean_spl - se_spl, ymax = mean_spl +se_spl), position = "dodge") +
  xlab("Treatment") +
  ylab("Mean Sound Pressure Level") +
  labs(fill='Distance (m)') +
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 30), 
        legend.title = element_text(size=30),legend.text=element_text(size=30), legend.position = 'bottom')
box_spl


######################### Modelling the data##############################################


### INVESTIGATING TOTAL BAT PASSES
# Use GLM with poisson
# filter out light treatment
exp_data_light <- filter(exp_data, Treatment == "Light"  | Treatment == "Control")
M1l <- glm(bat_passes~Treatment + Location + Location*Treatment, 
         data = exp_data_light, family="poisson")
summary(M1l)
par(mfrow=c(2,2))
testDispersion(M1l)
testQuantiles(M1l)
testUniformity(M1l)
testZeroInflation(M1l)

# filter out sound treatment
exp_data_sound <- filter(exp_data, Treatment == "Sound"  | Treatment == "Control")
M1s <- glm(bat_passes~Treatment + Location + Location*Treatment, 
          data = exp_data_sound, family="poisson")
summary(M1s)
testDispersion(M1s)
testQuantiles(M1s)
testUniformity(M1s)
testZeroInflation(M1s)
# filter out light_sound treatment
exp_data_light_sound <- filter(exp_data, Treatment == "Light_Sound"  | Treatment == "Control")
M1ls <- glm(bat_passes~Treatment + Location + Location*Treatment, 
           data = exp_data_light_sound, family="poisson")
summary(M1ls)
testDispersion(M1ls)
testQuantiles(M1ls)
testUniformity(M1ls)
testZeroInflation(M1ls)


# Use GLM with poisson - check if difference using start_date rather than treatment
# filter out light treatment
exp_data_light <- filter(exp_data, Treatment == "Light"  | Treatment == "Control")
M1l <- glm(bat_passes~Start_date + Location + Location*Start_date, 
           data = exp_data_light, family="poisson")
summary(M1l)
# filter out sound treatment
exp_data_sound <- filter(exp_data, Treatment == "Sound"  | Treatment == "Control")
M1s <- glm(bat_passes~Start_date + Location + Location*Start_date, 
           data = exp_data_sound, family="poisson")
summary(M1s)
# filter out light_sound treatment
exp_data_light_sound <- filter(exp_data, Treatment == "Light_Sound"  | Treatment == "Control")
M1ls <- glm(bat_passes~Start_date + Location + Location*Start_date, 
            data = exp_data_light_sound, family="poisson")
summary(M1ls)



################################## box plot by location ##################################
mean_data <- exp_data %>% 
  group_by(Treatment, Location, .add = TRUE) %>%
  dplyr::select(bat_passes) %>% 
  summarize_all(sum)
mean_data
mean_data <- as.data.frame(mean_data)
mean_data <- mutate(mean_data, mean_passes = mean_data$bat_passes/3)

var_data <- exp_data %>% 
  group_by(Treatment, Location, .add = TRUE) %>%
  dplyr::select(bat_passes) %>% 
  summarize_all(var)
var_data
var_data <- as.data.frame(var_data)
colnames(var_data)[colnames(var_data) == 'bat_passes'] <- 'var'
mean_data <- merge(mean_data, var_data,by.x=c('Location', 'Treatment'),
                   by.y=c('Location', 'Treatment'))
mean_data <- mutate(mean_data, se = sqrt(mean_data$var/length(unique(bat_data$Transect))))


box_passes <- ggplot(mean_data, aes(fill=Location, y=mean_passes, x=Treatment)) +
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, alpha=.75)   +
  geom_errorbar(aes(ymin = mean_passes - se, ymax = mean_passes +se), position = "dodge") +
  xlab("Treatment") +
  ylab("Mean bat passes") +
  labs(fill='Location') +
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20), 
        legend.position = 'bottom', legend.title = element_text(size=17),legend.text=element_text(size=17))
box_passes

################################ checking light and sound v control ########################

lux_100 <- filter(lux_data, Distance_from_treatment == 100)
M2 <- lm(lux~Treatment, data = lux_100)
summary(M2)
anova(M2)

lux_50 <- filter(lux_data, Distance_from_treatment == 50)
M3 <- lm(lux~Treatment, data = lux_50)
summary(M3)
anova(M3)

lux_20 <- filter(lux_data, Distance_from_treatment == 20)
M4 <- lm(lux~Treatment, data = lux_20)
summary(M4)
anova(M4)

lux_5 <- filter(lux_data, Distance_from_treatment == 5)
M5 <- lm(lux~Treatment, data = lux_5)
summary(M5)
anova(M5)

# spl

M6 <- lm(spl~Treatment, data = lux_100)
summary(M6)
anova(M6)






#################################################################################
#################################################################################
#######################   AUTO ANALYSIS    ######################################
#################################################################################
#################################################################################
list.files()
library(tidyverse)
library(viridis) 


bat_data <- read.csv("bat_data_auto_analysis_raw.csv", header=TRUE)

# summarise kaleidescope data into containing bats or not
bat_data <- mutate_all(bat_data, ~ifelse(is.na(.), 0, .))
bat_data <- mutate(bat_data, kal_sum = bat_data$Pippip + bat_data$Pippyg + bat_data$Myonat +
                     bat_data$Nycnoc + bat_data$Nyclas + bat_data$NoID + bat_data$Nyclei +
                     bat_data$Pleaus + bat_data$Pleaur + bat_data$Eptser + bat_data$Pipnat +
                     bat_data$Myotis +bat_data$Barbar) 
bat_data <- mutate(bat_data, kal_activity = if_else(bat_data$kal_sum==0,0,1))

# summarise manual data into containing bats or not
bat_data <- mutate(bat_data, manual_sum = bat_data$com_pip_passes + bat_data$sop_pip_passes +
                     bat_data$noid_pip_passes + bat_data$nat_pip_passes + bat_data$nls_passes +
                     bat_data$myotis_passes +bat_data$unknown_passes) 
bat_data <- mutate(bat_data, manual_activity = if_else(bat_data$manual_sum==0,0,1))

# check for false negatives
bat_data_activity <- filter(bat_data, bat_data$manual_activity==1)

bat_data_activity <- mutate(bat_data_activity, kal_false_neg = if_else(bat_data_activity$kal_activity==0,1,0))
sum(bat_data_activity$kal_false_neg)
sum(bat_data_activity$kal_false_neg)/sum(bat_data_activity$manual_activity)

bat_data_activity <- mutate(bat_data_activity, sonobat_false_neg = if_else(bat_data_activity$Sonobat_filtered_tolerant==0,1,0))
sum(bat_data_activity$sonobat_false_neg)
sum(bat_data_activity$sonobat_false_neg)/sum(bat_data_activity$manual_activity)

# check for false positives
bat_data_no_activity <- filter(bat_data, bat_data$manual_activity==0)

bat_data_no_activity <- mutate(bat_data_no_activity, kal_false_pos = if_else(bat_data_no_activity$kal_activity==1,1,0))
sum(bat_data_no_activity$kal_false_pos)
sum(bat_data_no_activity$kal_false_pos)/sum(bat_data$kal_activity)

bat_data_no_activity <- mutate(bat_data_no_activity, sonobat_false_pos = if_else(bat_data_no_activity$Sonobat_filtered_tolerant==1,1,0))
sum(bat_data_no_activity$sonobat_false_pos)
sum(bat_data_no_activity$sonobat_false_pos)/sum(bat_data$Sonobat_filtered_tolerant)



length(bat_data$Site)
sum(bat_data$Sonobat_filtered_tolerant)
sum(bat_data$Sonobat_filtered_tolerant)/length(bat_data$Site)
sum(bat_data$kal_activity)
sum(bat_data$kal_activity)/length(bat_data$Site)
sum(bat_data$manual_activity)
sum(bat_data$manual_activity)/length(bat_data$Site)


### AUTO-FILTERING RESULTS 
# manual true positive and negative
sum(bat_data$manual_activity)
length(bat_data$Site)-sum(bat_data$manual_activity)
# kaleidoscope true positive, false positive, true negative, false negative
sum(bat_data$kal_activity)
sum(bat_data$kal_activity)-sum(bat_data_no_activity$kal_false_pos)
sum(bat_data_no_activity$kal_false_pos)
length(bat_data$Site)-sum(bat_data$kal_activity)-sum(bat_data_activity$kal_false_neg)
sum(bat_data_activity$kal_false_neg)

# sonobat true positive, false positive, true negative, false negative
sum(bat_data$Sonobat_filtered_tolerant)
sum(bat_data$Sonobat_filtered_tolerant)-sum(bat_data_no_activity$sonobat_false_pos)
sum(bat_data_no_activity$sonobat_false_pos)
length(bat_data$Site)-sum(bat_data$Sonobat_filtered_tolerant)-sum(bat_data_activity$sonobat_false_neg)
sum(bat_data_activity$sonobat_false_neg)
sum(bat_data_no_activity$sonobat_false_pos)/length(bat_data$Site)

# create dataframe of results
df <- data.frame(Filter_method=rep(c("Manual", "Kaleidoscope", "Sonobat"), times=4),
                 Result=rep(c("D: True Positive", "B: False Positive", "A: True Negative", "C: False Negative"), each=3),
                 files=c((sum(bat_data$manual_activity)/length(bat_data$Site)),((sum(bat_data$kal_activity)-sum(bat_data_no_activity$kal_false_pos))/length(bat_data$Site)),  ((sum(bat_data$Sonobat_filtered_tolerant)-sum(bat_data_no_activity$sonobat_false_pos))/length(bat_data$Site)),
                         0, (sum(bat_data_no_activity$kal_false_pos)/length(bat_data$Site)),(sum(bat_data_no_activity$sonobat_false_pos)/length(bat_data$Site)), 
                         ((length(bat_data$Site)-sum(bat_data$manual_activity))/length(bat_data$Site)), ((length(bat_data$Site)-sum(bat_data$kal_activity)-sum(bat_data_activity$kal_false_neg))/length(bat_data$Site)),((length(bat_data$Site)-sum(bat_data$Sonobat_filtered_tolerant)-sum(bat_data_activity$sonobat_false_neg))/length(bat_data$Site)),
                         0,(sum(bat_data_activity$kal_false_neg)/length(bat_data$Site)),  (sum(bat_data_activity$sonobat_false_neg))/length(bat_data$Site)))
df

df <- data.frame(Filter_method=rep(c("Manual", "Kaleidoscope", "Sonobat"), times=4),
                 Result=rep(c("D: True Positive", "B: False Positive", "A: True Negative", "C: False Negative"), each=3),
                 files=c(8.6,6.5, 7.6,
                         0, 12.0,68.4, 
                         91.3, 79.4,23.0,
                         0,2.1,1.0)) 
df

ggplot(df, aes(fill=Result, y=files, x=Filter_method)) + 
  geom_bar(position='stack', stat='identity', width = 0.75)+
  theme_minimal() + 
  labs(x='Method of bat file filtering', y='Files') +
  theme(plot.title = element_text(hjust=0.5, size=19), axis.text=element_text(size=19),
        axis.title=element_text(size=19), legend.title = element_blank (),
        legend.text=element_text(size=17),legend.position="bottom") +
  scale_fill_viridis(discrete = TRUE, alpha=0.8, direction =-1)



### IMPACT ON SPECIES BREAKDOWN
## Manual breakdown
man_cp <- sum(bat_data_activity$com_pip_passes)/sum(bat_data_activity$manual_sum)
man_sp <- sum(bat_data_activity$sop_pip_passes)/sum(bat_data_activity$manual_sum)
man_nls <- sum(bat_data_activity$nls_passes)/sum(bat_data_activity$manual_sum)
man_myot <- sum(bat_data_activity$myotis_passes)/sum(bat_data_activity$manual_sum)
man_other <- (sum(bat_data_activity$unknown_passes)+sum(bat_data_activity$noid_pip_passes)+sum(bat_data_activity$nat_pip_passes))/sum(bat_data_activity$manual_sum)
sum(man_cp,man_sp,man_other,man_nls,man_myot)


## Breakdown using Kaleidoscope auto-filter and manual ID
# filter results by kaleidoscope activity
bat_data_kal_activity <- filter(bat_data, bat_data$kal_activity==1)
# filter for manual activity - ie the true positives
bat_data_kal_truepos <- filter(bat_data_kal_activity, bat_data_kal_activity$manual_activity==1)
# calculate breakdown
kal_cp <- sum(bat_data_kal_truepos$com_pip_passes)/sum(bat_data_kal_truepos$manual_sum)
kal_sp <- sum(bat_data_kal_truepos$sop_pip_passes)/sum(bat_data_kal_truepos$manual_sum)
kal_nls <- sum(bat_data_kal_truepos$nls_passes)/sum(bat_data_kal_truepos$manual_sum)
kal_myot <- sum(bat_data_kal_truepos$myotis_passes)/sum(bat_data_kal_truepos$manual_sum)
kal_other <- (sum(bat_data_kal_truepos$unknown_passes)+sum(bat_data_kal_truepos$noid_pip_passes)+sum(bat_data_kal_truepos$nat_pip_passes))/sum(bat_data_kal_truepos$manual_sum)
sum(kal_cp,kal_sp,kal_nls,kal_myot,kal_other)

# calculate % of true passes found
sum(bat_data_kal_truepos$com_pip_passes)/sum(bat_data_activity$com_pip_passes)
sum(bat_data_kal_truepos$sop_pip_passes)/sum(bat_data_activity$sop_pip_passes)
sum(bat_data_kal_truepos$nls_passes)/sum(bat_data_activity$nls_passes)
sum(bat_data_kal_truepos$myotis_passes)/sum(bat_data_activity$myotis_passes)


## Breakdown using Sonobat auto-filter and manual ID
# filter results by sonobat activity
bat_data_sono_activity <- filter(bat_data, bat_data$Sonobat_filtered_tolerant==1)
# filter for manual activity - ie the true positives
bat_data_sono_truepos <- filter(bat_data_sono_activity, bat_data_sono_activity$manual_activity==1)
# calculate breakdown
sono_cp <- sum(bat_data_sono_truepos$com_pip_passes)/sum(bat_data_sono_truepos$manual_sum)
sono_sp <- sum(bat_data_sono_truepos$sop_pip_passes)/sum(bat_data_sono_truepos$manual_sum)
sono_nls <- sum(bat_data_sono_truepos$nls_passes)/sum(bat_data_sono_truepos$manual_sum)
sono_myot <- sum(bat_data_sono_truepos$myotis_passes)/sum(bat_data_sono_truepos$manual_sum)
sono_other <- (sum(bat_data_sono_truepos$unknown_passes)+sum(bat_data_sono_truepos$noid_pip_passes)+sum(bat_data_sono_truepos$nat_pip_passes))/sum(bat_data_sono_truepos$manual_sum)
sum(sono_cp,sono_sp,sono_nls,sono_myot,sono_other)

# calculate % of true passes found
sum(bat_data_sono_truepos$com_pip_passes)/sum(bat_data_activity$com_pip_passes)
sum(bat_data_sono_truepos$sop_pip_passes)/sum(bat_data_activity$sop_pip_passes)
sum(bat_data_sono_truepos$nls_passes)/sum(bat_data_activity$nls_passes)
sum(bat_data_sono_truepos$myotis_passes)/sum(bat_data_activity$myotis_passes)


## Breakdown using auto-filter and auto-ID
bat_data_kal_activity <- mutate(bat_data_kal_activity, Myotis_all = bat_data_kal_activity$Myonat + bat_data_kal_activity$Myotis)
bat_data_kal_activity <- mutate(bat_data_kal_activity, NLS_all = bat_data_kal_activity$Nycnoc + bat_data_kal_activity$Nyclas + bat_data_kal_activity$Nyclei + bat_data_kal_activity$Eptser)
bat_data_kal_activity <- mutate(bat_data_kal_activity, Other_all = bat_data_kal_activity$Pipnat + bat_data_kal_activity$Pleaus + bat_data_kal_activity$Pleaur + bat_data_kal_activity$Barbar+ bat_data_kal_activity$NoID)

bat_data_kal_activity <- mutate(bat_data_kal_activity, Pippip_act = if_else(bat_data_kal_activity$Pippip==0,0,1))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Pippyg_act = if_else(bat_data_kal_activity$Pippyg==0,0,1))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Myotis_act = if_else(bat_data_kal_activity$Myotis_all==0,0,1))
bat_data_kal_activity <- mutate(bat_data_kal_activity, NLS_act = if_else(bat_data_kal_activity$NLS_all ==0,0,1))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Other_act = if_else(bat_data_kal_activity$Other_all ==0,0,1))

bat_data_kal_activity <- mutate(bat_data_kal_activity, kal_passes_sum = bat_data_kal_activity$Pippip_act +bat_data_kal_activity$Pippyg_act + bat_data_kal_activity$Myotis_act + bat_data_kal_activity$NLS_act + bat_data_kal_activity$Other_act)

# calculate breakdown
kal_auto_cp <- sum(bat_data_kal_activity$Pippip_act)/sum(bat_data_kal_activity$kal_passes_sum)
kal_auto_sp <- sum(bat_data_kal_activity$Pippyg_act)/sum(bat_data_kal_activity$kal_passes_sum)
kal_auto_other <- sum(bat_data_kal_activity$Other_act)/sum(bat_data_kal_activity$kal_passes_sum)
kal_auto_nls <- sum(bat_data_kal_activity$NLS_act)/sum(bat_data_kal_activity$kal_passes_sum)
kal_auto_myot <- sum(bat_data_kal_activity$Myotis_act)/sum(bat_data_kal_activity$kal_passes_sum)
sum(kal_auto_cp,kal_auto_sp,kal_auto_nls,kal_auto_myot,kal_auto_other)


## Breakdown using auto-filter and auto-ID over set level
set_level <- 0.85
bat_data_kal_activity <- mutate(bat_data_kal_activity, Pippip_act_70 = if_else(bat_data_kal_activity$Pippip>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Pippyg_act_70 = if_else(bat_data_kal_activity$Pippyg>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Myonat_act_70 = if_else(bat_data_kal_activity$Myonat>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Myotis_act_70 = if_else(bat_data_kal_activity$Myotis>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Nycnoc_act_70 = if_else(bat_data_kal_activity$Nycnoc>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Nyclas_act_70 = if_else(bat_data_kal_activity$Nyclas>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Nyclei_act_70 = if_else(bat_data_kal_activity$Nyclei>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Eptser_act_70 = if_else(bat_data_kal_activity$Eptser>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Pipnat_act_70 = if_else(bat_data_kal_activity$Pipnat>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Pleaus_act_70 = if_else(bat_data_kal_activity$Pleaus>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Pleaur_act_70 = if_else(bat_data_kal_activity$Pleaur>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, Barbar_act_70 = if_else(bat_data_kal_activity$Barbar>=set_level,1,0))
bat_data_kal_activity <- mutate(bat_data_kal_activity, NoID_act_70 = if_else(bat_data_kal_activity$NoID>=set_level,1,0))

bat_data_kal_activity <- mutate(bat_data_kal_activity, Myotis_all_sum_70 = bat_data_kal_activity$Myonat_act_70 +bat_data_kal_activity$Myotis_act_70)
bat_data_kal_activity <- mutate(bat_data_kal_activity, Myotis_all_act_70 = if_else(bat_data_kal_activity$Myotis_all_sum_70>0,1,0))

bat_data_kal_activity <- mutate(bat_data_kal_activity, NLS_all_sum_70 = bat_data_kal_activity$Nycnoc_act_70 +bat_data_kal_activity$Nyclas_act_70 +bat_data_kal_activity$Nyclei_act_70 + bat_data_kal_activity$Eptser_act_70)
bat_data_kal_activity <- mutate(bat_data_kal_activity, NLS_all_act_70 = if_else(bat_data_kal_activity$NLS_all_sum_70>0,1,0))

bat_data_kal_activity <- mutate(bat_data_kal_activity, Other_all_sum_70 = bat_data_kal_activity$Pipnat_act_70 +bat_data_kal_activity$Pleaus_act_70 +bat_data_kal_activity$Pleaur_act_70 + bat_data_kal_activity$Barbar_act_70 + bat_data_kal_activity$NoID_act_70)
bat_data_kal_activity <- mutate(bat_data_kal_activity, Other_all_act_70 = if_else(bat_data_kal_activity$Other_all_sum_70>0,1,0))

bat_data_kal_activity <- mutate(bat_data_kal_activity, kal_passes_sum_70 = bat_data_kal_activity$Pippip_act_70 +bat_data_kal_activity$Pippyg_act_70 + bat_data_kal_activity$Myotis_all_act_70 + bat_data_kal_activity$NLS_all_act_70 + bat_data_kal_activity$Other_all_act_70)

# calculate breakdown
kal_auto_70_cp <- sum(bat_data_kal_activity$Pippip_act_70)/sum(bat_data_kal_activity$kal_passes_sum_70)
kal_auto_70_sp <- sum(bat_data_kal_activity$Pippyg_act_70)/sum(bat_data_kal_activity$kal_passes_sum_70)
kal_auto_70_other <- sum(bat_data_kal_activity$Other_all_act_70)/sum(bat_data_kal_activity$kal_passes_sum_70)
kal_auto_70_nls <- sum(bat_data_kal_activity$NLS_all_act_70)/sum(bat_data_kal_activity$kal_passes_sum_70)
kal_auto_70_myot <- sum(bat_data_kal_activity$Myotis_all_act_70)/sum(bat_data_kal_activity$kal_passes_sum_70)
sum(kal_auto_70_cp,kal_auto_70_sp,kal_auto_70_nls,kal_auto_70_myot,kal_auto_70_other)



#### CREATE CHART

Species <- c('P. pipistrellus','P. pygmaeus', 'Nyctalus/Eptesicus', 'Myotis', 'Other')
Manual <- c(man_cp,man_sp,man_nls,man_myot,man_other)
Kaleidoscope <- c( kal_cp,kal_sp,kal_nls,kal_myot,kal_other)
Kaleidoscope_autoID <-  c( kal_auto_cp,kal_auto_sp,kal_auto_nls,kal_auto_myot,kal_auto_other)
Kaleidoscope_autoID_over70 <-  c(kal_auto_70_cp,kal_auto_70_sp,kal_auto_70_nls,kal_auto_70_myot,kal_auto_70_other)
df_sp <- data.frame(Species, Manual, Kaleidoscope, Kaleidoscope_autoID, Kaleidoscope_autoID_over70)
df_sp
library(ggplot2)

pie_man <- ggplot(df_sp, aes(x="", y=Manual, fill=Species)) + geom_bar(stat="identity", width=1) + 
  coord_polar("y", start=0) + geom_text(aes(x =  1.3,label = paste0(round(Manual*100), "%")), position = position_stack(vjust = 0.5)) + 
  scale_fill_viridis(discrete = TRUE, alpha=.75) +  
  labs(x = NULL, y = NULL, fill = NULL, title = "Manual") + 
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          plot.title = element_text(hjust = 0.5, size=20)) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.text=element_text(face="italic")) +
  theme(legend.position="bottom") #+
#theme(legend.position="none")

pie_kal <- ggplot(df_sp, aes(x="", y=Kaleidoscope, fill=Species)) + geom_bar(stat="identity", width=1) + 
  coord_polar("y", start=0) + geom_text(aes(x =  1.3,label = paste0(round(Kaleidoscope*100), "%")), position = position_stack(vjust = 0.5)) + 
  scale_fill_viridis(discrete = TRUE, alpha=.75) +  
  labs(x = NULL, y = NULL, fill = NULL, title = "Kaleidoscope") + 
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          plot.title = element_text(hjust = 0.5, size=20)) +
  theme(legend.position="none")

pie_kal_autoID <- ggplot(df_sp, aes(x="", y=Kaleidoscope_autoID, fill=Species)) + geom_bar(stat="identity", width=1) + 
  coord_polar("y", start=0) + geom_text(aes(x =  1.3,label = paste0(round(Kaleidoscope_autoID*100), "%")), position = position_stack(vjust = 0.5)) + 
  scale_fill_viridis(discrete = TRUE, alpha=.75) +  
  labs(x = NULL, y = NULL, fill = NULL, title = "Kaleidoscope with auto-ID") + 
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          plot.title = element_text(hjust = 0.5, size=20)) +
  theme(legend.position="none")

pie_kal_autoID_over70 <- ggplot(df_sp, aes(x="", y=Kaleidoscope_autoID_over70, fill=Species)) + geom_bar(stat="identity", width=1) + 
  coord_polar("y", start=0) + geom_text(aes(x =  1.3,label = paste0(round(Kaleidoscope_autoID_over70*100), "%")), position = position_stack(vjust = 0.5)) + 
  scale_fill_viridis(discrete = TRUE, alpha=.75) +  
  labs(x = NULL, y = NULL, fill = NULL, title = "Kaleidoscope with auto-ID over 0.85") + 
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          plot.title = element_text(hjust = 0.5, size=20)) +
  theme(legend.position="none")


######  CREATING MULTIPLOT FUNCTION
###### taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Plot all four with one legend
multiplot(pie_man, pie_kal, cols=2)
multiplot(box_lux, box_spl)
multiplot(M1glmm_plot, M2cp_plot, M2sp_plot, M2myot_plot, cols=2)
multiplot(M2nls_plot, M2fbsp_plot, M2fbcp_plot, M2soc_plot, cols=2)
multiplot(passes_plot, fb_plot, soc_plot)
