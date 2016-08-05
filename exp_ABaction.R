##################################################
# Seaching for fingerprints of antibiotic action
##################################################

# Author: Rene Niehus
# Date: 1 August 2016

# clear workspace, load libraries
rm(list=ls())
library(dplyr)
library(reshape)

# set the working directory
setwd("~/Dropbox/LOMHWRU_MORU/SATURN_ESBL/_R/R_git/qPCR/")

### My functions
# funcion that compares neighbouring elements of a data vector
Comp <- function(data)
{
  output <- vector()
  for(i in 1:(length(data) - 1))
  {
    if(data[i] == data[i + 1])
    {
      output[i] <- 1
    }
    else
    {
      output[i] <- 0
    }
  }
  return(output)
}

# Load Data file
SDATA <- read.csv ("./Cleaned_data/linked_qPCR_clin_abx.csv",
                   sep= ",", colClasses=c("character")) # Linked Data

# turn into numerical 
SDATA <- SDATA %>% 
  plyr::mutate(
         ratio_SD = as.numeric(ratio_SD),
         s_num = as.numeric(s_num),
         qu_ratio = as.numeric(qu_ratio),
         esbl_act = as.numeric(esbl_act),
         broad_spec = as.numeric(broad_spec),
         Piperacillin.tazobactam = as.numeric(Piperacillin.tazobactam),
         Vancomycin = as.numeric(Vancomycin),
         isoniazid = as.numeric(isoniazid),
         Levofloxacin = as.numeric(Levofloxacin),
         Meropenem = as.numeric(Meropenem),
         pyrazinamide = as.numeric(pyrazinamide),
         Rifampicin = as.numeric(Rifampicin),
         Ampicillin = as.numeric(Ampicillin),
         Oxacillin = as.numeric(Oxacillin),
         Ceftriaxone = as.numeric(Ceftriaxone),
         ESBL_Ecoli = as.numeric(ESBL_Ecoli),
         ESBL_PPM = as.numeric(ESBL_PPM),
         ESBL_KESC = as.numeric(ESBL_KESC)
  )

# turn into dates
SDATA <- SDATA %>%
  plyr::mutate(RectalDate = as.Date(as.character(SDATA$RectalDate), format="%Y-%m-%d"))

# order after rectal data FOR EACH patient
SDATA <- SDATA[order(SDATA$patient_id, SDATA$RectalDate),]

######## Some Descriptive analysis
dim(SDATA)[1:2] # 659 observations, 127 variables
head(SDATA, n = -657) # negative n means the length of cut-off tail
names(SDATA) # gives a nice list of a the variables
str(SDATA) # gives a summary like in the global environemnt
levels(SDATA$patient_id) # this is only for categorical data
summary(SDATA$qu_ratio) # gives a summery of a coloumn
# the median is 0.0077: majority of measures are close to zero. I expect that this ratio is mostly very small. 
# how accurate are small values of the ratio. Look at the error
summary(SDATA$ratio_SD) # gives a summery of a coloumn
# the SD of the ratio is around 1% of 
mean.quratio = mean(SDATA$qu_ratio, na.rm = T) # calculate the mean
max.quratio = max(SDATA$qu_ratio, na.rm = T) # calculate the mean
# make a histogram and density plot
png('./Figures/histogram and kernel density plot.png')
hist(SDATA$qu_ratio, breaks = 5, freq = F, xlab = 'CTXm-16s ratio', xlim = c(0,20), ylim = c(0, 1), ylab = 'Probability', main = 'Histogram of CTXm abundance with Kernel Density Plot')
lines(density(SDATA$qu_ratio, na.rm = T, from = 0, to = max.quratio))
dev.off()

########
# small table version of SDATA
#names(SDATA)[90:126] # those are all the antibiotics in action
DF1 <- SDATA %>%
  dplyr::select(num,qu_ratio,Tdiff,RectalDate,esbl_act,broad_spec,
                Piperacillin.tazobactam,Vancomycin,isoniazid,Levofloxacin,
                Meropenem,pyrazinamide,Rifampicin,Ampicillin,Oxacillin,
                Ceftriaxone,Ceftazidime)

# replace all NA by 0
DF1[is.na(DF1)] <- 0

# remove the qu_ratio over 100
DF1 <- DF1[-which(max(DF1$qu_ratio) == DF1$qu_ratio),]
DF1 <- DF1[-which(max(DF1$qu_ratio) == DF1$qu_ratio),]

# add the time since last measurement as predictor
# add a column with previous date
next.date <- c(DF1$RectalDate[2:(length(DF1$RectalDate))] , as.Date("2011-1-1"))
DF1$TdiffToNext <- next.date - DF1$RectalDate

# the next qu_ratio is what i want to predict
DF1$NextQuRatio = c(DF1$qu_ratio[2 : (nrow(DF1))], 0 )

# 
first.meas = which(DF1$num == 0)
first.meas <- first.meas[2 : length(first.meas)] - 1

# Filter away first measurements
DF2 <- DF1[-first.meas,]

# Make new DF for regression
DF3 <- DF2 %>% dplyr::select(num,qu_ratio,esbl_act,broad_spec,
                             Piperacillin.tazobactam,Vancomycin,isoniazid,Levofloxacin,
                             Meropenem,pyrazinamide,Rifampicin,Ampicillin,Oxacillin,
                             Ceftriaxone,Ceftazidime,
                             TdiffToNext,NextQuRatio)

### Now perform different regressions
# Output: NextQuRatio -- Predicting Features: qu_ratio, esbl_act, broad_spec
DF_R1 <-  DF3 %>% dplyr::select(num,qu_ratio,broad_spec,NextQuRatio)

# Fit the NextQuRatio using the previous and broad_spec
fit.1 <- lm(DF_R1$NextQuRatio ~ DF_R1$broad_spec + DF_R1$qu_ratio)
fit.1
# graphical display
colors <- ifelse(DF_R1$broad_spec==1, "black", "gray")
plot(DF_R1$qu_ratio, DF_R1$NextQuRatio, xlab="Previous QuRatio", ylab="Next Qu",
     col=colors, pch=20)
curve(cbind(1, 1, x) %*% coef(fit.1), add=TRUE, col="black")
curve(cbind(1, 0,x) %*% coef(fit.1), add=TRUE, col="gray")

# Output: NextQuRatio -- Predicting Features: qu_ratio, esbl_act, broad_spec
DF_R2 <-  DF3 %>% dplyr::select(qu_ratio,esbl_act,NextQuRatio)

# Fit the NextQuRatio using the previous and broad_spec
fit.2 <- lm(DF_R2$NextQuRatio ~ DF_R2$esbl_act + DF_R2$qu_ratio)
fit.2
# graphical display
colors <- ifelse(DF_R2$esbl_act==2, "black", "gray")
plot(DF_R2$qu_ratio, DF_R2$NextQuRatio, xlab="Previous QuRatio", ylab="Next Qu",
     col=colors, pch=20)
curve(cbind(1, 1, x) %*% coef(fit.2), add=TRUE, col="black")
curve(cbind(1, 0,x) %*% coef(fit.2), add=TRUE, col="gray")

# allow for interactions: important because we expect 2nd variable to only be effective at non-zero esbl!
fit.3 <- lm(DF_R1$NextQuRatio ~ DF_R1$broad_spec + DF_R1$qu_ratio + DF_R1$broad_spec:DF_R1$qu_ratio)
fit.3
# graphical display
colors <- ifelse(DF_R1$broad_spec==1, "black", "gray")
plot(DF_R1$qu_ratio, DF_R1$NextQuRatio, xlab="Previous QuRatio", ylab="Next Qu",
     col=colors, pch=20)
curve(cbind(1, 1, x, 1*x) %*% coef(fit.3), add=TRUE, col="black")
curve(cbind(1, 0,x, 0*x) %*% coef(fit.3), add=TRUE, col="gray")

# Output: NextQuRatio -- Predicting Features: qu_ratio, esbl_act, broad_spec, TdiffToNext

### play around with data cleaning
SDATA


