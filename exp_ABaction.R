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
# function that compares neighbouring elements of a data vector
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

### Mutations  
# turn into numerical
SDATA <- SDATA %>% 
  dplyr::mutate(
         num = as.numeric(num),
         ratio_SD = as.numeric(ratio_SD),
         s_num = as.numeric(s_num),
         qu_ratio = as.numeric(qu_ratio),
         ESBL_Ecoli = as.numeric(ESBL_Ecoli),
         ESBL_PPM = as.numeric(ESBL_PPM),
         ESBL_KESC = as.numeric(ESBL_KESC)
  )
# then also mutate the range of Antibiotics!
SDATA <- SDATA %>% dplyr::mutate_each(funs(as.numeric),esbl_act:Cefalexin)
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
# the median is 0.0077: majority of measures are close to zero.
# I expect that this ratio is mostly very small. 
# how accurate are small values of the ratio. Look at the error
summary(SDATA$ratio_SD) # gives a summery of a coloumn
# the SD of the ratio is around 1% of 
mean.quratio = mean(SDATA$qu_ratio, na.rm = T) # calculate the mean
max.quratio = max(SDATA$qu_ratio, na.rm = T) # calculate the mean

########
# small table version of SDATA
#names(SDATA)[90:126] # those are all the antibiotics in action
DF1 <- SDATA %>%
  dplyr::select(num,qu_ratio,Tdiff,RectalDate,ESBL_Ecoli,ESBL_KESC,ESBL_PPM,
                esbl_act:Cefalexin)

# replace all NA by 0
DF1[is.na(DF1)] <- 0

# remove 2 outliers [maybe outliers, why are they outliers?]
DF1 <- DF1[-which(max(DF1$qu_ratio) == DF1$qu_ratio),]
DF1 <- DF1[-which(max(DF1$qu_ratio) == DF1$qu_ratio),]
# turn esbl_act into binary [is this a good idea?]
DF1$esbl_act[DF1$esbl_act == 2] <- 1

# add the time since last measurement as predictor
# add a column with previous date
next.date <- c(DF1$RectalDate[2:(length(DF1$RectalDate))] , as.Date("2011-1-1"))
DF1$DaysToNextQuRatio <- next.date - DF1$RectalDate
DF1$DaysToNextQuRatio <- as.numeric(DF1$DaysToNextQuRatio)

# the next qu_ratio is what i want to predict
DF1$NextQuRatio = c(DF1$qu_ratio[2 : (nrow(DF1))], 0 )

# 
first.meas = which(DF1$num == 0)
first.meas <- first.meas[2 : length(first.meas)] - 1

# Filter away first measurements, de-select Tdiff
DF2 <- DF1[-first.meas,] %>% 
  dplyr::select(-Tdiff) %>%
  dplyr::filter(!DaysToNextQuRatio < 0)

# add a completely binary variable
DF2$RandBinary <- rbinom(nrow(DF2),1,0.5)

### Regression with interaction [important because most points are 
# small ratio -> small ratio, forces coefficient of binary to zero]
# fit.3: broad_spec, qu_ratio + interaction
DF_R1 <-  DF2 %>% dplyr::select(num,qu_ratio,broad_spec,NextQuRatio,RandBinary)
# centre and rescale
DF_R1 <- mutate(
  DF_R1,qu_ratio = (qu_ratio - mean(qu_ratio))/(2*(sd(qu_ratio))),
  broad_spec = (broad_spec - mean(broad_spec))/(2*(sd(broad_spec)))
)
fit.3 <- lm(DF_R1$NextQuRatio ~ DF_R1$broad_spec + DF_R1$qu_ratio + DF_R1$broad_spec:DF_R1$qu_ratio)
fit.3
# graphical display
colors <- ifelse(DF_R1$broad_spec>0, "black", "gray")
plot(DF_R1$qu_ratio, DF_R1$NextQuRatio, xlab="Previous QuRatio", ylab="Next Qu",
     col=colors, pch=20)
curve(cbind(1, 1, x, 1*x) %*% coef(fit.3), add=TRUE, col="black")
curve(cbind(1, 0,x, 0*x) %*% coef(fit.3), add=TRUE, col="gray")

# fit.4: esbl_act, qu_ratio + interaction
DF_R4 <-  DF2 %>% dplyr::select(num,qu_ratio,esbl_act,NextQuRatio,RandBinary) 
DF_R4$esbl_act[DF_R4$esbl_act == 2] <- 1 # make maybe = yes
fit.4 <- lm(DF_R4$NextQuRatio ~ DF_R4$esbl_act + DF_R4$qu_ratio + DF_R4$esbl_act:DF_R4$qu_ratio)
fit.4
# graphical display
colors <- ifelse(DF_R4$esbl_act==1, "black", "gray")
plot(DF_R4$qu_ratio, DF_R4$NextQuRatio, xlab="Previous QuRatio", ylab="Next Qu",
     col=colors, pch=20)
curve(cbind(1, 1, x, 1*x) %*% coef(fit.4), add=TRUE, col="black")
curve(cbind(1, 0,x, 0*x) %*% coef(fit.4), add=TRUE, col="gray")


# Output: NextQuRatio -- Predicting Features: 6 different ones
DF_R5 <-  DF2
#DF_R5 <- DF_R5 %>% dplyr::mutate_each(funs(as.numeric),esbl_act:Cefalexin)
#(qu_ratio - mean(qu_ratio))/(2*(sd(qu_ratio)))
fit.5 <- lm(formula = NextQuRatio ~ DaysToNextQuRatio + 
     qu_ratio + esbl_act + qu_ratio:RandBinary, data = DF_R5) 
summary(fit.5)
# make output readable


### Regressions without interaction term
# Output: NextQuRatio -- Predicting Features: qu_ratio, esbl_act, broad_spec
DF_R1 <-  DF2 %>% dplyr::select(num,qu_ratio,broad_spec,NextQuRatio)

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
