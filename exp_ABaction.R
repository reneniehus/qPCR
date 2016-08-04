##################################################
# Seaching for fingerprints of antibiotic action
##################################################

# Author: Rene Niehus
# Date: 1 August 2016

# clear workspace, load libraries
rm(list=ls())
library(plyr)
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

# Output: NextQuRatio -- Predicting Features: qu_ratio, esbl_act, broad_spec, TdiffToNext




# CHECK if there are multiple measures from the same date
id.meas.a.p = Comp(DF1$s_num)*Comp(DF1$patient_id) # idential measure point and patient
doub.m = c(which(id.meas.a.p == 1), (which(id.meas.a.p == 1) + 1)) # 
doub.m <- doub.m[order(doub.m)]
# After correcting s_num, no more doublicates on same date [before: IT317, RM2199, RM2337, RM2654, RM3913, RM4083]

# add num which always starts at 0
DF1$num <- DF1$s_num
DF1$num[1] <- 0
for (i in 2:length(DF1$num)) {
  if (DF1$patient_id[i] != DF1$patient_id[i - 1]) {
    DF1$num[i] <- 0
  } else {
    DF1$num[i] <- DF1$num[i - 1] + 1 
  }
}
# check that each patients num starts at 0
#ddply(DF1, .(patient_id), summarise, MinNum=min(num))


#### FOR DIFFQU: remove first measument for each patient, because 
DF1 <- DF1 %>%
  dplyr::filter(!(DF1$s_num == 0))



# reorder data frame
DF1 <- DF1 %>%
  select(cnt,patient_id,diffqu,esbl_act,RectalDate,PrevDate)

DF1$Tdiff2 <- DF1$Tdiff





# make diffqu into -1 0 and +1
DFdiff <- DF1 %>%
  select(diffqu) # make a new table DFdiff that contains the steps to simplify

# select where to make cut-off for no change in CTX-m abundance
summary_of_quDiff = summary(DFdiff$diffqu)
small_vals = (-0.02 <= DFdiff$diffqu & DFdiff$diffqu <= 0.02) # introduce small value-variable

DFdiff$qu3groups <- DFdiff$diffqu # make a simple version of qu_diff
DFdiff$qu3groups[small_vals] <- 0
## 1: increase -1:decrease 0: stays
DFdiff$qu3groups[DFdiff$qu3groups > 0] <- 1
DFdiff$qu3groups[DFdiff$qu3groups < 0] <- -1

# add the simple version od qu_diff to DF1
DF1$qu3groups <- DFdiff$qu3groups

# fit a model using esbl active and broad spec as X variables
model1 <- lm(DF1$qu3groups ~ DF1$esbl_act + DF1$broad_spec)
summary(model1)

# calculate the Pear
cor(DF1$esbl_act, DF1$broad_spec, method="pearson")

# conf intervals for model coefficients
confint(model1, conf.level = 0.95)

### Model using all X variables
# small table version of SDATA
DF2 <- SDATA %>%
  select(cnt,patient_id,sample_name,diffqu,esbl_act,broad_spec,
         Piperacillin.tazobactam,Vancomycin,isoniazid,Levofloxacin,
         Meropenem,pyrazinamide,Rifampicin,Ampicillin,Oxacillin,
         Ceftriaxone)

# replace all NA by 0
DF2[is.na(DF2)] <- 0

# make treaments numeric

# remove the first measument for each patient
DF2 <- DF2 %>%
  filter(!(DF2$cnt == 1))

# make diffqu into -1 0 and +1
DFdiff <- DF2 %>%
  select(diffqu) # make a new table DFdiff that contains the steps to simplify

# select where to make cut-off for no change in CTX-m abundance
summary_of_quDiff = summary(DFdiff$diffqu)
small_vals = (-0.02 <= DFdiff$diffqu & DFdiff$diffqu <= 0.02) # introduce small value-variable

DFdiff$qu3groups <- DFdiff$diffqu # make a simple version of qu_diff
DFdiff$qu3groups[small_vals] <- 0
## 1: increase -1:decrease 0: stays
DFdiff$qu3groups[DFdiff$qu3groups > 0] <- 1
DFdiff$qu3groups[DFdiff$qu3groups < 0] <- -1

# add the simple version od qu_diff to DF2
DF2$qu3groups <- DFdiff$qu3groups

# fit a model using esbl active and broad spec as X variables
model2 <- lm(DF2$qu3groups ~ DF2$broad_spec +
               DF2$Piperacillin.tazobactam + DF2$Vancomycin +
               DF2$Levofloxacin + DF2$Meropenem +
               DF2$Rifampicin + DF2$Ampicillin +
               DF2$Oxacillin + DF2$Ceftriaxone)
summary(model2)

# calculate the Pearson correlation
cor(DF1$esbl_act, DF1$broad_spec, method="pearson")

# conf intervals for model coefficients
confint(model1, conf.level = 0.95)


### Model using all X variables and stratify
# small table version of SDATA
DF3 <- SDATA %>%
  select(ESBL_Ecoli,ESBL_KESC,ESBL_PPM,cnt,patient_id,sample_name,diffqu,esbl_act,broad_spec,
         Piperacillin.tazobactam,Vancomycin,isoniazid,Levofloxacin,
         Meropenem,pyrazinamide,Rifampicin,Ampicillin,Oxacillin,
         Ceftriaxone)

# make treaments numeric

# replace all NA by 0
DF3[is.na(DF3)] <- 0

# remove the first measument for each patient
DF3 <- DF3 %>%
  filter(!(DF3$cnt == 1))

# make diffqu into -1 0 and +1
DFdiff <- DF3 %>%
  select(diffqu) # make a new table DFdiff that contains the steps to simplify

# select where to make cut-off for no change in CTX-m abundance
summary_of_quDiff = summary(DFdiff$diffqu)
small_vals = (-0.02 <= DFdiff$diffqu & DFdiff$diffqu <= 0.02) # introduce small value-variable

DFdiff$qu3groups <- DFdiff$diffqu # make a simple version of qu_diff
DFdiff$qu3groups[small_vals] <- 0
## 1: increase -1:decrease 0: stays
DFdiff$qu3groups[DFdiff$qu3groups > 0] <- 1
DFdiff$qu3groups[DFdiff$qu3groups < 0] <- -1

# add the simple version od qu_diff to DF2
DF3$qu3groups <- DFdiff$qu3groups

# fit a model using esbl active and broad spec as X variables
model3 <- lm(DF3$qu3groups ~ DF3$broad_spec +
               DF3$Piperacillin.tazobactam + DF3$Vancomycin +
               DF3$Levofloxacin + DF3$Meropenem +
               DF3$Rifampicin + DF3$Ampicillin +
               DF3$Oxacillin + DF3$Ceftriaxone)
summary(model3)

# calculate the Pearson correlation
cor(DF1$esbl_act, DF1$broad_spec, method="pearson")

# conf intervals for model coefficients
confint(model1, conf.level = 0.95)
