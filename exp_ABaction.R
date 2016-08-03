##################################################
# Seaching for fingerprints of antibiotic action
##################################################

# Author: Rene Niehus
# Date: 1 August 2016

# clear workspace, load libraries
rm(list=ls())
library(dplyr)
library(reshape)
library(plyr)

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
  mutate(qu_ratio = as.numeric(qu_ratio),
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
  mutate(RectalDate = as.Date(as.character(SDATA$RectalDate), format="%Y-%m-%d"))

# order after rectal data FOR EACH patient
SDATA <- SDATA[order(SDATA$patient_id, SDATA$RectalDate),]

# calculate differences in the qu_ratio (first measure per patient will become meaningless)
diffqu <- SDATA$qu_ratio[2 : length(SDATA$sample_name)] - 
  SDATA$qu_ratio[1 : (length(SDATA$sample_name) - 1)]
diffqu <- c(0,diffqu)

# add diffqu to SDATA
SDATA$diffqu <- diffqu

# create numbers in order to remove the first samples
# create a counter for the samples from each patient, to identify first sample
cnt <- SDATA$patient_id[2:length(SDATA$patient_id)] ==
  SDATA$patient_id[1:(length(SDATA$patient_id)-1)]
cnt[cnt == "TRUE"] <- 1 # convert the counter to 0 and 1
cnt <- c(0, cnt) # add 0 for first element because this is first timepoint for the patient

# create new counter 0, 1, ..., sample_n of the same length
cntN <- vector(mode="numeric", length=length(cnt))
for (i in 1:length(cnt)){
  if (cnt[i] == 0){
    cntN[i] <- 0
  } else {
    cntN[i] <- cntN[i - 1] + 1
  }
}
SDATA$cnt <- cntN

###
# small table version of SDATA
DF1 <- SDATA %>%
  select(cnt,s_num,patient_id,RectalDate,diffqu,qu_ratio,esbl_act)

# CHECK if there are multiple measures from the same date
id.meas.a.p = Comp(DF1$s_num)*Comp(DF1$patient_id) # idential measure point and patient

doub.m = c(which(id.meas.a.p == 1), (which(id.meas.a.p == 1) + 1)) # 
doub.m <- doub.m[order(doub.m)]
# After correcting Italian IT317, here S1 and S14 (!!) have same rectalDate 2011-05-02 
DF1_dbls <- DF1[doub.m,]
# 5 Romanian samples RM2199, RM2337, RM2654, RM3913, RM4083. 
# Always timestamp 10 is double, where SD was taken. I suspect that at patient release (SD) 
# an extra sample was taken. And then sometimes 3 days later a follow up??
ddply(DF1_dbls, .(patient_id), summarise, mean=mean(qu_ratio),sigma=sd(qu_ratio),sigma/mean)

DF1_dbls$qu_r_diff <- c(diff(DF1_dbls$qu_ratio),0)/DF1_dbls$qu_ratio
DF1_dbls$zo <- rep(c(0,1),5)


DF1[, ] # there are quite a few 
  

# replace all NA by 0
DF1[is.na(DF1)] <- 0

# remove the first measument for each patient
DF1 <- DF1 %>%
  filter(!(DF1$cnt == 0))

# get to the times between measurements in days
DF1$PrevDate <- c("2011-1-1", DF1$RectalDate[1:(length(DF1$RectalDate) - 1)])

# reorder data frame
DF1 <- DF1 %>%
  select(cnt,patient_id,diffqu,esbl_act,RectalDate,PrevDate)

#
DF1$Tdiff <- as.Date(as.character(DF1$RectalDate), format="%Y-%m-%d") -
  as.Date(as.character(DF1$PrevDate), format="%Y-%m-%d")

DF1$Tdiff2 <- DF1$Tdiff

for (i in 1:length(DF1$cnt)) {
  if (DF1$cnt[i] == 1) {
    DF1$Tdiff[i] <- 0
  } else {
    DF1$Tdiff[i] <- DF1$Tdiff[i] + DF1$Tdiff[i - 1]
  }
}



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
