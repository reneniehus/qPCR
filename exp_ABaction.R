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

# calculate differences in the qu_ratio (first measure per patient will become meaningless)
diffqu <- SDATA$qu_ratio[2 : length(SDATA$sample_name)] - 
  SDATA$qu_ratio[1 : (length(SDATA$sample_name) - 1)]
diffqu <- c(0,diffqu)

# add diffqu to SDATA
SDATA$diffqu <- diffqu

########
# small table version of SDATA
DF1 <- SDATA %>%
  dplyr::select(s_num,patient_id,RectalDate,diffqu,qu_ratio,esbl_act,broad_spec)

# replace all NA by 0
DF1[is.na(DF1)] <- 0

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

## Get days since first measurement
# add a column with previous date
DF1$PrevDate <- c(as.Date("2011-1-1"), DF1$RectalDate[1:(length(DF1$RectalDate) - 1)])
DF1$Tdiff <- DF1$RectalDate - DF1$PrevDate
# add up the differences in days
for (i in 1:length(DF1$Tdiff)) {
  if (DF1$num[i] == 0) {
    DF1$Tdiff[i] <- 0
  } else {
    DF1$Tdiff[i] <- DF1$Tdiff[i] + DF1$Tdiff[i - 1]
  }
}
# check that for each patient first measurement is at 0 days
#ddply(DF1, .(patient_id), summarise, MinNum=min(Tdiff))





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
