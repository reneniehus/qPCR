setwd("~/Dropbox/RGNOSIS/Data/SATURN/Yehuda/")

rm(list=ls())
# load those libraries
library(dplyr)
library(car)
library(tidyr)
library(ggplot2)
# load the data file [this is the one from Pieter], put NA where there is emptiness or undetermined
DF <- read.csv("./Output16_CTX.csv",na.strings = c("", 'Undetermined'))
ex_dup <-read.csv("./marked_sample_ex.csv") # Samples which are re-run and should be excluded (on top of the ones Rene has
                                            # manually exluded at lines 32-77). Selection based on 1) whether marked yellow by
                                            # agatha. 2) when none of the re-runned samples marked yellow, the first run with
                                            # valid quantities is taken (i.e. quantity != NA)
#
# great habit to look at the first and the last couple of rows
head(DF)
tail(DF)
#View(DF)
# to be consistent, change all names in the header to small-case letters
names(DF) <- tolower(names(DF))
#
# Create DF2: only samples, no control samples, turn factor into character or numeric, adjust hyphones in sample names 
DF2 <- DF %>% 
  select(-pcr_type) %>%
  filter(class == "UNKNOWN",
         !(sample_name %in% c("quantity control sample","quantity control Rm1299"))) %>%
  mutate(sample_name=as.character(sample_name),
         sample_name=toupper(sample_name),
         sample_name=gsub('-','_',sample_name),
         type=as.character(type),
         class=as.character(class),
         ct = as.numeric(ct))
# Manually remove sample that she named invalid
DF2 <- DF2 %>%
  filter(!(DF2$sample_name == "SE_019_S1" & DF2$run_name == "20140426_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_019_S2" & DF2$run_name == "20140426_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_019_S3" & DF2$run_name == "20140426_16S_ATR_data.tmp"),
         !(DF2$sample_name == "IT_4119_S0" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "IT_4119_S1" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "IT_4119_S2" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "RM_3704_S3" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "RM_3704_S4" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "RM_3704_S5" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "RM_3704_SD" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_1922_S0" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_1922_S1" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_317_S0" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_317_S1" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_317_S2" & DF2$run_name == "20140506_16S_ATR_data.tmp"),
         !(DF2$sample_name == "SE_019_S1" & DF2$run_name == "20140320_CTX-M-A6-A8.tmp"),
         !(DF2$sample_name == "SE_019_S2" & DF2$run_name == "20140320_CTX-M-A6-A8.tmp"),
         !(DF2$sample_name == "SE_019_S3" & DF2$run_name == "20140320_CTX-M-A6-A8.tmp"),
         !(DF2$sample_name == "IT_4119_S0" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "IT_4119_S1" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "IT_4119_S2" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "RM_3704_S3" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "RM_3704_S4" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "RM_3704_S5" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "RM_3704_SD" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "SE_1922_S0" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "SE_1922_S1" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "SE_317_S0" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "SE_317_S1" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "SE_317_S2" & DF2$run_name == "20140507_CTX-M-A6-A8_ATR.tmp"),
         !(DF2$sample_name == "IT_3060_S7" & DF2$run_name == "20140718_CTX-M_LB_data.tmp"),
         !(DF2$sample_name == "RM_2148_S0" & DF2$run_name == "20140729_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_2654_S2" & DF2$run_name == "20140902_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_2654_S5" & DF2$run_name == "20140902_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_3864_S2" & DF2$run_name == "20140904_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_2864_S0" & DF2$run_name == "20140905_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_2864_S8" & DF2$run_name == "20140905_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_2148_S0" & DF2$run_name == "20140908_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_2654_S8" & DF2$run_name == "20140908_CTX-M_LB.tmp"),
         !(DF2$sample_name == "SE_275_S3" & DF2$run_name == "20140926_CTX-M_LB.tmp"),
         !(DF2$sample_name == "SE_537_S3" & DF2$run_name == "20140930_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_4332_S2" & DF2$run_name == "20141001_CTX-M_LB.tmp"),
         !(DF2$sample_name == "RM_4586_S1" & DF2$run_name == "20141001_CTX-M2_LB.tmp"))
# about above: i cannot spot the run "20140922_CTX-M_LB (2)" in the imported data
DF2 = DF2[order(DF2$sample_name,DF2$type),]

# Samples that have been diluted, multiply quantity by 10 or 100 (as diluted 1:10 or 1:100)
unique(DF2$sample_name[grep("1:10", DF2$sample_name)])
DF2$quantity_adj = DF2$quantity

DF2$quantity_adj[grep("1:10", DF2$sample_name)] = DF2$quantity[grep("1:10", DF2$sample_name)]*10
DF2$quantity_adj[grep("1:100", DF2$sample_name)] = DF2$quantity[grep("1:100", DF2$sample_name)]*100

# Which samples have been re-run
too_many_samp = (DF2 %>% group_by(sample_name, type) %>% summarise(reps = n()) %>% filter(reps>3))
#View(DF2[DF2$sample_name%in%too_many_samp$sample_name,])

DF_test = DF2 %>%
  filter(!(DF2$run_name%in%ex_dup$run & DF2$sample_name%in%ex_dup$sample_name ))

# Check with all unique samples have been included
length(unique(DF2$sample_name)) ; length(unique(DF_test$sample_name))
DF2$sample_name[which(!DF2$sample_name%in%(unique(DF_test$sample_name)))] # These samples are all excluded, where actually two shouldn't, as across present samples, all run and sample combinations are all true

# Add first sample runs of these two
DF2 = rbind(DF_test, DF2[which(!DF2$sample_name%in%(unique(DF_test$sample_name)))[c(1,2,5,6)],])
length(unique(DF2$sample_name))


# Create DF3: remove mean and sd collums and class, make groups based on sample and type, create my own means and standard deviations 
DF3 <- DF2 %>%
  select(-ct_mean, -ct_sd, -quantity_mean, - quantity_sd, -class) %>%
  group_by(sample_name, type) %>%
  summarise(ct_mean2=mean(ct,na.rm=TRUE),ct_sd2=sd(ct,na.rm=TRUE),qu_mean2=mean(quantity_adj,na.rm=TRUE),qu_sd2=sd(quantity_adj,na.rm=TRUE),
            qu_CV=sd(quantity_adj,na.rm=TRUE)/mean(quantity_adj,na.rm=TRUE), reps =n() ) %>%
  ungroup()

# Compare DF3 with DF2 in order to see where she calculated different means
names(DF3)
names(DF2)
# This is where I join the 2 DFs to compare
# test <- left_join(DF3, DF2, by=c('sample_name', 'type'))
test <- left_join(DF2, DF3)
test$diffctmean <- abs(test$ct_mean2 - test$ct_mean) > 1e-4
table(test$diffctmean)  # No difference now the extra samples have been removed

# put the 16s and the CTX-M data for one measure into same row, missing data will be NA
DF16s <- DF3 %>% filter(type=="16S")
DFCTX <- DF3 %>% filter(type=="CTX-M")
DF4 <- full_join(DF16s,DFCTX,by = "sample_name")
names(DF4) =  c("sample_name","type1","ct_mean_16s","ct_sd_16s","qu_mean_16s","qu_sd_16s","qu_CV_16s",
                "reps_16s","type2","ct_mean_CTX","ct_sd_CTX","qu_mean_CTX","qu_sd_CTX","qu_CV_CTX","reps_CTX")
DF4 <- DF4 %>%
  select(-type1, -type2)
# Add the quantity ratios
DF4$qu_ratio <- DF4$qu_mean_CTX/DF4$qu_mean_16s
# Add the error for the quantity ratio, relevant paper: doi:10.1016/j.clinbiochem.2006.12.014
DF4$ratio_CV <- sqrt(DF4$qu_CV_CTX^2 + DF4$qu_CV_16s^2 + 3*DF4$qu_CV_16s^2*DF4$qu_CV_CTX^2 + 8*DF4$qu_CV_16s^4)/(1 + DF4$qu_CV_16s^2)
#
# Data detectability: set ratio to zero when Ct_mean of CTXm > 30 
# my CV cutoff needs to 30% based on paper: doi:10.1016/j.clinbiochem.2006.12.014
# remove all rows with ratio NA
DF4 <- DF4 %>%
  filter(!is.na(qu_ratio))

DF5 <- DF4
ex1 = table((DF5$ct_mean_CTX) > 30 & (!is.na(DF5$ct_mean_CTX)))["TRUE"]
DF5$qu_ratio[(DF5$ct_mean_CTX) > 30 & (!is.na(DF5$ct_mean_CTX))] <- 0 # set the ratio to zero when CTXm is not detectable
#
ex2 = table((DF5$ct_mean_16s) > 30 & (!is.na(DF5$ct_mean_CTX)))["TRUE"]
DF5$qu_ratio[(DF5$ct_mean_16s) > 30 & (!is.na(DF5$ct_mean_CTX))] <- NA # set ratio to NA when 16s is not detectable
#
ex3 = table((DF5$qu_CV_16s > 0.30))["TRUE"]
DF5$qu_ratio[(DF5$qu_CV_16s > 0.30)] <- NA # if the 16s has a too high variance then we don't believe result
ex4 = table((DF5$qu_CV_CTX > 0.30) & !(DF5$ct_mean_CTX > 30))["TRUE"]
DF5$qu_ratio[(DF5$qu_CV_CTX > 0.30) & !(DF5$ct_mean_CTX > 30)] <- NA # if the CTX has too high varaince and is not set to zero
#
ex5 = table(is.na(DF5$ratio_CV))["TRUE"]
DF5$qu_ratio[is.na(DF5$ratio_CV)] <- NA # if the CTX has too high varaince and is not set to zero

tot_excluded = ex1 + ex2 + ex3 + ex4

# Add variable with patient id
DF5 <- DF5 %>%
  filter(!is.na(qu_ratio))

DF5$patient_id = gsub("S[D, 1,2,3,4,5,6,7,8,9,10, =D ]*","",DF5$sample_name)

# Add variable with sample number
DF6=NULL  
for(i in unique(DF5$patient_id)){
  d = DF5[DF5$patient_id == i,]
  d$s_num = c(1:length(d$sample_name))
  DF6 = rbind(DF6,d)
}

# Check for strange outliers
ggplot(DF6, aes(x=s_num, y=as.numeric(qu_ratio), group=patient_id))+geom_point()+geom_line()+facet_wrap(~patient_id,ncol=10)+ylim(0,120)

# Patient IT_3294 is only one with ratio >100 (sample IT_3294_S1). This sample has been tested 3 times, which might suggest something strange happening. 
# Also the 16s quantity is pretty low in this sample.


# export excel table
write.csv(DF6, file="./CleanedCTX_M16sRatioErr.csv")

# What distribution does the within PCR error follow
# Get to the distribution of the quantitiy control sample
DFQC1 <- DF %>%
  filter(sample_name == "quantity control sample",type == "16S")

DFQC2 <- DF %>%
  filter(sample_name == "quantity control sample",type == "CTX-M")

#check this out! 




