##################################################
# Linkage of qPCR data to clinical and abx data
##################################################

# Author: Esther van Kleef with changes by Rene Niehus
# Date: 25 July 2016

rm(list=ls())
require(dplyr)
library(ggplot2);

# set working directory
#setwd("~/Documents/RGNOSIS/qPCR/") # for Esther
setwd("~/Dropbox/LOMHWRU_MORU/SATURN_ESBL/_R/R_git/qPCR/") # for Rene

# Load data files
clin_main <- read.csv ("./Raw_data/Main Data.csv", sep= ",", colClasses=c("character")) # Main clinical data, one row per patient
clin_fu <- read.csv("./Raw_data/Main FU Data.csv", sep=",", colClasses = c("character")) # Follow up clinical data (for follow up samples), multiple row per patients
lab_main <- read.csv("./Raw_data/Lab Main Data.csv", sep=",", colClasses = c("character")) # Main lab data, one row per patient. Actually, labmain and clinmain contain the same information I think
lab_fu <- read.csv ("./Raw_data/Lab FU Data.csv", sep= ",", colClasses=c("character")) # Lab follow up data, multiple rows per patient
pcr <- read.csv ("./Cleaned_data/CleanedCTX_M16sRatioErr.csv", sep= ",", colClasses=c("character")) # cleaned qPCR data, multiple rows per patient 
abx <- read.csv ("./Raw_data/Lab Ant Data.csv", sep= ",", colClasses=c("character")) # Abx use per patient
abxcat <- read.csv ("./Raw_data/Antibiotics_list.csv", sep= ",", colClasses=c("character")) # Used to categorise the abx used
dates_pcr <-read.csv("./Raw_data/SATURN-WP5-sample list for qPCR.csv",sep= ",", colClasses=c("character"))

# Create patient sample ID in lab_fu
lab_fu$Patient..ID = gsub("-", "", lab_fu$Patient..ID)
lab_fu$sample <- with(lab_fu, paste0(Country.Code, sep= "_", Patient..ID, sep= "_S", as.numeric(ScreeningNumber)-1)) # ScreeningNumber - 1 as I think that in lab_fu the count starts with 1 whereas in pcr data with 0
lab_fu$s_num_clean = as.numeric(lab_fu$ScreeningNumber)-1
lab_fu$ESBL.16S <- gsub("%", "", lab_fu$ESBL.16S)
lab_fu$patient_id <- paste0(lab_fu$Country.Code,lab_fu$Patient..ID)

# Create sample ID in dates_PCR
dates_pcr$sample = paste0(dates_pcr$Country,"_",dates_pcr$patient_id_num,"_",dates_pcr$S_num)
dates_pcr$patient_id = paste0(dates_pcr$Country,dates_pcr$patient_id_num)

# Create s_num_clean from the sample names
screen_n = c("S0","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14")
dates_pcr$s_num_clean=rep(NA,1,length(dates_pcr$sample)) # initiate dates_pcr
for(i in unique(screen_n)){
  dates_pcr$s_num_clean[grep(i,dates_pcr$sample)] = which(screen_n==i)-1
}
# For the ones with sample number = SD NAs where produce, as these are discharge samples, samplenumber would be the follow up to the last one
dates_pcr$s_num_clean[which(is.na(dates_pcr$s_num_clean))] = dates_pcr$s_num_clean[which(is.na(dates_pcr$s_num_clean))-1]+1

# Clean Patient ID variable clin_main
clin_main$CountryCode_full = ifelse(clin_main$CountryCode == 1,"NA",
                                    ifelse(clin_main$CountryCode==2,"IT",
                                           ifelse(clin_main$CountryCode==3,"RM","SE")))

clin_main$patient_id = clin_main$PatientStudyID
clin_main$patient_id = gsub(c("Se|se|SE "), "SE", clin_main$patient_id)
clin_main$patient_id = gsub("\\(|\\)", "", clin_main$patient_id)
clin_main$patient_id = gsub("/", "", clin_main$patient_id)
clin_main$patient_id = gsub(" ", "", clin_main$patient_id)
clin_main$Country_mis =ifelse(grepl(c("SE|RM|IT"),clin_main$patient_id)==0,1,0)
clin_main$patient_id = ifelse(clin_main$Country_mis==1,paste0(clin_main$CountryCode_full,clin_main$patient_id),clin_main$patient_id)

clin_main = clin_main[order(clin_main$patient_id),]
# Check number of unique id's
#length(unique(clin_main$PatientStudyID));length(unique(clin_main$patient_id))
#View(clin_main[clin_main$patient_id%in%unique(clin_main$patient_id[which(duplicated(clin_main$patient_id))]),])

# 7 unique IDs less, as not recognised as same sample in original data due to brackets. Probably want 
# to include the observations without the brackets as the brackets are probably referring to the duplication


# Clean Patient ID variable lab_main
lab_main$PatientStudyID = gsub("-", "", lab_main$PatientStudyID)
lab_main$patient_id = lab_main$PatientStudyID
lab_main$patient_id = gsub(c("Se|se|SE "), "SE", lab_main$patient_id)
lab_main$patient_id = gsub("\\(|\\)", "", lab_main$patient_id)
lab_main$patient_id = gsub("/", "", lab_main$patient_id)
lab_main$patient_id = gsub(" ", "", lab_main$patient_id)

lab_main$Country_mis =ifelse(grepl(c("SE|RM|IT"),lab_main$patient_id)==0,1,0)
lab_main$patient_id = ifelse(lab_main$Country_mis==1,paste0(lab_main$CountryCode,lab_main$patient_id),lab_main$patient_id)

lab_main = lab_main[order(lab_main$patient_id),]
length(unique(lab_main$patient_id));length(unique(lab_main$PatientStudyID))

# Create abx patient ID variable
abx$patient_id = paste0(abx$Country.Code,abx$Patient.ID)

# In pcr remove =D
pcr$sample_name2 = gsub("_S[D, 1,2,3,4,5,6,7,8,9,10, =D ]*","",pcr$sample_name)
pcr$sample_name2 = paste0(pcr$sample_name2,"_S",pcr$s_num)
#################################
# MERGE DATA FRAMES

# Merge pcr data with dates_pcr
length(pcr$sample_name[which(!pcr$sample_name %in%unique(dates_pcr$sample) )])
DF = merge(pcr,dates_pcr[,which(names(dates_pcr)%in%c("patient_id","s_num_clean","RectalDate"))],by.x=c("patient_id","s_num"), by.y=c("patient_id","s_num_clean"), all.x=T)
# Check merge
length(DF$RectalDate[is.na(DF$RectalDate)])

# Merge pcr data with lab_fu 
length(unique(lab_fu$sample));length(unique(pcr$sample_name))
length(pcr$sample_name[which(!pcr$sample_name %in%unique(lab_fu$sample) )])# 129 samples not in the lab_follow up dates
pcr$sample_mis = ifelse(!pcr$sample_name %in%unique(lab_fu$sample),1,0)

length(unique(pcr$patient_id[which(!pcr$patient_id %in%unique(lab_fu$patient_id))])) # 11 unique patients not in the lab_fu data

# Check which variables of relevance
sapply(lab_fu, function(x) unique(x))
sapply(lab_fu, function(x) length(x[x==""]))

DF = merge(DF, lab_fu[,-c(4,15:45)], by.x=c("s_num","patient_id"), by.y=c("s_num_clean","patient_id"),all.x = T)

# convert qu_ratio to numeric
DF <- DF %>% mutate(qu_ratio = as.numeric(qu_ratio))

length(unique(lab_fu$sample)) # Duplicates present in the lab_fu file, therefore multiple merges
#View(lab_fu[lab_fu$sample%in%unique(lab_fu$sample[which(duplicated(lab_fu$sample))]),])

no_dup = which(!duplicated(DF$sample_name))

DF = data.frame(cbind(patient_id=DF$patient_id,country =DF$Country.Code,
                      sample_name=DF$sample_name,sample_name2=DF$sample_name2,RectalDate=as.character(DF$RectalDate),DischargeDate=as.character(DF$PatientDischargeDate),as.data.frame(sapply(DF[,which(!names(DF)%in%c("patient_id","Country.Code","sample_name","RectalDate", 
                                                                                                                                                                                                                        "sample_name2","PatientDischargeDate"))]
                                                                                                             ,function(x) as.numeric(x)))))
# ignore errors, is due to applying as.numeric to empty cells

DF2 = DF %>%
  dplyr::select(-patient_id,-RectalDate,-DischargeDate,-country,-sample,-sample_name2) %>%
  group_by(sample_name) %>% summarise_each(funs(median)) 
a = as.data.frame(cbind(patient_id=as.character(DF$patient_id[no_dup]),sample_name=as.character(DF$sample_name[no_dup]),
                        sample_name2=as.character(DF$sample_name2[no_dup]),RectalDate=as.character(DF$RectalDate[no_dup])))

a = a[order(a$sample_name),]
DF2 = merge(DF2,a,by.x=c("sample_name"),by.y=c("sample_name"))

length(unique(DF2$patient_id[which(is.na(DF2$RectalDate))]))  

# Merge clin_main and pcr

# Check how similar clin_main and lab_main are 
which(!names(clin_main)%in%names(lab_main))
which(!names(lab_main)%in%names(clin_main))

sapply(lab_main[,-c(1)],function(x) length(x[x==""]))-sapply(clin_main[,-c(1,62)],function(x) length(x[x==""]))
# lab main is more complete 

length(unique(pcr$patient_id[which(!pcr$patient_id %in%unique(lab_main$patient_id) )])) # 11 are not in the main_lab file
length(unique(pcr$patient_id[which(!pcr$patient_id %in%unique(clin_main$patient_id) )])) # 1 is not in the clin_main file

# Get the 11 patients missing in the lab data from the clinical data
add_on = clin_main[clin_main$patient_id%in%unique(pcr$patient_id[which(!pcr$patient_id %in%unique(lab_main$patient_id) )]),-c(1,62)]
 
clinNOTlab.names <- names(lab_main)[!names(lab_main)%in%names(add_on)]
add_on[,clinNOTlab.names ] <- NA
add_on = add_on[,order(colnames(add_on))]
lab_main = lab_main[,order(colnames(lab_main))]  

lab_main = rbind(lab_main, add_on)

DF3 = merge(DF2, lab_main[,-c(6,58)], by.x=c("patient_id"), by.y=c("patient_id"), all.x=T)

# Merge Abx with Abx categories
colnames(abxcat)[which(names(abxcat)== "Antibiotic.ID")] <- "AntibioticID"
abx <- merge(abx,abxcat, by="AntibioticID")
abx$StartTreatmentDate = as.Date(abx$StartTreatmentDate,format="%d-%b-%y")
abx$EndTreatmentDate = as.Date(abx$EndTreatmentDate,format="%d-%b-%y")

# Merge Abx with pcr data
abx=abx[order(abx$patient_id),]
rows = length(DF3$sample_name[which(DF3$patient_id%in%unique(abx$patient_id))])

ab_matrix = data.frame(matrix(NA,nrow=rows,ncol= length(unique(abx$Antibiotic.Name))+5))
colnames(ab_matrix) = c("patient_id","sample_name", "RectalDate","esbl_act","broad_spec",unique(abx$Antibiotic.Name))

pat_id = which(DF3$patient_id%in%unique(abx$patient_id))
ab_matrix[,1] = as.character(DF3$patient_id[pat_id])
ab_matrix[,2] = as.character(DF3$sample_name[which(DF3$patient_id%in%unique(abx$patient_id))])
ab_matrix[,3] = as.Date(DF3$RectalDate[pat_id],format="%d-%b-%y")
ab_matrix[,c(4,5)] = 0

for(i in unique(abx$patient_id)){
  d = abx[abx$patient_id==i,]
  for(a in unique(d$Antibiotic.Name)){
    a_rows = d[which(d$Antibiotic.Name==a),]
    if(length(a_rows$StartTreatmentDate)>1){
      dif_t = NULL
      for(p in 1:length(a_rows$StartTreatmentDate)){
      dif = c(a_rows$StartTreatmentDate[p]+0:(a_rows$EndTreatmentDate[p] - a_rows$StartTreatmentDate[p]))
      dif_t = append(dif_t,dif, after=length(dif_t))
      }
      treat = sort(dif_t)
    }
    else{
      treat = c(a_rows$StartTreatmentDate[1]+0:(a_rows$EndTreatmentDate[1] - a_rows$StartTreatmentDate[1]))
    }
    b = ab_matrix$RectalDate[which(ab_matrix$patient_id==i)]
    ab_matrix[which(ab_matrix$patient_id==i),which(names(ab_matrix)==a)] = ifelse(b%in%treat,1,0)
  }
}

#View(ab_matrix[ab_matrix$sample_name%in%unique(ab_matrix$sample_name[which(duplicated(ab_matrix$sample_name))]),])

ab_matrix[is.na(ab_matrix)] = 0

esbl = which(colnames(ab_matrix)%in%unique(abx$Antibiotic.Name[abx$betalactamase_activity=="Yes"]));length(esbl)
esbl_p = which(colnames(ab_matrix)%in%unique(abx$Antibiotic.Name[abx$betalactamase_activity=="Possibly"]));length(esbl_p)
b_broad = which(colnames(ab_matrix)%in%unique(abx$Antibiotic.Name[abx$Spectrum=="Broad"])) ;length(b_broad)

esbl_rows <- as.vector(unique(unlist(sapply(ab_matrix[,esbl], function(x) which(x==1), simplify=TRUE))))
esbl_p_rows <- as.vector(unique(unlist(sapply(ab_matrix[,esbl_p], function(x) which(x==1), simplify=TRUE))))
b_broad_rows <- as.vector(unique(unlist(sapply(ab_matrix[,b_broad], function(x) which(x==1), simplify=TRUE))))

ab_matrix$esbl_act[esbl_p_rows] = 2 # No = 0; Yes = 1; Possibly = 2
ab_matrix$esbl_act[esbl] = 1
ab_matrix$broad_spec[b_broad_rows] = 1 # No = 0; Yes = 1

DF4 = merge(DF3, ab_matrix[,-c(1,3)], by.x=c("sample_name"),by.y=c("sample_name"), all.x=T)

DF4[,which(names(DF4)%in%c("RectalDate","StartTreatmentDate", 
                           "EndTreatmentDate","AdmittanceDateInHospital"))] = lapply(DF4[,which(names(DF4)%in%c("RectalDate","StartTreatmentDate", 
                                                                                                                "EndTreatmentDate","AdmittanceDateInHospital"))],
                                                                                     function(x) as.Date(x,format="%d-%b-%y"))

DF4 <- DF4 %>%
  dplyr::select(-X,-Patient..ID,-ScreeningNumber,-reps_16s,-reps_CTX,-sample_name2,-PatientStudyID, -SamplePatID)

# I've left the abx fields for the patient that presumably haven't taken abx at all (i.e. were not present in the abx data file) is NA, as I suppose we're not 
# certain whether their usage is missing or really no abx was taken

# sort patient samples per patient after Rectal date
DF4 <- DF4[order(DF4$patient_id, DF4$RectalDate),]

# add num which always starts at 0
DF4$num <- DF4$s_num
DF4$num[1] <- 0
for (i in 2:length(DF4$num)) {
  if (DF4$patient_id[i] != DF4$patient_id[i - 1]) {
    DF4$num[i] <- 0
  } else {
    DF4$num[i] <- DF4$num[i - 1] + 1 
  }
}
# check that each patients num starts at 0
#ddply(DF4, .(patient_id), summarise, MinNum=min(num))

## Get days since first measurement
# add a column with previous date
prev.date <- c(as.Date("2011-1-1"), DF4$RectalDate[1:(length(DF4$RectalDate) - 1)])
DF4$Tdiff <- DF4$RectalDate - prev.date
# add up the differences in days
for (i in 1:length(DF4$Tdiff)) {
  if (DF4$num[i] == 0) {
    DF4$Tdiff[i] <- 0
  } else {
    DF4$Tdiff[i] <- DF4$Tdiff[i] + DF4$Tdiff[i - 1]
  }
}
# check that for each patient first measurement is at 0 days
#ddply(DF4, .(patient_id), summarise, MinDays=min(Tdiff), MaxDays=max(Tdiff))

# calculate between measurement differences in the qu_ratio (first measure per patient will become meaningless)
diff.qu <- DF4$qu_ratio[2 : nrow(DF4)] - DF4$qu_ratio[1 : (nrow(DF4) - 1)]
diff.qu <- c(0,diff.qu)
# set diff.qu for first measurement to NA
diff.qu[DF4$num == 0] <- NA

# add diffqu to DF
DF4$DiffRatio <- diff.qu

DF4_check = merge(DF4, abx,by.x=c("patient_id"), by.y=c("patient_id"))

######################
# Check with plot (this is only plotting the ones with abx usage present)

png(filename="./Figures/ratio_with_abx.png", width=1500, height=1200)

i = ggplot(DF4_check, aes(x=RectalDate,y=as.numeric(qu_ratio), group=patient_id))+ geom_line()+geom_point()+geom_hline(yintercept=0,linetype=2)+
geom_segment(aes(x = StartTreatmentDate, y = -5, xend = EndTreatmentDate, yend = -5, size=1,colour = betalactamase_activity),data=DF4_check[DF4_check$betalactamase_activity%in%"Yes",])+ 
  geom_segment(aes(x = StartTreatmentDate, y = -7, xend = EndTreatmentDate, yend = -7,size=1, colour = betalactamase_activity),data=DF4_check[DF4_check$betalactamase_activity%in%"No",])+
  geom_segment(aes(x = StartTreatmentDate, y = -9, xend =EndTreatmentDate, yend = -9, size=1,colour = betalactamase_activity),data=DF4_check[DF4_check$betalactamase_activity%in%"Possibly",])+
  scale_colour_manual(name="Antibiotic ESBL activity", values=c("deepskyblue","deeppink","green"))+
  facet_wrap(~patient_id, scales="free",ncol=10)+ylab("% abundance ESBL to 16s")
print(i)
dev.off()

png(filename="./Figures/ratio_with_abx_below1only.png", width=1500, height=1200)

i = ggplot(DF4_check[DF4_check$qu_ratio<1,], aes(x=RectalDate,y=as.numeric(qu_ratio), group=patient_id))+ geom_line()+geom_point()+geom_hline(yintercept=0,linetype=2)+
  geom_segment(aes(x = StartTreatmentDate, y = -0.5, xend = EndTreatmentDate, yend = -0.5, size=0.5,colour = betalactamase_activity),data=DF4_check[DF4_check$betalactamase_activity%in%"Yes",])+ 
  geom_segment(aes(x = StartTreatmentDate, y = -1, xend = EndTreatmentDate, yend = -1,size=0.5, colour = betalactamase_activity),data=DF4_check[DF4_check$betalactamase_activity%in%"No",])+
  geom_segment(aes(x = StartTreatmentDate, y = -1.5, xend =EndTreatmentDate, yend = -1.5, size=0.5,colour = betalactamase_activity),data=DF4_check[DF4_check$betalactamase_activity%in%"Possibly",])+
  scale_colour_manual(name="Antibiotic ESBL activity", values=c("deepskyblue","deeppink","green"))+
  facet_wrap(~patient_id, scales="free_x",ncol=10)+ylab("% abundance ESBL to 16s")+ylim(-2,1)
print(i)
dev.off()

DF_plot <- DF4_check[DF4_check$qu_ratio<1,]
png(filename="./Figures/ratio_with_abx_below1only_broadSpec.png", width=3000, height=2400)
i = ggplot(DF_plot, aes(x=RectalDate, y=qu_ratio, group=patient_id)) + 
  geom_line(size=1) + geom_point(size=1) + geom_hline(yintercept=0,linetype=2) +
  geom_segment(aes(x = StartTreatmentDate, y = -0.5, xend = EndTreatmentDate, yend = -0.5, size=0.5, colour = as.character(broad_spec)),data=DF_plot[DF_plot$broad_spec==1,]) + 
  geom_segment(aes(x = StartTreatmentDate, y = -1, xend = EndTreatmentDate, yend = -1,size=0.5, colour = as.character(broad_spec)),data=DF_plot[DF_plot$broad_spec==0,]) +
  scale_colour_brewer(name="Broad Spectrum", palette="Pastel2") + theme(axis.title = element_text(size = 9), 
    plot.title = element_text(size = 15), 
    panel.background = element_rect(fill = "gray95"))+
  facet_wrap(~patient_id, scales="free_x",ncol=10)+ylab("% abundance ESBL to 16s")+ylim(-2,1)
print(i)
dev.off()



write.csv(DF4,file="./Cleaned_data/linked_qPCR_clin_abx.csv")
