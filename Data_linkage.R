##################################################
# Linkage of qPCR data to clinical and abx data
##################################################

# Author: Esther van Kleef
# Date: 25 July 2016
rm(list=ls())

setwd("~/Documents/RGNOSIS/qPCR/")


# Load data files
clin_main <- read.csv ("./Raw_data/Main Data.csv", sep= ",", colClasses=c("character")) # Main clinical data, one row per patient
clin_fu <- read.csv("./Raw_data/Main FU Data.csv", sep=",", colClasses = c("character")) # Follow up clinical data (for follow up samples), multiple row per patients
lab_main <- read.csv("./Raw_data/Lab Main Data.csv", sep=",", colClasses = c("character")) # Main lab data, one row per patient. Actually, labmain and clinmain contain the same information I think
lab_fu <- read.csv ("./Raw_data/Lab FU Data.csv", sep= ",", colClasses=c("character")) # Lab follow up data, multiple rows per patient
pcr <- read.csv ("./Cleaned_data/CleanedCTX_M16sRatioErr.csv", sep= ",", colClasses=c("character")) # cleaned qPCR data, multiple rows per patient 
abx <- read.csv ("./Raw_data/Lab Ant Data.csv", sep= ",", colClasses=c("character")) # Abx use per patient
abxcat <- read.csv ("./Raw_data/Antibiotics_list.csv", sep= ",", colClasses=c("character")) # Used to categorise the abx used

# Create patient sample ID in lab_fu
lab_fu$Patient..ID = gsub("-", "", lab_fu$Patient..ID)
lab_fu$sample <- with(lab_fu, paste0(Country.Code, sep= "_", Patient..ID, sep= "_S", ScreeningNumber)) 
lab_fu$ESBL.16S <- gsub("%", "", lab_fu$ESBL.16S)
lab_fu$patient_id <- paste0(lab_fu$Country.Code,lab_fu$Patient..ID)

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


#################################
# MERGE DATA FRAMES

# Merge pcr data with lab_fu to obtain rectal dates
# Merge lab data with cleaned 16s values
length(pcr$sample_name[which(!pcr$sample_name %in%unique(lab_fu$sample) )])

# Check which variables of relevance
sapply(lab_fu, function(x) unique(x))
sapply(lab_fu, function(x) length(x[x==""]))

DF = merge(pcr, lab_fu[,-c(13:45)], by.x=c("sample_name","s_num","patient_id"), by.y=c("sample","ScreeningNumber","patient_id"),all.x = T)
View(DF[DF$sample_name%in%unique(DF$sample_name[which(duplicated(DF$sample_name))]),])
 
no_dup = which(!duplicated(DF$sample_name))

DF = data.frame(cbind(patient_id=DF$patient_id,country =DF$Country.Code,
                      sample_name=DF$sample_name,RectalDate=as.character(DF$RectalDate),as.data.frame(sapply(DF[,which(!names(DF)%in%c("patient_id","Country.Code","sample_name","RectalDate"))]
                                                                                                             ,function(x) as.numeric(x)))))

DF2 = DF %>%
  dplyr::select(-patient_id,-RectalDate,-country) %>%
  group_by(sample_name) %>% summarise_each(funs(mean))%>%
  mutate(patient_id=DF$patient_id[no_dup],RectalDate=DF$RectalDate[no_dup],
         country=DF$country[no_dup])
  
# Merge clin_main and pcr

# Check how similar clin_main and lab_main are 
which(!names(clin_main)%in%names(lab_main))
which(!names(lab_main)%in%names(clin_main))

sapply(lab_main[,-c(1)],function(x) length(x[x==""]))-sapply(clin_main[,-c(1,62)],function(x) length(x[x==""]))
# lab main is more complete 

length(unique(pcr$patient_id[which(!pcr$patient_id %in%unique(lab_main$patient_id) )])) # 11 are not in the main_lab file ?!
length(unique(pcr$patient_id[which(!pcr$patient_id %in%unique(clin_main$patient_id) )])) # 1 is not in the clin_main file

# Get the 11 patients missing in the lab data from the clinical data
add_on = clin_main[clin_main$patient_id%in%unique(pcr$patient_id[which(!pcr$patient_id %in%unique(lab_main$patient_id) )]),-c(1,62)]
 
clinNOTlab.names <- names(lab_main)[!names(lab_main)%in%names(add_on)]
add_on[,clinNOTlab.names ] <- NA
add_on = add_on[,order(colnames(add_on))]
lab_main = lab_main[,order(colnames(lab_main))]  

lab_main = rbind(lab_main, add_on)

DF3 = merge(pcr, lab_main, by.x=c("patient_id"), by.y=c("patient_id"), all.x=T)

# Merge Abx with Abx categories
colnames(abxcat)[which(names(abxcat)== "Antibiotic.ID")] <- "AntibioticID"
abx <- merge(abx,abxcat, by="AntibioticID")
abx$StartTreatmentDate = as.Date(abx$StartTreatmentDate,format="%d-%b-%y")
abx$EndTreatmentDate = as.Date(abx$EndTreatmentDate,format="%d-%b-%y")

# Merge Abx with pcr data
DF4 = merge(DF3, abx[,-c(1,2)], by.x=c("patient_id"),by.y=c("patient_id"), all.x=T)

DF4[,which(names(DF4)%in%c("RectalDate","StartTreatmentDate", 
                           "EndTreatmentDate","AdmittanceDateInHospital"))] = lapply(DF4[,which(names(DF4)%in%c("RectalDate","StartTreatmentDate", 
                                                                                                                "EndTreatmentDate","AdmittanceDateInHospital"))],
                                                                                     function(x) as.Date(x,format="%d-%b-%y"))


######################
# Check with plot NEED TO CHANGE RECTAL DATE STILL!
png(filename="./Figures/ratio_with_abx.png", width=1000, height=1000)

i = ggplot(DF4, aes(x=RectalDate,y=as.numeric(qu_ratio)))+ geom_segment(aes(x = StartTreatmentDate, y = -5, xend = EndTreatmentDate, yend = -5, size=1,colour = betalactamase_activity),data=DF4[DF4$betalactamase_activity%in%"Yes",])+ 
  geom_segment(aes(x = StartTreatmentDate, y = -7, xend = EndTreatmentDate, yend = -7,size=1, colour = betalactamase_activity),data=DF4[DF4$betalactamase_activity%in%"No",])+
  geom_segment(aes(x = StartTreatmentDate, y = -9, xend =EndTreatmentDate, yend = -9, size=1,colour = betalactamase_activity),data=DF4[DF4$betalactamase_activity%in%"Possibly",])+geom_line()+geom_point()+
  scale_colour_manual(name="Antibiotic ESBL activity", values=c("deepskyblue","deeppink","green"))+
  facet_wrap(~patient_id, scales="free",ncol=6)+ylab("% abundance ESBL to 16s")
print(i)
dev.off()


