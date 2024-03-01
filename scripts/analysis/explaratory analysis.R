---
title: "Exploratory Analysis"
output: html_notebook
---

#install these too
```{r, message=FALSE}
library(survival)
library(survminer)

library(tidyverse)
library(ggpubr)
library(dplyr)
library(survival); library(survminer)

setwd("~/Documents/GGCT_BRCA_assocation")

```

#SW's code
```{r}
patient.ICI <- patient.WT %>% subset(ici_YN == "1")
patient.ICI$streoid <- relevel(as.factor(patient.ICI$steroid), ref="no")

fit = survfit(Surv(OS, status_OS) ~ ici_cat, data= patient.WT)
plot(fit)
ggsurvplot(fit, data=patient.WT, pval=T)
coxph()


fit.ICI <- survfit(Surv(OS, status_OS) ~ steroid, data=patient.ICI)
ggsurvplot(fit.ICI, facet.by="ici_cat", data=patient.ICI, pval=T)

cox.ICI <- coxph(Surv(OS, status_OS) ~ age + sex+ EOR+ ici_cat + steroid + KPS, data=patient.ICI)
ggforest(cox.ICI, data=patient.ICI)

```


#read in bed coverage output
```{r, quiet=T}
file_list <- list.files(path="~/Documents/GGCT_BRCA_assocation/bedtools_coverage_out", pattern = "*.ggct.bed", full.names = F)

#read in bedfile, mutate and add extra column for name of bed file
mutDF = data.frame(matrix(nrow=1231, ncol=29))

for (startFile in ((file_list))) { 
  
  #remove .ggct.bed to get the file name
  myFile=str_split(startFile, ".ggct.bed")[[1]][1]
  
  # Read the file and create a data frame
  thisFile =  paste0("/Users/agalianese/Documents/GGCT_BRCA_assocation/bedtools_coverage_out/", startFile)
  myData <- data.frame(read_tsv(file=thisFile, col_names = F, show_col_types = FALSE))
  
  i=which(str_starts(file_list, myFile))

  mutDF[i, 1:4] = myData[7,7:10]
  mutDF[i, 5:8] = myData[6,7:10]
  mutDF[i, 9:12] = myData[5,7:10]
  mutDF[i, 13:16] = myData[4,7:10]
  mutDF[i, 17:20]  <- myData[3,7:10]
  mutDF[i, 21:24]  <- myData[1,7:10]
  mutDF[i, 25:28] <- myData[2,7:10]
  mutDF[i,29] <- myFile
}

colnames(mutDF) <- c("E1.1.N.Features", "E1.1.A.Bases", "E1.1.A.Len", "E1.1.Coverage.Ratio", 
                     "E1.2.N.Features", "E1.2.A.Bases", "E1.2.A.Len", "E1.2.Coverage.Ratio", 
                     "E2.N.Features", "E2.A.Bases", "E2.A.Len", "E2.Coverage.Ratio", 
                     "E3.1.N.Features", "E3.1.A.Bases", "E3.1.A.Len", "E3.1.Coverage.Ratio", 
                     "E3.2.N.Features", "E3.2.A.Bases", "E3.2.A.Len", "E3.2.Coverage.Ratio", 
                     "E4.1.N.Features", "E4.1.A.Bases", "E4.1.A.Len", "E4.1.Coverage.Ratio", 
                     "E4.2.N.Features",  "E4.2.A.Bases", "E4.2.A.Len", "E4.2.Coverage.Ratio", "fileId")

mutDF


```

#create dataframe for easy ggplot (ggDF)
```{r}

ggDF = data.frame()
for (startFile in ((file_list))) { 
  
  #remove .ggct.bed to get the file name
  myFile=str_split(startFile, ".ggct.bed")[[1]][1]
  
  myFile
  
  # Read the file and create a data frame
  thisFile =  paste0("/Users/agalianese/Documents/GGCT_BRCA_assocation/bedtools_coverage_out/", startFile)
  myData <- data.frame(read_tsv(file=thisFile, col_names = F, show_col_types = FALSE))

  myData$fileId <- myFile
  ggDF <- rbind(myData, ggDF)

}

ggDF




```


# Quietly read in Clinical Metadata
```{r, message=FALSE}

#TCGA clinical metadata includes family_history.tsv, follow_up.tsv, pathology_detail.tsv > All were empty and removed from the analysis

TCGA_full_clinicalDF <- data.frame(read_delim("~/Documents/GGCT_BRCA_assocation/clinical.cases_selection.2023-10-17/brca_tcga_pan_can_atlas_2018_clinical_data.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
na_TCGA_cols <- data.frame(sapply(TCGA_full_clinicalDF, function(x) all(is.na(x))))
TCGA_clinicalDF <- TCGA_full_clinicalDF[, !na_TCGA_cols, drop = FALSE]

clinicalFullDF <- data.frame(read_tsv("~/Documents/GGCT_BRCA_assocation/clinical.cases_selection.2023-10-17/clinical.tsv", na = c("'--", "not reported")))
na_clinical_cols <- data.frame(sapply(clinicalFullDF, function(x) all(is.na(x))))
clinicalDF <- clinicalFullDF[, !na_clinical_cols, drop = FALSE]

#This will be used to match my data with clinical
matchFullDF <- data.frame(read.delim("~/Downloads/gdc_sample_sheet.2023-10-31.tsv", na = c("'--", "not reported")))
na_match_cols <- data.frame(sapply(matchFullDF, function(x) all(is.na(x))))
matchDF <- matchFullDF[, !na_match_cols, drop = FALSE]
```

#print colnames of all clinical metadata
```{r}
print("TCGA_clinical_data")
head(TCGA_clinical_data)

print("clinicalDF")
head(clinicalDF)

print("matchDF")
head(matchDF)

```

#merge
```{r}

#merge my Sample.ID column into testDF matching using File.ID
myMatchDF <- data.frame(matchDF$File.ID, matchDF$Case.ID, matchDF$Sample.ID)
colnames(myMatchDF) <- c("fileId", "Patient.ID", "Sample.ID")

mergedDF <- merge(mutDF, myMatchDF, by = "fileId", no.dups=F)
myDF <- merge(mergedDF, TCGA_clinicalDF, by="Patient.ID", no.dups=F)

plotDF1 <- merge(ggDF, myMatchDF, by = "fileId", no.dups=F)
plotDF <- merge(plotDF1, TCGA_clinicalDF, by="Patient.ID", no.dups=F)
```


#more SW code
```{r}

finalDF

survfit(Surv())

fit = survfit(Surv(OS, status_OS) ~ ici_cat, data= patient.WT)

Disease.Free.Status


survfit(Surv(OS, status_OS))

ggsurvplot(fit, data=patient.WT, pval=T)

cox.ICI
fit.ICI <- survfit(Surv(OS, status_OS) ~ steroid, data=patient.ICI)
ggsurvplot(fit.ICI, facet.by="ici_cat", data=patient.ICI, pval=T)

cox.ICI <- coxph(Surv(OS, status_OS) ~ age + sex+ EOR+ ici_cat + steroid + KPS, data=patient.ICI)
ggforest(cox.ICI, data=patient.ICI)

dist(myDF$Coverage_Ratio)


plot(myDF$Coverage_Ratio)

head(plotDF)
#X10 = coverage ratio, X4 = exon name
plotDF$Binary.Overall.Survival.Status
ggplot(plotDF, aes(x=X10, y=Months.of.disease.specific.survival, color=Binary.Overall.Survival.Status)) + geom_point()


survfit(Coverage_Ratio ~ Binary.Disease.Free.Status, data = myDF)

cox.ICI <- coxph(Surv(OS, status_OS) ~ age + sex+ EOR+ ici_cat + steroid + KPS, data=patient.ICI)

mutDF


hist(finalDF$Chr)

table(exon2$Coverage_Ratio)
```
#Histograms of Exons
```{r}
hist(as.double(myDF$E1.1.Coverage.Ratio))

myDF

hist(myDF$E1.2.Coverage.Ratio)

hist(myDF$E2.Coverage.Ratio)

hist(myDF$E3.1.Coverage.Ratio)

hist(myDF$E3.2.Coverage.Ratio)

hist(myDF$E4.1.Coverage.Ratio)

hist(myDF$E4.2.Coverage.Ratio)

```


#Survival variables to binary
#survival status is given as a string, instead of the double that R needs
```{r}

table(myDF$Overall.Survival.Status) # 0:LIVING (7119) 1:DECEASED (1400) > n=8519
#myDF$Overall.Survival..Months.

table(myDF$Progression.Free.Status) #0:CENSORED (7294) 1:PROGRESSION (1218) > n=8512
#myDF$Progress.Free.Survival..Months.

table(myDF$Disease.specific.Survival.status) #0:ALIVE OR DEAD TUMOR FREE (7553) 1:DEAD WITH TUMOR (749) > n=8302
#myDF$Months.of.disease.specific.survival

table(myDF$Disease.Free.Status) # 0:DiseaseFree (6629) 1:Recurred/Progressed (651) > n=7280
#myDF$Disease.Free..Months.


myDF <- myDF %>%
  #mutate(Progress.Free.Survival.Years=(Progress.Free.Survival..Months./12)) %>%
  mutate(Binary.Overall.Survival.Status = if_else(str_detect(Overall.Survival.Status, "DECEASED"), 1, 0)) %>%
  mutate(Binary.Progression.Free.Status = if_else(str_detect(Progression.Free.Status, "PROGRESSION"), 1, 0)) %>%
  mutate(Binary.Disease.specific.Survival.status = if_else(str_detect(Disease.specific.Survival.status, "DEAD WITH TUMOR"), 1, 0)) %>%
  mutate(Binary.Disease.Free.Status = if_else(str_detect(Disease.Free.Status, "Recurred/Progressed"), 1, 0)) 

TCGA_clinicalDF <- TCGA_clinicalDF %>%
  #mutate(Progress.Free.Survival.Years=(Progress.Free.Survival..Months./12)) %>%
  mutate(Binary.Overall.Survival.Status = if_else(str_detect(Overall.Survival.Status, "DECEASED"), 1, 0)) %>%
  mutate(Binary.Progression.Free.Status = if_else(str_detect(Progression.Free.Status, "PROGRESSION"), 1, 0)) %>%
  mutate(Binary.Disease.specific.Survival.status = if_else(str_detect(Disease.specific.Survival.status, "DEAD WITH TUMOR"), 1, 0)) %>%
  mutate(Binary.Disease.Free.Status = if_else(str_detect(Disease.Free.Status, "Recurred/Progressed"), 1, 0))

plotDF <- plotDF %>%
  #mutate(Progress.Free.Survival.Years=(Progress.Free.Survival..Months./12)) %>%
  mutate(Binary.Overall.Survival.Status = if_else(str_detect(Overall.Survival.Status, "DECEASED"), 1, 0)) %>%
  mutate(Binary.Progression.Free.Status = if_else(str_detect(Progression.Free.Status, "PROGRESSION"), 1, 0)) %>%
  mutate(Binary.Disease.specific.Survival.status = if_else(str_detect(Disease.specific.Survival.status, "DEAD WITH TUMOR"), 1, 0)) %>%
  mutate(Binary.Disease.Free.Status = if_else(str_detect(Disease.Free.Status, "Recurred/Progressed"), 1, 0))
```

#Plot Associations between Exons and Survivability
```{r}

#E.1 plots

plotDF
ggplot(myDF, aes(x=E1.1.Coverage.Ratio, y=Months.of.disease.specific.survival, color=Binary.Overall.Survival.Status)) + geom_point()

#X4=exon, x10=ratio

ggplot(plotDF, aes(x = X10, y = Binary.Overall.Survival.Status, color = X4)) +
  geom_point() + stat_compare_means() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plotDF, aes(x = Binary.Overall.Survival.Status, y = X10, color = X4)) +
  geom_violin() + stat_compare_means() 

ggplot(plotDF, aes(x = X4, y = A_Bases, color = Exon)) +
  geom_violin() + stat_compare_means() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot(plotDF, aes(x = X4, y = A_Len, color = Exon)) +
  geom_boxplot() + stat_compare_means() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



lows <- myDF %>%
  filter(Coverage_Ratio < 1) %>%
  arrange(Coverage_Ratio)

hist(lows$Coverage_Ratio)

ggplot(lows, aes(x = Exon, y = Coverage_Ratio, color = Exon)) +
  geom_boxplot() + stat_compare_means() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

```

```{r}
out1 <- coxph(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ Diagnosis.Age + Fraction.Genome.Altered + Sex, data=myDF)

out2 <- coxph(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ E1.1.Coverage.Ratio, data=myDF)

E2.Coverage.Ratio + E3.1.Coverage.Ratio + E3.2.Coverage.Ratio + E4.1.Coverage.Ratio + E4.2.Coverage.Ratio

#patients with high GGCT expression had a poor prognosis.
out1 <- coxph(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ Diagnosis.Age + Fraction.Genome.Altered + Sex, data=myDF)

#
ggforest(coxph(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ E1.1.Coverage.Ratio + E1.2.Coverage.Ratio + E2.Coverage.Ratio + E3.1.Coverage.Ratio + E3.2.Coverage.Ratio + E4.1.Coverage.Ratio + E4.2.Coverage.Ratio
, data=myDF))

E2.Coverage.Ratio + E3.1.Coverage.Ratio + E3.2.Coverage.Ratio + E4.1.Coverage.Ratio + E4.2.Coverage.Rat



Li et al. also reported that GGCT was upregulated in ovarian cancers and associated with advanced FIGO (International Federation of Gynecology and Obstetrics) stage, lymph node metastases, and ascitic fluid volume in high-grade serous ovarian cancers (HGSCs)

ggforest(out1, data=myDF)
ggforest(out2, data=myDF)

myDF
```


```{r}
myDF$Sex

table(myDF$Overall.Survival.Status) # 0:LIVING (7119) 1:DECEASED (1400) > n=8519
#myDF$Overall.Survival..Months.
TCGA_clinicalDF

fit <- survfit(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ Sex, data = TCGA_clinicalDF)



ggsurvplot(fit, data = TCGA_clinicalDF,
           title = "Survival Curves",
           pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "Sex",               # Change legend titles
           legend.labs = c("Male", "female"),  # Change legend labels
           palette = "jco",                    # Use JCO journal color palette
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
)

```
 

```{r}

fit <- survfit(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ Sex, data = TCGA_clinicalDF)


# Drawing curves
ggsurvplot(fit, data = TCGA_clinicalDF, pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
          surv.median.line = "hv") +           # Add median survival lines
          labs(x = "Time (in Months)", y = "Overall survival probability", title="Title")               # Change legend titles)

survfit2(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ Sex, data = TCGA_clinicalDF) %>% 
  ggsurvfit() +
  labs(x = "Time (in Months)", y = "Overall survival probability") + 
  add_confidence_interval() +
  add_risktable()
```

#Hypothesis from 
```{r}

#patients with high GGCT expression had a poor prognosis.

hist(plotDF$X10)

out1 <- coxph(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ Diagnosis.Age + Fraction.Genome.Altered + Sex, data=myDF)

out2 <- coxph(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ E1.1.Coverage.Ratio + E1.2.Coverage.Ratio + E2.Coverage.Ratio + E3.1.Coverage.Ratio + E3.2.Coverage.Ratio + E4.1.Coverage.Ratio + E4.2.Coverage.Ratio, data=myDF)

E2.Coverage.Ratio + E3.1.Coverage.Ratio + E3.2.Coverage.Ratio + E4.1.Coverage.Ratio + E4.2.Coverage.Ratio
ggtree(out2)


Li et al. also reported that GGCT was upregulated in ovarian cancers and associated with advanced FIGO (International Federation of Gynecology and Obstetrics) stage, lymph node metastases, and ascitic fluid volume in high-grade serous ovarian cancers (HGSCs)

fit <- survfit(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ Sex + Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code, data = TCGA_clinicalDF)

hist(myDF$Aneuploidy.Score)
hist(myDF$TMB..nonsynonymous.)
hist(myDF$Primary.Lymph.Node.Presentation.Assessment)
table(myDF$ICD.10.Classification)


# Drawing curves
ggsurvplot(fit, data = TCGA_clinicalDF, pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
          surv.median.line = "hv") +           # Add median survival lines
          labs(x = "Time (in Months)", y = "Overall survival probability", title="Title")               # Change legend titles)

survfit2(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ E1.1.Coverage.Ratio, data = myDF) %>% 
  ggsurvfit() +
  labs(x = "Time (in Months)", y = "Overall survival probability") + 
  add_confidence_interval() +
  add_risktable()


Their multivariate analysis showed that FIGO stage, lymph node metastasis and GGCT expression were independent prognostic factors for overall and progression free survival

showed that the region located at 371 to +14 bp of the 5 0 end of GGCT is important for the activation of its transcription in both cancer (HeLa, MCF7) and non-cancer cells (IMR-90).

high expression of GGCT was associated with poor survival.
```

