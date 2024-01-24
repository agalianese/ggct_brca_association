---
title: "R Notebook"
output: html_notebook
---


```{r}
library(survival)
library(survminer)
library(ggpmisc)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(survival); library(survminer)

setwd("~/Documents/GGCT_BRCA_assocation")

```

#read in bed coverage output
```{r, quiet=T}
file_list <- list.files(path="~/Documents/GGCT_BRCA_assocation/bedtools_coverage_out", pattern = "*.ggct.bed", full.names = F)

#read in bedfile, mutate and add extra column for name of bed file
mutDF0 = data.frame(matrix(nrow=1231, ncol=3))

problematic_files <- c()

for (startFile in file_list) { 
   
  #remove .ggct.bed to get the file nam
  myFile=str_split(startFile, ".ggct.bed")[[1]][1]
    
  # Read the file and create a data frame
  thisFile =  paste0("/Users/agalianese/Documents/GGCT_BRCA_assocation/bedtools_coverage_out/", startFile)
    
    myData <- data.frame(read_tsv(file=thisFile, col_names = F, show_col_types = FALSE))
  
    tryCatch({
      # Check for parsing issues
      if (any(problems(thisFile))) {
        print(paste0("Parsing issues in file:", thisFile, "\n"))
        print(problems(thisFile))
        problematic_files <- c(problematic_files, thisFile)
      
      } else {
        # Continue processing if no parsing issues
        i = which(str_starts(file_list, myFile))
        mutDF0[i, 1] = as.double(myData[4, 7])
        mutDF0[i, 2] = as.double(myData[3, 7])
        mutDF0[i, 3] = myFile
      }
    }, error = function(e) {
      cat("Error reading file:", thisFile, "\n")
      problematic_files <- c(problematic_files, thisFile)
    })
} 

print(problematic_files)
```


```{r, quiet=T}
colnames(mutDF0) <- c("E3.202.Features","E3.204.Features", "fileId")

mutDF <- mutDF0 %>%
  filter(!(E3.202.Features==0 & E3.204.Features == 0)) %>%
  mutate(
    Isoform.Ratio.202 = as.double(E3.202.Features) / (as.double(E3.202.Features) + as.double(E3.204.Features)),
    Isoform.Ratio.204 = as.double(E3.204.Features) / (as.double(E3.202.Features) + as.double(E3.204.Features)),
    
    Majorly.Expressed.Isoform = ifelse(Isoform.Ratio.202 > Isoform.Ratio.204, "Isoform.202", "Isoform.204"),
  )
mutDF

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
print("TCGA_clinicalDF")
head(TCGA_clinicalDF)

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
```

```{r}

myDF %>% filter(Majorly.Expressed.Isoform == "Isoform.204")

write_csv(myDF, file="~/Desktop/mine.csv")
```
