---
title: "R Notebook"
output: html_notebook
---

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
  
  if (myData[3,10] != myData[4,10]) {
    print(myFile)
  }
  
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

#Calculate isoform ratios
                                
```{r}
myDF <- myDF %>%
  mutate(
    
    Isoform.Ratio.202 = as.double(E3.1.N.Features) / (as.double(E3.1.N.Features) + as.double(E3.2.N.Features)),
    Isoform.Ratio.204 = as.double(E3.2.N.Features) / (as.double(E3.1.N.Features) + as.double(E3.2.N.Features)),
    
    Main.Isoform = ifelse(Isoform.Expression.202 > Isoform.Expression.204, "Isoform.202", "Isoform.204"),
  )
```
