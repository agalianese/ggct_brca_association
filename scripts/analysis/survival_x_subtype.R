---
title: "R Notebook"
output: html_notebook
---

# DSS x ILC
```{r}

ilc <- myDF %>%
    filter(Oncotree.Code == "ILC")

idc <- myDF %>%
    filter(Oncotree.Code == "IDC")


subtypes <- readxl::read_xlsx("/ag4555/RStudio/code/myGGCT/myGGCT/ILC_subtypes.xlsx", col_names = T, skip=1)
colnames(subtypes) <- c("oldSampleId", "ILC_Subtype", "Patient.ID")

mySubs <- left_join(x=ilc, y=subtypes, by="Patient.ID") %>% filter(!is.na(ILC_Subtype))

```

# Disease Free Survival by Subtype 
```{r}
selected_columns <- c("Subtype", "Disease.Free..Months.", "Months.of.disease.specific.survival", "Overall.Survival..Months.", "Progress.Free.Survival..Months.")

selected_data <- dplyr::select(myDF, all_of(selected_columns))
longSurvival <- reshape2::melt(selected_data, id.vars="Subtype")

#add anova to plot
good <- ggplot(longSurvival, aes(x = Subtype, y = value, group_by=variable, fill=Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Survival in Months", title = "Survival against Subtype") +
  theme_minimal() +  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) + facet_wrap("variable", nrow=1) 

out <- good + stat_anova_test(aes(group = Subtype)) 
out

#ggsave(file="alive_overall.png", out, height=10, width=8)
```


#Dead Disease Specific
```{r}

selectedDead <- dplyr::select(deadDiseaseSpecific, all_of(selected_columns))
deadSubtype <- reshape2::melt(selectedDead, id.vars="Subtype")

deadSubtype
deadPlot <- ggplot(deadSubtype, aes(x = Subtype, y = value, group_by=variable, fill=Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Survival in Months", title = "Survival against Subtype") +
  theme_minimal() + facet_wrap("variable", nrow=1) +  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))


deadPlot <- ggplot(deadSubtype, aes(x = Subtype, y = value, group_by=variable, fill=Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Survival in Months", title = "Survival against Subtype") +
  theme_minimal() + facet_wrap("variable", nrow=1) +  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))


good <- deadPlot  + stat_kruskal_test(aes(group = Subtype))
good

#ggsave(file="dead_DS.png", good, height=10, width=10)

```



#Dead Disease Free
```{r}

selectedDead <- dplyr::select(deadDiseaseFree, all_of(selected_columns))
deadSubtype <- reshape2::melt(selectedDead, id.vars="Subtype")

deadPlot <- ggplot(deadSubtype, aes(x = Subtype, y = value, group_by=variable, fill=Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Survival in Months", title = "Survival against Subtype") +
  theme_minimal() + facet_wrap("variable", nrow=2) +  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

ready <- deadPlot + stat_kruskal_test(aes(group = Subtype))
ready
#ggsave(file="/ag4555/RStudio/code/myGGCT/myGGCT/dead_DF.png", ready, height=10, width=8)

```


#Dead Progression 
```{r}
selectedDead <- dplyr::select(deadProgressionSpecific, all_of(selected_columns))
deadSubtype <- reshape2::melt(selectedDead, id.vars="Subtype")

deadPlot <- ggplot(deadSubtype, aes(x = Subtype, y = value, group_by=variable, fill=Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Survival in Months", title = "Survival against Subtype") +
  theme_minimal() + facet_wrap("variable", nrow=1) +  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

final <- deadPlot + stat_kruskal_test(aes(group = Subtype))
final
```


```{r}
ggplot(myDF, aes(x=Subtype, y=Isoform.Ratio.202, fill=Subtype)) + geom_boxplot() + stat_anova_test(aes(group = Subtype)) + scale_y_continuous(trans = "log10")


ggplot(myDF, aes(x=Subtype, y=Isoform.Ratio.202, fill=Subtype)) + geom_boxplot() + stat_anova_test(aes(group = Subtype))

ggplot(myDF, aes(x=Subtype, y=Isoform.Ratio.202, fill=Subtype)) + geom_boxplot() + stat_anova_test(aes(group = Subtype))

myDF


#ggsave(file="/ag4555/RStudio/code/myGGCT/myGGCT/dead_PF.png", final, height=10, width=8)
```


```{r}
table(ilc$Subtype)
table(idc$Subtype)
```


# KM - subtype
```{r}

ggplot(myDF, aes(x=Subtype, y=Months.of.disease.specific.survival, fill=Subtype)) + geom_boxplot() + theme_pubr() + stat_anova_test(aes(group=Subtype))

fit = survfit(Surv(Months.of.disease.specific.survival, Binary.Disease.specific.Survival.status) ~ Subtype, data= myDF)
test <- ggsurvplot(fit, data=myDF, pval=T, risk.table = T, surv.median.line="hv")


test
#ggsave(file="/ag4555/RStudio/code/myGGCT/myGGCT/plots/KM_alive_subtype.png", plot=(test), height=10, width=10)
```


```{r}

deadOverall

ggplot(deadOverall, aes(x=Subtype, y=Months.of.disease.specific.survival, fill=Subtype)) + geom_boxplot() + theme_pubr() + stat_anova_test(aes(group=Subtype))

fit = survfit(Surv(Months.of.disease.specific.survival, Binary.Disease.specific.Survival.status) ~ Subtype, data= deadOverall)
test <- ggsurvplot(fit, data=deadOverall, pval=T, risk.table = T, surv.median.line="hv")

#ggsave(file="/ag4555/RStudio/code/myGGCT/myGGCT/plots/KM_dead_subtypes.png", (test), height=10, width=10)
```



#Dead Overall
```{r}

selectedDead <- dplyr::select(deadOverall, all_of(selected_columns))
deadSubtype <- reshape2::melt(selectedDead, id.vars="Subtype")

deadPlot <- ggplot(deadSubtype, aes(x = Subtype, y = value, group_by=variable, fill=Subtype)) +
  geom_boxplot() +
  labs(x = "Subtype", y = "Survival in Months", title = "Survival against Subtype") +
  facet_wrap("variable") +  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

final <- deadPlot + stat_kruskal_test(aes(group = Subtype))

final

#ggsave(file="/ag4555/RStudio/code/myGGCT/myGGCT/plots/dead_overall.png", final, height=10, width=10)

```





#we determined that reactive-like ILC patients had a significantly better disease-specific (DSS) (p = 0.038, HR: 0.47) and overall survival (OS) (p = 0.023, HR: 0.50) compared to proliferative ILC patients in the METABRIC dataset, which has a median follow-up of 7.2 years (compared to the TCGA median follow-up of less than 2 years
```{r}


ggplot(mySubs, aes(x=ILC_Subtype, y=Months.of.disease.specific.survival, fill=ILC_Subtype)) + geom_boxplot() + theme_pubr() + stat_kruskal_test()

fit = survfit(Surv(Months.of.disease.specific.survival, Binary.Disease.specific.Survival.status) ~ ILC_Subtype, data= mySubs)
ggsurvplot(fit, data=mySubs, pval=T)
```


```{r}
ggplot(mySubs, aes(x=ILC_Subtype, y=Overall.Survival..Months., fill=ILC_Subtype)) + geom_boxplot() + theme_pubr()  + stat_kruskal_test()

fit = survfit(Surv(Overall.Survival..Months., Binary.Overall.Survival.Status) ~ ILC_Subtype, data= mySubs)
ggsurvplot(fit, data=mySubs, pval=T)
```


```{r}
ggplot(mySubs, aes(x=ILC_Subtype, y=Progress.Free.Survival..Months., fill=ILC_Subtype)) + geom_boxplot() + theme_pubr() + stat_kruskal_test()

fit = survfit(Surv(Progress.Free.Survival..Months., Binary.Progression.Free.Status) ~ ILC_Subtype, data= mySubs)
ggsurvplot(fit, data=mySubs, pval=T)

```


```{r}
fit = survfit(Surv(Disease.Free..Months., Binary.Disease.Free.Status) ~ ILC_Subtype, data= mySubs)
sig <- ggsurvplot(fit, data=mySubs, pval=T,risk.table = T)

sig


#ggsave("/ag4555/RStudio/code/myGGCT/myGGCT/plots/survplot.png", plot = print(sig), height=10, width=8)
```


```{r}
cox.fit <- coxph(Surv(Disease.Free..Months., Binary.Disease.Free.Status) ~ ILC_Subtype, data= mySubs)

ggforest(cox.fit, data=mySubs)

```


```{r}
deadDiseaseSpecific <- disease.specific.df %>%
    filter(Binary.Disease.specific.Survival.status == 1)

nrow(deadDiseaseSpecific)

#fit logistic regression model
selected_columns <- c("ILC_Subtype", "Disease.Free..Months.", "Months.of.disease.specific.survival", "Overall.Survival..Months.", "Progress.Free.Survival..Months.")
littleSubs <- dplyr::select(mySubs, all_of(selected_columns))

fit = survfit(Surv(mySubs$Binary.Disease.specific.Survival.status, mySubs$Months.of.disease.specific.survival) ~ ILC_Subtype, data= mySubs)



plot(fit)
ggsurvplot(fit, data=mySubs, pval=T)

```



#No sig-nificant differences in DSS or OS were identified between the immune-related subgroup and either the proliferative or reactive-like subgroup. These results are consistent with previ-ous studies reporting that the reactive stromal phenotype is associated with a good prognosis in breast cancer while prolifer-ation is one of the strongest indicators of worse outcome in luminal/ER+ breast cancers
```{r}

```

