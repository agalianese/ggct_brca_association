---
title: "Dead patients only"
output: html_notebook
---


```{r}

deadOverall <- overall.df %>%
    filter(Binary.Overall.Survival.Status == 1)

deadDiseaseFree <- disease.free.df %>%
    filter(Binary.Disease.Free.Status == 1)

deadDiseaseSpecific <- disease.specific.df %>%
    filter(Binary.Disease.specific.Survival.status == 1)

deadProgressionSpecific <- progression.free.df %>%
    filter(Binary.Progression.Free.Status == 1)

```



# Dead Overall  - Isoform.Ratio.202
```{r}
nrow(deadOverall)

#fit logistic regression model
model <- glm(cbind(deadOverall$Overall.Survival..Months., deadOverall$Binary.Overall.Survival.Status) ~ Isoform.Ratio.202, family = binomial(link = "logit"), data=deadOverall)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))

```


# Dead Disease Free - Isoform.Ratio.202 
```{r}

nrow(deadDiseaseFree)

#fit logistic regression model
model <- glm(cbind(deadDiseaseFree$Disease.Free..Months., deadDiseaseFree$Binary.Disease.Free.Status) ~ Isoform.Ratio.202, family = binomial(link = "logit"), data=deadDiseaseFree)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))
```



# Dead Isoform.Ratio.202 - Disease Specific
```{r}



nrow(deadDiseaseSpecific)

#fit logistic regression model
model <- glm(cbind(deadDiseaseSpecific$Months.of.disease.specific.survival, deadDiseaseSpecific$Binary.Disease.specific.Survival.status) ~ Isoform.Ratio.202, family = binomial(link = "logit"), data=deadDiseaseSpecific)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))
```


# Dead Isoform.Ratio.202 - Progression Free
```{r}




nrow(deadProgressionSpecific)

#fit logistic regression model
model <- glm(cbind(deadProgressionSpecific$Progress.Free.Survival..Months., deadProgressionSpecific$Binary.Progression.Free.Status) ~ Isoform.Ratio.202, family = binomial(link = "logit"), data=deadProgressionSpecific)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))
```





# Dead - Overall
```{r}

nrow(deadOverall)

#fit logistic regression model
model <- glm(cbind(deadOverall$Overall.Survival..Months., deadOverall$Binary.Overall.Survival.Status) ~ Isoform.Ratio.202 + E3.202.Features + E3.204.Features, family = binomial(link = "logit"), data=deadOverall)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))

```



# Dead - Disease Free
```{r}
nrow(deadDiseaseFree)

#fit logistic regression model
model <- glm(cbind(deadDiseaseFree$Disease.Free..Months., deadDiseaseFree$Binary.Disease.Free.Status) ~ Isoform.Ratio.202 + E3.202.Features + E3.204.Features, family = binomial(link = "logit"), data=deadDiseaseFree)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))
```


# Dead People - Disease Specific
```{r}

nrow(deadDiseaseSpecific)

#fit logistic regression model
model <- glm(cbind(deadDiseaseSpecific$Months.of.disease.specific.survival, deadDiseaseSpecific$Binary.Disease.specific.Survival.status) ~ Isoform.Ratio.202 + E3.202.Features + E3.204.Features, family = binomial(link = "logit"), data=deadDiseaseSpecific)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))
```


# Dead People - Progression Free
```{r}




nrow(deadProgressionSpecific)

#fit logistic regression model
model <- glm(cbind(deadProgressionSpecific$Progress.Free.Survival..Months., deadProgressionSpecific$Binary.Progression.Free.Status) ~ Isoform.Ratio.202 + E3.202.Features + E3.204.Features, family = binomial(link = "logit"), data=deadProgressionSpecific)

#disable scientific notation for model summary
options(scipen=999)

#view model summary
summary(model)

#no predictive power
pscl::pR2(model)["McFadden"]

exp(coef(model))
```

