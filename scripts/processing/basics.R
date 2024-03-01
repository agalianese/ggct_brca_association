
```{r}
.libPaths("/R")

library(survival)
library(survminer)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggpmisc) 
library(MASS)
library(GGally)
library(ggrepel)
library(car)
library(pscl)
library(ggforce)

#myDF <- read_csv("/home/myGGCT/mine.csv")
myDF <- read_csv("/ag4555/RStudio/code/myGGCT/mine.csv", show_col_types = FALSE)

# replace subtype BRCA_Normal with NA because Normal isn't a real subtype
myDF <- myDF %>%
    mutate(Subtype= ifelse(Subtype == "BRCA_Normal", NA, Subtype))

```


```{r}
table(myDF$Subtype)
table(is.na(myDF$Subtype))

```


# method so i can ggsave my survplots
```{r}
# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
```


#gender piechart
```{r}

data <- data.frame(
  gender = c("Female", "Male"),
  count = c(1161, 13)
)

data <- data %>% mutate(percent = count*100/(1161+13))

# Create the pie chart
pie_chart <- ggplot(data, aes(x = "", y = count, fill = gender)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(label = paste(gender, ": ", round(percent, 2), "%")), position = position_stack(vjust = 0.5))

pie_chart

ggsave(file="/ag4555/RStudio/code/myGGCT/myGGCT/plots/gender_pie.png", (pie_chart), height=10, width=10)
```

#age hist
```{r}

age <- ggplot(myDF, aes(x=Diagnosis.Age, fill=Diagnosis.Age)) + geom_bar() +  theme_minimal() 

age

ggsave(file="/ag4555/RStudio/code/myGGCT/myGGCT/plots/age_hist.png", (age), height=10, width=10)
```


#overall surv numbers
```{r}
selected_columns <- c("Binary.Disease.Free.Status", "Binary.Disease.specific.Survival.status", "Binary.Progression.Free.Status", "Binary.Overall.Survival.Status")
selected_data <- dplyr::select(myDF, all_of(selected_columns))
```

```{r}
table(myDF$Binary.Disease.Free.Status)
table(myDF$Binary.Disease.specific.Survival.status)
table(myDF$Binary.Progression.Free.Status)
table(myDF$Binary.Overall.Survival.Status)
```



#Survival by Subtype 
```{r}
head(myDF)
table(myDF$Oncotree.Code)
```


```{r}
#binary.disease.free is missing the most amount of data (about 160 cases when compared w overall survival)

disease.free.df <- myDF[complete.cases(myDF$Binary.Disease.Free.Status), ]
disease.specific.df <- myDF[complete.cases(myDF$Binary.Disease.specific.Survival.status), ]
progression.free.df <- myDF[complete.cases(myDF$Binary.Progression.Free.Status), ]
overall.df <- myDF[complete.cases(myDF$Binary.Overall.Survival.Status), ]

nrow(disease.free.df)
nrow(disease.specific.df)
nrow(progression.free.df)
nrow(overall.df)
```




# Categorical Variables
```{r}
table(myDF$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code) #12 cats, n=1174
table(myDF$American.Joint.Committee.on.Cancer.Metastasis.Stage.Code) #12 cats, n=1174
table(myDF$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code) 
table(myDF$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
```



```{r}
#binary.disease.free is missing the most amount of data (about 160 cases when compared w overall survival)

disease.free.df <- myDF[complete.cases(myDF$Binary.Disease.Free.Status), ]
disease.specific.df <- myDF[complete.cases(myDF$Binary.Disease.specific.Survival.status), ]
progression.free.df <- myDF[complete.cases(myDF$Binary.Progression.Free.Status), ]
overall.df <- myDF[complete.cases(myDF$Binary.Overall.Survival.Status), ]

nrow(disease.free.df)
nrow(disease.specific.df)
nrow(progression.free.df)
nrow(overall.df)
```
