outie <- glm(Binary.Disease.Free.Status ~ Isoform.Ratio.202 + Diagnosis.Age + Fraction.Genome.Altered + Majorly.Expressed.Isoform + Subtype + E3.202.Features + Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code, data=myDF, family = binomial(link = "logit"))
summary(outie)
hist(outie$residuals)
plot(outie$residuals)


#full backwards regression
