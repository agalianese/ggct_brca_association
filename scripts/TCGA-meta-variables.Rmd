# Continuous Variables
Diagnosis.Age 
Aneuploidy.Score	
Buffa.Hypoxia.Score	| Ragnum.Hypoxia.Score	| Winter.Hypoxia.Score	
MSI.MANTIS.Score	
MSIsensor.Score	
Mutation.Count	
Fraction.Genome.Altered
TMB..nonsynonymous.	

# Categorical Variables
Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code	> 12 cats, n=1174	
Ethnicity.Category | Genetic.Ancestry.Label | Race.Category		
Oncotree.Code	| Tumor.Type	|	ICD.10.Classification | Cancer.Type.Detailed > only using cancer.type.detailed
American.Joint.Committee.on.Cancer.Metastasis.Stage.Code > 4 cats, n=1174	
Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code > cat = 16, n=1168	
American.Joint.Committee.on.Cancer.Tumor.Stage.Code	> 13 cats, n=1174
Prior.Diagnosis	#no=1108, yes=65
Sex	

# Maybe
Neoadjuvant.Therapy.Type.Administered.Prior.To.Resection.Text	 > no=1166, yes=7
New.Neoplasm.Event.Post.Initial.Therapy.Indicator > no=833, yes=103, n=936
Subtype	> 4 cats, n=1066 but great data
Radiation.Therapy	#no=450, yes=562, n=1042
Number.of.Samples.Per.Patient	[1,4]
Person.Neoplasm.Cancer.Status	#tumor_free=941, with_tumor=103, n=1044

# Outcome Variables
Overall.Survival..Months.	+ Binary.Overall.Survival.Status
Binary.Progression.Free.Status	+ Progress.Free.Survival..Months.	
Binary.Disease.specific.Survival.status + Months.of.disease.specific.survival	
Binary.Disease.Free.Status + Disease.Free..Months.

# Discarded Variables
Birth.from.Initial.Pathologic.Diagnosis.Date	
Last.Alive.Less.Initial.Pathologic.Diagnosis.Date.Calculated.Day.Value	
Sample.Type	> all="Primary"
Somatic.Status> all ="Matched"
Tissue.Source.Site.Code > where the tissue is from ex: Columbia, Texas Anderson etc.
Other.Patient.ID	
Primary.Lymph.Node.Presentation.Assessment > only 794 cases, removed for NA
Tissue.Prospective.Collection.Indicator	| Tissue.Retrospective.Collection.Indicator	 > whether the tissue was 
American.Joint.Committee.on.Cancer.Publication.Version.Type > n=1002, too many missing vars
