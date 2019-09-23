# ############ For PV RNA samples only ########################
# dataDir = paste0(getwd(), "/data/")
# ResultsPath = "ResultsPV_Final"
# SampleFilter =  function(x) x %>% filter(!is.na(sample_id_rna), cohort == "PV") %>% droplevels()
 
# CasePeakLoc = "data/Peaks/Case_PV_RNA.narrowPeakClean.gz"
# ControlPeakLoc = "data/Peaks/Control_PV_RNA.narrowPeakClean.gz"
# CombinedPeakLoc = "data/Peaks/PV_ALL_RNA.narrowPeakClean.gz"
 
# CountMatrixLoc = "data/PV_ALL_RNA_counts.tsv.gz"
# CellTypePeakCountLoc = "CellTypeH3K27ac/data/PV_ALL_RNA_counts.tsv"
# Cohort = "PV"

############ For NBB samples only ########################
dataDir = paste0(getwd(), "/data/")
ResultsPath = "ResultsNBB_Final"
SampleFilter =  function(x) x %>% filter(cohort == "NBB") %>% droplevels()

CasePeakLoc = "data/Peaks/Case_NBB.narrowPeakClean.gz"
ControlPeakLoc = "data/Peaks/Control_NBB.narrowPeakClean.gz"
CombinedPeakLoc = "data/Peaks/NBB_ALL.narrowPeakClean.gz"

CountMatrixLoc = "data/NBB_ALL_counts.tsv.gz"
CellTypePeakCountLoc = "CellTypeH3K27ac/data/NBB_ALL_counts.tsv"
Cohort = "NBB"

############ Aging ########################
# dataDir = paste0(getwd(), "/data/")
# ResultsPath = "ResultsAging"
# SampleFilter =  function(x) x %>% filter(condition == "Cont", cohort == "PV", age > 20) %>% droplevels()
# 
# ControlPeakLoc = paste0(dataDir, "Peaks/Control_PV_above_20.narrowPeakClean.gz")
# CombinedPeakLoc = paste0(dataDir, "Peaks/Control_PV_above_20.narrowPeakClean.gz")
# 
# CountMatrixLoc = paste0(dataDir, "PV_ALL_above_20_counts.tsv.gz")
# CellTypePeakCountLoc = "CellTypeH3K27ac/data/PV_ALL_above_20_counts.tsv"
# Cohort = "AgingControl"


############# All adult parkome samples ########################
# dataDir = paste0(getwd(), "/data/")
# ResultsPath = "ResultsAllPVadult"
# SampleFilter =  function(x) x %>% filter(age > 41, cohort == "PV") %>% droplevels()
# 
# CasePeakLoc = "data/Peaks/Case_PV_above_20.narrowPeakClean.gz"
# ControlPeakLoc = "data/Peaks/Control_PV_above_20.narrowPeakClean.gz"
# CombinedPeakLoc = "data/Peaks/PV_ALL_above_20.narrowPeakClean.gz"
# 
# CountMatrixLoc = "data/PV_ALL_above_20_counts.tsv.gz"
# CellTypePeakCountLoc = "CellTypeH3K27ac/data/PV_ALL_above_20_counts.tsv"
# Cohort = "PV all adults"

############# All adult parkome samples (including 20-40) ########################
# dataDir = paste0(getwd(), "/data/")
# ResultsPath = "ResultsPVabove20"
# SampleFilter =  function(x) x %>% filter(age > 19, cohort == "PV") %>% droplevels()
# 
# CasePeakLoc = "data/Peaks/Case_PV_above_20.narrowPeakClean.gz"
# ControlPeakLoc = "data/Peaks/Control_PV_above_20.narrowPeakClean.gz"
# CombinedPeakLoc = "data/Peaks/PV_ALL_above_20.narrowPeakClean.gz"
# 
# CountMatrixLoc = "data/PV_ALL_above_20_counts.tsv.gz"
# CellTypePeakCountLoc = "CellTypeH3K27ac/data/PV_ALL_above_20_counts.tsv"
# Cohort = "PV all adults above 20"
