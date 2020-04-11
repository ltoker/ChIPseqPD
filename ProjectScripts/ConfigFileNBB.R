########### For NBB samples only ########################
dataDir = paste0(getwd(), "/data/")
ResultsPath = "ResultsNBB_Final"
SampleFilter =  function(x) x %>% filter(cohort == "NBB") %>% droplevels()
 
CasePeakLoc = "data/Peaks/Case_NBB.narrowPeakClean.gz"
ControlPeakLoc = "data/Peaks/Control_NBB.narrowPeakClean.gz"
CombinedPeakLoc = "data/Peaks/NBB_ALL.narrowPeakClean.gz"
 
CountMatrixLoc = "data/NBB_ALL_counts.tsv.gz"
CellTypePeakCountLoc = "CellTypeH3K27ac/data/NBB_ALL_counts.tsv"
Cohort = "NBB"
