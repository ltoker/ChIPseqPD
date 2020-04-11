############ For PV RNA samples only ########################
dataDir = paste0(getwd(), "/data/")
ResultsPath = "ResultsPV_Final"
SampleFilter =  function(x) x %>% filter(!is.na(sample_id_rna), cohort == "PV") %>% droplevels()
 
CasePeakLoc = "data/Peaks/Case_PV_RNA.narrowPeakClean.gz"
ControlPeakLoc = "data/Peaks/Control_PV_RNA.narrowPeakClean.gz"
CombinedPeakLoc = "data/Peaks/PV_ALL_RNA.narrowPeakClean.gz"
 
CountMatrixLoc = "data/PV_ALL_RNA_counts.tsv.gz"
CellTypePeakCountLoc = "CellTypeH3K27ac/data/PV_ALL_RNA_counts.tsv"
Cohort = "PV"


