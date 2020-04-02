# This is a wrapper script for the main analysis. In order to run the CellTypeMPPvalidation script, the region specific data needs to
# be accessed through https://www.synapse.org/#!Synapse:syn5613797

#Basic analysis in the PW cohort 
source("ProjectScripts/ConfigFilePV.R") #This files contains all the relevant file locations and filters samples relevant to PW cohort
source("ProjectScripts/AnalysisCaseControl.R")

#Basic analysis in the PW cohort
source("ProjectScripts/ConfigFilePW.R") #This files contains all the relevant file locations and filters samples relevant to NBB cohort
source("ProjectScripts/AnalysisCaseControl.R")


# Comparing both cohorts
source("ProjectScripts/ComparePVNBB.R")  #this is the main analysis for comparing the two cohorts

source("ProjectScripts/MoreComparison.R") #this is mostly for the MA and manhattan plots. Also has the statistics for the overlaps
source("ProjectScripts/AdditionalAnalyses.R") #this part does enrichment analysis based on decoupling and performs analyses related to
                                              #association with p300 and HDAC binding sites