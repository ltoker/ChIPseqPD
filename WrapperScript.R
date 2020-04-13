# This is a wrapper script for the main analysis. In order to run the CellTypeMPPvalidation script, the region specific data needs to
# be accessed through https://www.synapse.org/#!Synapse:syn5613797. Of note, enrichment analysis performed as part of the AdditionalAnalyses.R script,
# requires Java and might not run on sertain Java vrsions on windows machine. If error occurs, this part can be run on Linux or outside of R on 
# ermineJ software https://erminej.msl.ubc.ca/

#Comparison of ChIP-seq normalization methods
source("ProjectScripts/NormalizationComparison.R")
rm(list = ls(all.names = TRUE))


#Basic analysis in the PW cohort 
source("ProjectScripts/ConfigFilePV.R") #This files contains all the relevant file locations and filters samples relevant to PW cohort
source("ProjectScripts/AnalysisCaseControl.R")
rm(list = ls(all.names = TRUE))

#Basic analysis in the NBB cohort
source("ProjectScripts/ConfigFileNBB.R") #This files contains all the relevant file locations and filters samples relevant to NBB cohort
source("ProjectScripts/AnalysisCaseControl.R")
rm(list = ls(all.names = TRUE))


# Comparing both cohorts
source("ProjectScripts/ComparePVNBB.R")  #this is the main analysis for comparing the two cohorts
rm(list = ls(all.names = TRUE))

source("ProjectScripts/MoreComparison.R") #This is for the MA and manhattan plots.
rm(list = ls(all.names = TRUE))


source("ProjectScripts/AdditionalAnalyses.R") #this part does enrichment analysis based on decoupling and performs analyses related to
                                              #association with p300 and HDAC binding site.
rm(list = ls(all.names = TRUE))

source("ProjectScripts/WetLabPlotting.R") #this parts plots all the wet lab-related figures
