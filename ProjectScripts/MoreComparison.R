source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
ResultsPath = "GeneralResults"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")
packageF("qqman")


NBBresults <- readRDS("ResultsNBB_Final//NBBDEresults.Rds") %>% arrange(padj)
PVresults <- readRDS("ResultsPV_Final//PVDEresults.Rds") %>% arrange(padj) 

ggMA(resultObject = NBBresults %>% filter(!duplicated(PeakName)), lims = c(-1.4,1.4), geneColName = PeakName, contours = T, colour_pvalue = T)
ggsave(paste0(ResultsPath,"NBB_MAplot.pdf"), device = "pdf", width = 10, height = 4, dpi = 300, useDingbats = F)

ggMA(resultObject = PVresults %>% filter(!duplicated(PeakName)), lims = c(-1,1), geneColName = PeakName, contours = T, colour_pvalue = T)
ggsave(paste0(ResultsPath,"PV_MAplot.pdf"), device = "pdf", width = 10, height = 4, dpi = 300, useDingbats = F)

#Identification of genes with DARs in the same direction
TotalGenesNBB <- NBBresults %>% .[!grepl("^MIR|SNOR", .$symbol),] %>%
  filter(!is.na(symbol)) %>% .$symbol %>% unique
TotalGenesNBBup <- NBBresults %>% filter(padj < 0.05, log2FoldChange > 0) %>%
  .[!grepl("^MIR|SNOR", .$symbol),] %>% filter(!is.na(symbol)) %>% .$symbol %>% unique
TotalGenesNBBdown <- NBBresults %>% filter(padj < 0.05, log2FoldChange < 0) %>%
  .[!grepl("^MIR|SNOR", .$symbol),] %>% filter(!is.na(symbol)) %>% .$symbol %>% unique 

TotalGenesPV <- PVresults %>% .[!grepl("^MIR|SNOR", .$symbol),] %>% filter(!is.na(symbol)) %>% .$symbol %>% unique
TotalGenesPVup <- PVresults %>% filter(padj < 0.05, log2FoldChange > 0) %>%
  .[!grepl("^MIR|SNOR", .$symbol),] %>% filter(!is.na(symbol)) %>% .$symbol %>% unique
TotalGenesPVdown <- PVresults %>% filter(padj < 0.05, log2FoldChange < 0) %>%
  .[!grepl("^MIR|SNOR", .$symbol),] %>% filter(!is.na(symbol)) %>% .$symbol %>% unique
  

CommonUp <- intersect(TotalGenesNBBup, TotalGenesPVup)
CommonDown <- intersect(TotalGenesNBBdown, TotalGenesPVdown)

pCommonUp = dhyper(x = length(CommonUp),
                   m = length(TotalGenesPVup),
                   n = c(length(TotalGenesPV) - length(TotalGenesPVup)),
                   k = length(TotalGenesNBBup))

pCommonDown = dhyper(x = length(CommonDown),
                   m = length(TotalGenesPVdown),
                   n = c(length(TotalGenesPV) - length(TotalGenesPVdown)),
                   k = length(TotalGenesNBBdown))

ManhattanList <- list(PV = PVresults %>% filter(!duplicated(PeakName)),
                      NBB = NBBresults %>% filter(!duplicated(PeakName)))

ManhattanList <- lapply(ManhattanList, function(Cohort){
  Cohort$CHR <- sapply(Cohort$Peak.CHR, function(x){
  x = gsub("chr", "", x)
  x = gsub("x", 23, x, ignore.case = T)
  gsub("y", 24, x, ignore.case = T)
  }) %>% as.numeric
  pThresh = Cohort %>% filter(padj < 0.05) %>% arrange(desc(pvalue)) %>% .$pvalue %>% .[1] %>% -log10(.)
  list(Results = Cohort,
       pThresh = pThresh)
})


sapply(names(ManhattanList), function(CohortName){
  Cohort = ManhattanList[[CohortName]]
  pdf(paste0(ResultsPath, "ManhattanPlot", CohortName, ".pdf"), useDingbats = F, width = 10, height = 4)
  manhattan2(Cohort$Results, chr = "CHR", bp = "Peak.START", p = "pvalue", snp = "PeakName", col = c("gray48", "gray73"),
            chrlabs = c(1:22, "X", "Y"), genomewideline = FALSE, suggestiveline = Cohort$pThresh)
  dev.off()
})

