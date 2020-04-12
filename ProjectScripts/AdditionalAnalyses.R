source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
packageF("ggpubr")

ResultsPath = "GeneralResults/"

load("GeneralResults/CombinedData.Rda")
resultsNBB <- readRDS("ResultsNBB_Final//NBBDEresults.Rds") %>% arrange(padj)
resultsPV <- readRDS("ResultsPV_Final//PVDEresults.Rds") %>% arrange(padj) 
PDgenes <- read.table(paste0(ResultsPath, "PDgeneStat.tsv"), sep = "\t", header = TRUE)

install_github("https://github.com/PavlidisLab/ermineR")
library(ermineR)
HumanAnno <- fread('https://gemma.msl.ubc.ca/annots/Generic_human_noParents.an.txt.gz', header = T) %>% data.frame()

PVenrich <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)
NBBenrich <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)

write.table(PVenrich$results %>% data.frame()
            %>% select(-GeneMembers, -Same.as, -ID, -NumProbes, -RawScore, -Multifunctionality) %>%
              arrange(Pval) %>% filter(CorrectedPvalue < 0.05),
            "GeneralResults/EnrichDeltaPV.tsv", row.names = F, col.names = T, sep = "\t")

write.table(NBBenrich$results %>% data.frame()
            %>% select(-GeneMembers, -Same.as, -ID, -NumProbes, -RawScore, -Multifunctionality) %>%
              arrange(Pval) %>% filter(CorrectedPvalue < 0.05),
            "GeneralResults/EnrichDeltaNBB.tsv", row.names = F, col.names = T, sep = "\t")

GeneMembersPV <- PVenrich$results %>% data.frame() %>% filter(CorrectedPvalue < 0.05) %>% arrange(CorrectedMFPvalue) %>% .$GeneMembers %>% strsplit(.,"\\|")
names(GeneMembersPV) <- PVenrich$results %>% data.frame() %>% filter(CorrectedPvalue < 0.05) %>% arrange(CorrectedMFPvalue) %>% .$Name

GeneMembersNBB <- NBBenrich$results %>% data.frame() %>% filter(CorrectedPvalue < 0.05) %>% arrange(CorrectedMFPvalue) %>% .$GeneMembers %>% strsplit(.,"\\|")
names(GeneMembersNBB) <- NBBenrich$results %>% data.frame() %>% filter(CorrectedPvalue < 0.05) %>% arrange(CorrectedMFPvalue) %>% .$Name

MergedChIPrnaMelt <- gather(MergedChIPrna, key = Group, value = "Cor", matches("Cor"))

MergedChIPrnaMelt$Cohort <-  sapply(MergedChIPrnaMelt$Group, function(x){
  strsplit(x, "\\.")[[1]][2]
})
MergedChIPrnaMelt$Group <-  sapply(MergedChIPrnaMelt$Group, function(x){
  strsplit(x, "\\.")[[1]][1]
})

MergedChIPrnaMelt$Cohort <- factor(MergedChIPrnaMelt$Cohort, levels = c("PV", "NBB"))

p1 <- ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersPV$`fatty acid beta-oxidation`), aes(Group, Cor, color = Group)) +
  theme_minimal() +
  labs(title = "fatty acid beta-oxidation", x = "") +
  geom_violin(show.legend = F) +
  geom_boxplot(outlier.shape = NA, width = 0.3,aes(fill = Group), alpha = 0.4, show.legend = F) +
  geom_jitter(width = 0.2, alpha = 0.8, show.legend = F)+
  geom_hline(yintercept = 0, color = "red") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~Cohort)

p2 <- ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersPV$`mitochondrial respiratory chain complex I assembly`), aes(Group, Cor, color = Group)) +
  theme_minimal() +
  labs(title = "mitochondrial respiratory chain complex I assembly", x = "") +
  geom_violin(show.legend = F) +
  geom_boxplot(outlier.shape = NA, width = 0.3,aes(fill = Group), alpha = 0.4, show.legend = F) +
  geom_jitter(width = 0.2, alpha = 0.8, show.legend = F)+
  geom_hline(yintercept = 0, color = "red") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~Cohort)


p3 <- ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersNBB$`mitochondrial large ribosomal subunit`), aes(Group, Cor, color = Group)) +
  theme_minimal() +
  labs(title = "mitochondrial large ribosomal subunit", x = "") +
  geom_violin(show.legend = F) +
  geom_boxplot(outlier.shape = NA, width = 0.3,aes(fill = Group), alpha = 0.4, show.legend = F) +
  geom_jitter(width = 0.2, alpha = 0.8, show.legend = F)+
  geom_hline(yintercept = 0, color = "red") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~Cohort)

ggarrange(p2,p1,p3, nrow = 1)
ggsave(paste0("MitoGenesCor.pdf"), device = "pdf", width = 10, height = 2, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()


PVenrichPDdown <- gsr(scores = MergedChIPrna, scoreColumn = "CorPD.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)
NBBenrichPDdown <- gsr(scores = MergedChIPrna, scoreColumn = "CorPD.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)

write.table(PVenrichPDdown$results %>% data.frame()
            %>% select(-GeneMembers, -Same.as, -ID, -NumProbes, -RawScore, -Multifunctionality) %>%
              arrange(Pval) %>% filter(CorrectedPvalue < 0.05),
            "GeneralResults/PVenrichPDdown.tsv", row.names = F, col.names = T, sep = "\t")

write.table(NBBenrichPDdown$results %>% data.frame()
            %>% select(-GeneMembers, -Same.as, -ID, -NumProbes, -RawScore, -Multifunctionality) %>%
              arrange(Pval) %>% filter(CorrectedPvalue < 0.05),
            "GeneralResults/NBBenrichPDdown.tsv", row.names = F, col.names = T, sep = "\t")


#Look at p300 and HDAC binding sites
PVnarrowPeak <- read.table("data/Peaks/PV_ALL_RNA.narrowPeakClean.gz", header = F, sep = "\t")
names(PVnarrowPeak)[1:5] <- c("CHR", "START", "END", "PeakName", "Score")
PVnarrowPeak <- PVnarrowPeak[!grepl("GL|hs",PVnarrowPeak$CHR),]
PVnarrowPeak %<>%  mutate(CHR = paste0("chr", CHR))
PVnarrowPeak$DApvalue <- resultsPV$pvalue[match(PVnarrowPeak$PeakName, resultsPV$PeakName)]

NBBnarrowPeak <- read.table("data/Peaks/NBB_ALL.narrowPeakClean.gz", header = F, sep = "\t")
names(NBBnarrowPeak)[1:5] <- c("CHR", "START", "END", "PeakName", "Score")
NBBnarrowPeak <- NBBnarrowPeak[!grepl("GL|hs",NBBnarrowPeak$CHR),]
NBBnarrowPeak %<>%  mutate(CHR = paste0("chr", CHR))
NBBnarrowPeak$DApvalue <- resultsNBB$pvalue[match(NBBnarrowPeak$PeakName, resultsNBB$PeakName)]

UCSCrescale <- data.frame(OrgPvalue = c(PVnarrowPeak$DApvalue, NBBnarrowPeak$DApvalue)) %>%
  mutate(pPvalue = -log10(OrgPvalue))
UCSCrescale$RescaledVal <- rescale(UCSCrescale$pPvalue, c(100, 1000))
UCSCrescale %<>% filter(!duplicated(OrgPvalue))

PVnarrowPeak <- merge(PVnarrowPeak, UCSCrescale, by.x = "DApvalue", by.y = "OrgPvalue", sort = F)
PVnarrowPeak %<>% mutate(Score = RescaledVal)  %>% filter(!is.na(DApvalue))
PVnarrowPeakOrg <- PVnarrowPeak

NBBnarrowPeak <- merge(NBBnarrowPeak, UCSCrescale, by.x = "DApvalue", by.y = "OrgPvalue", sort = F)
NBBnarrowPeak %<>% mutate(Score = RescaledVal) %>% filter(!is.na(DApvalue))
NBBnarrowPeakOrg <- NBBnarrowPeak

write.table(NBBnarrowPeak %>% select(-DApvalue, -pPvalue,  -RescaledVal), "data/Peaks/NBB_ALL.narrowPeakClean2b", row.names = F, col.names = F, sep = "\t")
write.table(PVnarrowPeak %>% select(-DApvalue, -pPvalue,  -RescaledVal), "data/Peaks/PV_ALL_RNA.narrowPeakClean2b", row.names = F, col.names = F, sep = "\t")

download.file("https://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz",
              destfile = "data/wgEncodeRegTfbsClusteredV3.bed.gz")

ENCODEdata <- read.table("data/wgEncodeRegTfbsClusteredV3.bed.gz", header = F, sep = "\t")
names(ENCODEdata) <- c("CHR", "START", "END", "Name", "SCORE", "BlockCount", "BlockSizes", "BlockStarts")
ENCODEdata <- ENCODEdata[grepl("EP300|^HDAC|^SIRT", ENCODEdata$Name),]
ENCODEdata %<>% droplevels()


for(name in levels(ENCODEdata$Name)){
  temp <- ENCODEdata[ENCODEdata$Name == name,]
  temp2 <- findOverlaps(query = as(NBBnarrowPeak, "GRanges"), subject = as(temp, "GRanges")) %>%
    data.frame %>% group_by(queryHits) %>% summarise(n = n()) %>% data.frame()
  names(temp2) <- c("RowNum", name)
  NBBnarrowPeak[[name]] <- 0
  NBBnarrowPeak[[name]][temp2$RowNum] <- temp2[[name]]
}

for(name in levels(ENCODEdata$Name)){
  temp <- ENCODEdata[ENCODEdata$Name == name,]
  temp2 <- findOverlaps(query = as(PVnarrowPeak, "GRanges"), subject = as(temp, "GRanges")) %>%
    data.frame %>% group_by(queryHits) %>% summarise(n = n()) %>% data.frame()
  names(temp2) <- c("RowNum", name)
  PVnarrowPeak[[name]] <- 0
  PVnarrowPeak[[name]][temp2$RowNum] <- temp2[[name]]
}



NBBnarrowPeak %<>% select(matches("PeakN|EP|HDAC|SIRT"))
names(NBBnarrowPeak) <- paste0(names(NBBnarrowPeak), ".NBB")
PVnarrowPeak %<>% select(matches("PeakN|EP|HDAC|SIRT"))
names(PVnarrowPeak) <- paste0(names(PVnarrowPeak), ".PV")

CommonRegions <- merge(CommonRegions, NBBnarrowPeak, by = "PeakName.NBB")
CommonRegions <- merge(CommonRegions, PVnarrowPeak, by = "PeakName.PV")

NBBnarrowPeak <- merge(NBBnarrowPeak, resultsNBB %>% select(PeakName, symbol), by.x = "PeakName.NBB", by.y = "PeakName", sort = F)
PVnarrowPeak <-  merge(PVnarrowPeak, resultsPV %>% select(PeakName, symbol), by.x = "PeakName.PV", by.y = "PeakName", sort = F)

NBBnarrowPeakByGene <- NBBnarrowPeak %>% mutate(Filter = paste0(PeakName.NBB, symbol)) %>%
  filter(!duplicated(Filter)) %>%
  select(-PeakName.NBB, -Filter) %>% group_by(symbol) %>% summarise_all(sum)
                                                  
PVnarrowPeakByGene <- PVnarrowPeak %>% mutate(Filter = paste0(PeakName.PV, symbol)) %>%
  filter(!duplicated(Filter)) %>%
  select(-PeakName.PV, -Filter) %>% group_by(symbol) %>% summarise_all(sum)

  
CommonRegions <- merge(CommonRegions, NBBnarrowPeakByGene, by = "symbol", all.x = T, sort = F, suffixes = c("_Region", "_Gene"))
CommonRegions <- merge(CommonRegions, PVnarrowPeakByGene, by = "symbol", all.x = T, sort = F, suffixes = c("_Region", "_Gene"))

resultsNBB <- merge(resultsNBB, NBBnarrowPeakByGene, by = "symbol", all.x = T, sort = F)
names(resultsNBB) <- sapply(names(resultsNBB), function(x) gsub(".NBB", "_Gene", x))

resultsNBB %<>% mutate(HDAC_Binding_Gene = HDAC1_Gene + HDAC2_Gene + HDAC6_Gene + HDAC8_Gene + SIRT6_Gene,
                       EP300_DeltaBinding_Gene = EP300_Gene - HDAC_Binding_Gene,
                       pPvalue = -log10(pvalue))



resultsNBB$pPvalue2 <- NA
resultsNBB$pPvalue2[resultsNBB$log2FoldChange < 0] <- 1-resultsNBB$pvalue[resultsNBB$log2FoldChange < 0]/2
resultsNBB$pPvalue2[resultsNBB$log2FoldChange > 0] <- resultsNBB$pvalue[resultsNBB$log2FoldChange > 0]/2
resultsNBB %<>% mutate(pPvalue2 = -log10(pPvalue2))

resultsNBB <- merge(resultsNBB, NBBnarrowPeak %>% filter(!duplicated(PeakName.NBB)) %>% select(-symbol), by.x = "PeakName", by.y = "PeakName.NBB", all.x = T, sort = F)
resultsNBB %<>% mutate(HDAC_Binding.NBB = HDAC1.NBB + HDAC2.NBB + HDAC6.NBB + HDAC8.NBB + SIRT6.NBB)
names(resultsNBB) <- sapply(names(resultsNBB), function(x) gsub(".NBB", "_Peak", x))


resultsPV <- merge(resultsPV, PVnarrowPeakByGene, by = "symbol", all.x = T, sort = F)
names(resultsPV) <- sapply(names(resultsPV), function(x) gsub(".PV", "_Gene", x))

resultsPV %<>% mutate(HDAC_Binding_Gene = HDAC1_Gene + HDAC2_Gene + HDAC6_Gene + HDAC8_Gene + SIRT6_Gene,
                      EP300_DeltaBinding_Gene = EP300_Gene - HDAC_Binding_Gene,
                      pPvalue = -log10(pvalue))
                       
resultsPV$pPvalue2 <- NA
resultsPV$pPvalue2[resultsPV$log2FoldChange < 0] <- 1-resultsPV$pvalue[resultsPV$log2FoldChange < 0]/2
resultsPV$pPvalue2[resultsPV$log2FoldChange > 0] <- resultsPV$pvalue[resultsPV$log2FoldChange > 0]/2
resultsPV %<>% mutate(pPvalue2 = -log10(pPvalue2))

resultsPV <- merge(resultsPV, PVnarrowPeak %>% filter(!duplicated(PeakName.PV)) %>% select(-symbol), by.x = "PeakName", by.y = "PeakName.PV", all.x = T, sort = F)
resultsPV %<>% mutate(HDAC_Binding.PV = HDAC1.PV + HDAC2.PV + HDAC6.PV + HDAC8.PV + SIRT6.PV)
names(resultsPV) <- sapply(names(resultsPV), function(x) gsub(".PV", "_Peak", x))

#Add information about whether gene is a PDgene or not
resultsPV$PDgene <- "No"
resultsPV$PDgene[as.character(resultsPV$symbol) %in% as.character(PDgenes$Gene)] <- "Yes"

resultsNBB$PDgene <- "No"
resultsNBB$PDgene[as.character(resultsNBB$symbol) %in% as.character(PDgenes$Gene)] <- "Yes"

#Gene level
lm(pPvalue2 ~ EP300_Gene + HDAC1_Gene + HDAC2_Gene + HDAC6_Gene + HDAC8_Gene + SIRT6_Gene  + PDgene,
   data = resultsPV) %>% summary

lm(pPvalue2 ~ EP300_Gene + HDAC_Binding_Gene  + PDgene,
   data = resultsPV) %>% summary %>% .$coef %>%
  write.table(paste0(ResultsPath,"LmTFPV_Gene.tsv"), row.names = T, col.names = T, sep = "\t")

lm(pPvalue2 ~ EP300_Gene + HDAC1_Gene + HDAC2_Gene + HDAC6_Gene + HDAC8_Gene + SIRT6_Gene  + PDgene,
   data = resultsNBB) %>% summary

lm(pPvalue2 ~ EP300_Gene + HDAC_Binding_Gene  + PDgene,
   data = resultsNBB) %>% summary %>% .$coef %>%
  write.table(paste0(ResultsPath,"LmTFNBB_Gene.tsv"), row.names = T, col.names = T, sep = "\t")

#Peak level
lm(pPvalue2 ~ EP300_Peak + HDAC1_Peak + HDAC2_Peak + HDAC6_Peak + HDAC8_Peak + PDgene,
   data = resultsPV %>% filter(!duplicated(Peak_Gene))) %>% summary

lm(pPvalue2 ~ EP300_Peak + HDAC_Binding_Peak  + PDgene,
   data = resultsPV %>% filter(!duplicated(Peak_Gene))) %>% summary %>% .$coef %>%
  write.table(paste0(ResultsPath,"LmTFPV_Peak.tsv"), row.names = T, col.names = T, sep = "\t")

lm(pPvalue2 ~ EP300_Peak + HDAC1_Peak + HDAC2_Peak + HDAC6_Peak + HDAC8_Peak + PDgene,
   data = resultsNBB %>% filter(!duplicated(Peak_Gene))) %>% summary 

lm(pPvalue2 ~ EP300_Peak + HDAC_Binding_Peak  + PDgene,
   data = resultsNBB %>% filter(!duplicated(Peak_Gene))) %>% summary %>% .$coef %>%
  write.table(paste0(ResultsPath,"LmTFNBB_Peak.tsv"), row.names = T, col.names = T, sep = "\t")



PDgeneSignif <- PDgenes %>% filter(AdjMetaP < 0.05) %>% .$Gene %>% as.character()
PDgeneSignif <- PDgenes %>% filter(AdjMetaP < 0.05) %>% .$Gene %>% as.character()

resultsPV$MetaSignif <- "No"
resultsPV$MetaSignif[resultsPV$PeakName %in% (CommonRegions %>% filter(AdjmetaP < 0.05) %>% .$PeakName.PV)] <- "Yes"
resultsPV$MetaSignif[!resultsPV$PeakName %in% CommonRegions$PeakName.PV] <- NA


resultsNBB$MetaSignif <- "No"
resultsNBB$MetaSignif[resultsNBB$PeakName %in% (CommonRegions %>% filter(AdjmetaP < 0.05) %>% .$PeakName.NBB)] <- "Yes"
resultsNBB$MetaSignif[!resultsNBB$PeakName %in% CommonRegions$PeakName.NBB] <- NA


HDACdataPlotData <- rbind(resultsPV %>% filter(!duplicated(Peak_Gene)) %>%
                            select(pPvalue2, padj, PDgene, EP300_Peak, HDAC_Binding_Peak, MetaSignif) %>%
                            mutate(Cohort = "PW"),
                          resultsNBB %>% filter(!duplicated(Peak_Gene)) %>%
                            select(pPvalue2, padj, PDgene, EP300_Peak, HDAC_Binding_Peak, MetaSignif) %>%
                            mutate(Cohort = "NBB"))

HDACdataPlotData$CohortSignif <- sapply(HDACdataPlotData$padj, function(x){
  if(is.na(x)){
    NA
  } else if(x<0.05){
    "Yes"
  } else {
    "No"
  }
})


KSstatNBB <- ks.test(HDACdataPlotData %>% filter(CohortSignif == "No", Cohort == "NBB") %>% 
                       filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0) %>% mutate(Val = log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))) %>% .$Val,
                     HDACdataPlotData %>% filter(CohortSignif == "Yes", Cohort == "NBB") %>% 
                       filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0) %>% mutate(Val = log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))) %>% .$Val)

KSstatPV <- ks.test(HDACdataPlotData %>% filter(CohortSignif == "No", Cohort == "PW") %>% 
                      filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0) %>% mutate(Val = log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))) %>% .$Val,
                    HDACdataPlotData %>% filter(CohortSignif == "Yes", Cohort == "PW") %>% 
                      filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0) %>% mutate(Val = log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))) %>% .$Val)

KSdf <- data.frame(Cohort = c("NBB", "PW"),
                   Text = c(paste0("D = ", round(KSstatNBB$statistic, digits = 2), ", p = ", signif(KSstatNBB$p.value, digits = 2)),
                            paste0("D = ", round(KSstatPV$statistic, digits = 2), ", p = ", signif(KSstatPV$p.value, digits = 2))),
                   x = 0, y = 2.55)

ggplot(HDACdataPlotData %>% filter(!is.na(CohortSignif)) %>% 
         filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0),
       aes(log((EP300_Peak + 0.01)/(HDAC_Binding_Peak+0.01)), fill = CohortSignif)) +
  theme_minimal(base_size = 14) +
  labs(x = "log(p300_binding/HDAC_binding)", title = "All Genes") +
  scale_fill_manual(values = c("grey", "darkred")) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Cohort) +
  geom_text(data = KSdf, aes(x, y, label = Text), inherit.aes = F, size = 4)
ggsave(paste0("KStest.pdf"), device = "pdf", width = 10, height = 5, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()


WilcoxNBBCohort <- wilcox.test(log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))~CohortSignif,
                             data = HDACdataPlotData %>% filter(!is.na(CohortSignif), Cohort == "NBB", PDgene == "Yes") %>%
                               filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0))
WilcoxPVCohort <- wilcox.test(log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))~CohortSignif,
                               data = HDACdataPlotData %>% filter(!is.na(CohortSignif), Cohort == "PW", PDgene == "Yes") %>%
                                 filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0))
WilcoxCohortdf <- data.frame(Cohort = c("NBB", "PW"),
                             Text = c(paste0("p = ", signif(WilcoxNBBCohort$p.value, digits = 2)),
                                      paste0("p = ", signif(WilcoxPVCohort$p.value, digits = 2))),
                             x = c(1,2), y = 2.1)


WilcoxNBBMeta <- wilcox.test(log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))~MetaSignif,
                               data = HDACdataPlotData %>% filter(!is.na(MetaSignif), Cohort == "NBB", PDgene == "Yes") %>%
                                 filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0))
WilcoxPVMeta <- wilcox.test(log((EP300_Peak+0.01)/log(HDAC_Binding_Peak + 0.01))~MetaSignif,
                              data = HDACdataPlotData %>% filter(!is.na(MetaSignif), Cohort == "PW", PDgene == "Yes") %>%
                                filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0))

WilcoxMetadf <- data.frame(Cohort = c("NBB", "PW"),
                             Text = c(paste0("p = ", signif(WilcoxNBBMeta$p.value, digits = 2)),
                                      paste0("p = ", signif(WilcoxPVMeta$p.value, digits = 2))),
                             x = c(1,2), y = 2.1)



ggplot(HDACdataPlotData %>% filter(!is.na(MetaSignif), PDgene == "Yes") %>% 
         filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0),
       aes(Cohort, log(EP300_Peak + 0.01)-log(HDAC_Binding_Peak+0.01),fill = MetaSignif, color = MetaSignif)) +
  theme_minimal() +
  labs(y = "log(p300_binding/HDAC_binding)", x = "") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05)) +
  scale_fill_manual(values = c("grey", "darkred")) +
  scale_color_manual(values = c("grey40", "darkred")) +
  geom_text(data = WilcoxMetadf, aes(x, y, label = Text), inherit.aes = F, size = 4)
ggsave(paste0("MetaSignif_TF.pdf"), device = "pdf", width = 8, height = 5, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()

ggplot(HDACdataPlotData %>% filter(!is.na(CohortSignif), PDgene == "Yes") %>% 
         filter(EP300_Peak > 0 & HDAC_Binding_Peak > 0),
       aes(Cohort, log(EP300_Peak + 0.01)-log(HDAC_Binding_Peak+0.01),fill = CohortSignif, color = CohortSignif)) +
  theme_minimal() +
  labs(y = "log(p300_binding/HDAC_binding)", x = "") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05)) +
  scale_fill_manual(values = c("grey", "darkred")) +
  scale_color_manual(values = c("grey40", "darkred")) +
  geom_text(data = WilcoxCohortdf, aes(x, y, label = Text), inherit.aes = F, size = 4)
ggsave(paste0("CohortSignif_TF.pdf"), device = "pdf", width = 8, height = 5, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()

saveRDS(sessionInfo(), file = "SessionInfo.RDS")



