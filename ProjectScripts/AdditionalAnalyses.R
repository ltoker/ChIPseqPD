install_github("https://github.com/PavlidisLab/ermineR")
library(ermineR)
HumanAnno <- fread('https://gemma.msl.ubc.ca/annots/Generic_human_noParents.an.txt.gz', header = T) %>% data.frame()

PVenrich <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)
NBBenrich <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)

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

ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersPV$`mitochondrial respiratory chain complex I assembly`), aes(Group, Cor)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "mitochondrial respiratory chain complex I assembly", x = "") +

  ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersPV$`fatty acid beta-oxidation`), aes(Group, Cor, color = Group)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "fatty acid beta-oxidation", x = "") +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4)+
  geom_hline(yintercept = 0, color = "red") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~Cohort, scales = "free")


ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersPV$`mitochondrial respiratory chain complex I assembly`), aes(Group, Cor, color = Group)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "mitochondrial respiratory chain complex I assembly", x = "") +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4)+
  geom_hline(yintercept = 0, color = "red") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~Cohort, scales = "free")

ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersNBB$`organellar large ribosomal subunit`), aes(Group, Cor, color = Group)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "organellar large ribosomal subunit", x = "") +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4)+
  geom_hline(yintercept = 0, color = "red") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~Cohort, scales = "free")


PVenrichCont <- ora(threshold = 0.8, scores = MergedChIPrna, scoreColumn = "CorCont.PV", aspects = c("M", "B"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrichCont <- ora(threshold = 0.8, scores = MergedChIPrna, scoreColumn = "CorCont.NBB", aspects = c("M", "B"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)



PVenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)





PVenrichDown <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)
NBBenrichDown <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)

SequencedDeep <- read.table("15794-1446155273_Covered.bed", header = F, sep= "\t", skip = 2)
names(SequencedDeep) <- c("CHR", "START", "END", "GeneSymbol")

SequencedDeep %<>% as("GRanges")

SignifMetaRegions <- CommonRegions %>% filter(AdjmetaP < 0.05) %>%
  select(UniqueRegionGeneType, AdjmetaP, padj.PV, padj.NBB) 

SignifMetaRegions$CHR <- sapply(SignifMetaRegions$UniqueRegionGeneType, function(x){
  strsplit(x, ":")[[1]][1]
})

SignifMetaRegions$START <- sapply(SignifMetaRegions$UniqueRegionGeneType, function(x){
  strsplit(x, ":|-")[[1]][2] 
})

SignifMetaRegions$END <- sapply(SignifMetaRegions$UniqueRegionGeneType, function(x){
  strsplit(x, ":|-|_")[[1]][3] 
})

SignifMetaRegions %<>% as("GRanges")

temp <- mergeByOverlaps(SequencedDeep, SignifMetaRegions) %>%
  data.frame() %>% select(-c(5,6, 8:16))

write.table(temp, "DeepSequencedSignifMetaP.tsv", row.names = F, col.names = T, sep = "\t")

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


ENCODEdata <- read.table("data/wgEncodeRegTfbsClusteredV3.bed", header = F, sep = "\t")
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
names(resultsNBB) <- sapply(names(resultsNBB), function(x) gsub(".NBB", "", x))
resultsNBB %<>% mutate(HDAC_Binding = HDAC1 + HDAC2 + HDAC6 + HDAC8 + SIRT6,
                       EP300_DeltaBinding = EP300 - HDAC_Binding,
                       EP300_PropBinding = EP300/(HDAC_Binding + 0.1),
                       pPvalue = -log10(pvalue))

resultsNBB$pPvalue2[resultsNBB$DirectionChange == "Hypoacetylated"] <- 1-resultsNBB$pvalue[resultsNBB$DirectionChange == "Hypoacetylated"]/2
resultsNBB$pPvalue2[resultsNBB$DirectionChange == "Hyperacetylated"] <- resultsNBB$pvalue[resultsNBB$DirectionChange == "Hyperacetylated"]/2
resultsNBB %<>% mutate(pPvalue2 = -log10(pPvalue2))

resultsNBB <- merge(resultsNBB, NBBnarrowPeak %>% filter(!duplicated(PeakName.NBB)) %>% select(-symbol), by.x = "PeakName", by.y = "PeakName.NBB", all.x = T, sort = F)
resultsNBB %<>% mutate(HDAC_Binding.NBB = HDAC1.NBB + HDAC2.NBB + HDAC6.NBB + HDAC8.NBB + SIRT6.NBB)


resultsPV <- merge(resultsPV, PVnarrowPeakByGene, by = "symbol", all.x = T, sort = F)
names(resultsPV) <- sapply(names(resultsPV), function(x) gsub(".PV", "", x))

resultsPV %<>% mutate(HDAC_Binding = HDAC1 + HDAC2 + HDAC6 + HDAC8 + SIRT6,
                      EP300_DeltaBinding = EP300 - HDAC_Binding,
                      EP300_PropBinding = EP300/(HDAC_Binding + 0.1),
                      pPvalue = -log10(pvalue))
                       
resultsPV$pPvalue2[resultsPV$DirectionChange == "Hypoacetylated"] <- 1-resultsPV$pvalue[resultsPV$DirectionChange == "Hypoacetylated"]/2
resultsPV$pPvalue2[resultsPV$DirectionChange == "Hyperacetylated"] <- resultsPV$pvalue[resultsPV$DirectionChange == "Hyperacetylated"]/2
resultsPV %<>% mutate(pPvalue2 = -log10(pPvalue2))

resultsPV <- merge(resultsPV, PVnarrowPeak %>% filter(!duplicated(PeakName.PV)) %>% select(-symbol), by.x = "PeakName", by.y = "PeakName.PV", all.x = T, sort = F)
resultsPV %<>% mutate(HDAC_Binding.PV = HDAC1.PV + HDAC2.PV + HDAC6.PV + HDAC8.PV + SIRT6.PV)

lm(pPvalue2 ~ EP300 + HDAC1 + HDAC2 + HDAC6 + HDAC8 + SIRT6  + PDgene,
   data = resultsPV) %>% summary

lm(pPvalue2 ~ EP300.PV + HDAC1.PV + HDAC2.PV + HDAC6.PV + HDAC8.PV + SIRT6.PV  + PDgene,
   data = resultsPV  %>% filter(!duplicated(PeakName))) %>% summary

lm(pPvalue2 ~ EP300.PV + HDAC_Binding.PV  + PDgene,
   data = resultsPV  %>% filter(!duplicated(PeakName))) %>% summary

lm(pPvalue2 ~ EP300 + HDAC1 + HDAC2 + HDAC6 + HDAC8 + SIRT6  + PDgene,
   data = resultsNBB) %>% summary

lm(pPvalue2 ~ EP300.NBB + HDAC1.NBB + HDAC2.NBB + HDAC6.NBB + HDAC8.NBB + SIRT6.NBB  + PDgene,
   data = resultsNBB  %>% filter(!duplicated(PeakName))) %>% summary

lm(pPvalue2 ~ EP300.NBB + HDAC_Binding.NBB  + PDgene,
   data = resultsNBB  %>% filter(!duplicated(PeakName))) %>% summary


lm(-log10(metaP) ~ EP300.NBB_Region + HDAC1.NBB_Region + HDAC2.NBB_Region + HDAC6.NBB_Region + HDAC8.NBB_Region + SIRT6.NBB_Region  + PDgene +
   EP300.PV_Region + HDAC1.PV_Region + HDAC2.PV_Region + HDAC6.PV_Region + HDAC8.PV_Region + SIRT6.PV_Region + log10(RegionLength),
   data = CommonRegions %>% filter(DirectionChange == "Hypoacetylated", RegionLength > 0)) %>% summary


