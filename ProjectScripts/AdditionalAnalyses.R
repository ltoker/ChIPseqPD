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

MergedChIPrna %>% filter(GeneSymbol %in% GeneMembersPV$`mitochondrial translation`) %>% arrange(desc(CorCont.PV))



PVenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)





PVenrich2 <- gsr(scores = MergedChIPrna[MergedChIPrna$CorCont.PV > 0.5,], scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrich2 <- gsr(scores = MergedChIPrna %>% filter(CorCont.NBB > 0.6), scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)

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
PVnarrowPeak <- read.table("data/Peaks/PV_ALL_RNA.narrowPeakClean", header = F, sep = "\t")
names(PVnarrowPeak)[1:5] <- c("CHR", "START", "END", "PeakName", "Score")
PVnarrowPeak <- PVnarrowPeak[!grepl("GL|hs",PVnarrowPeak$CHR),]
PVnarrowPeak %<>%  mutate(CHR = paste0("chr", CHR))
PVnarrowPeak$DApvalue <- resultsPV$pvalue[match(PVnarrowPeak$PeakName, resultsPV$PeakName)]

NBBnarrowPeak <- read.table("data/Peaks/NBB_ALL.narrowPeakClean", header = F, sep = "\t")
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
                       pPvalue = -log10(pvalue),
                       pPvalue2 = pPvalue)
resultsNBB$pPvalue2[resultsNBB$DirectionChange == "Hypoacetylated"] <- -1*resultsNBB$pPvalue2[resultsNBB$DirectionChange == "Hypoacetylated"]


resultsPV <- merge(resultsPV, PVnarrowPeakByGene, by = "symbol", all.x = T, sort = F)
names(resultsPV) <- sapply(names(resultsPV), function(x) gsub(".PV", "", x))

resultsPV %<>% mutate(HDAC_Binding = HDAC1 + HDAC2 + HDAC6 + HDAC8 + SIRT6,
                      EP300_DeltaBinding = EP300 - HDAC_Binding,
                      EP300_PropBinding = EP300/(HDAC_Binding + 0.1),
                      pPvalue = -log10(pvalue),
                      pPvalue2 = pPvalue)
                       
resultsPV$pPvalue2[resultsPV$DirectionChange == "Hypoacetylated"] <- -1*resultsPV$pPvalue2[resultsPV$DirectionChange == "Hypoacetylated"]

lm(pPvalue2 ~ EP300 + HDAC_Binding + PDgene + PeakNum + log10(EffectiveLength),
   data = resultsPV %>% arrange(pvalue) %>% filter(!duplicated(symbol))) %>% summary

lm(pPvalue2 ~ EP300 + HDAC_Binding + PDgene,
   data = resultsNBB) %>% summary




