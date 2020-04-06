source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
plotMA = DESeq2::plotMA
packageF("org.Hs.eg.db")


if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

load("meta/Metadata.Rda")

Metadata %<>% SampleFilter

Metadata %<>% mutate(SampleID = paste0("X", activemotif_id))
Metadata$RNAseq <- sapply(Metadata$sample_id_rna, function(x){
  if(is.na(x)){
    "No"
    } else {
      "Yes"
      }
  }) %>% factor(levels = c("Yes", "No"))

#Filter only the relevant samples
MetaChipSeq %<>% filter(SampleID %in% Metadata$activemotif_id) %>% droplevels()

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")


################## Update ChIP Metadata with file locations  ##############################################

MetaChipSeqGroupCalled <- MetaChipSeq
MetaChipSeqGroupCalled$Peaks[grepl("Case", MetaChipSeqGroupCalled$Condition)] <- CasePeakLoc
MetaChipSeqGroupCalled$Peaks[grepl("Control", MetaChipSeqGroupCalled$Condition)] <- ControlPeakLoc

MetaChipSeqAllcalled <- MetaChipSeq

MetaChipSeqAllcalled$Peaks <- CombinedPeakLoc

###################### Analysis of peaks called based on grouped samples##############
SampleGroup <- Metadata %>% group_by(condition) %>% summarise(n = n()) %>% data.frame()

InputPeakGroupCalled <- list()

InputPeakGroupCalled$PeakData <- as.list(unique(MetaChipSeqGroupCalled$Peaks))
names(InputPeakGroupCalled$PeakData) <- unique(MetaChipSeqGroupCalled$Condition) %>% as.character
names(InputPeakGroupCalled$PeakData) <- sapply(names(InputPeakGroupCalled$PeakData), function(x){
  if(x == "Case"){
    "PD"
  } else if(x == "Control"){
    "Cont"
  }
})

InputPeakGroupCalled$PeakData <- mclapply(InputPeakGroupCalled$PeakData, function(grp){
  temp <- read.table(grp, header = F, sep = "\t")
  names(temp) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")[1:ncol(temp)]
  temp$CHR <- sapply(as.character(temp$CHR), function(x) paste0("chr", x))
  temp <- temp[!grepl("GL|hs", temp$CHR),] %>% droplevels()
  temp %<>% mutate(pValue = 10^(-pPvalue)) %>% filter(pValue < 10^(-7))
  temp %>% arrange(CHR, START) %>% as(.,"GRanges")
}, mc.cores = detectCores()) 

GroupPeakOverlap <- findOverlaps(InputPeakGroupCalled$PeakData$PD, InputPeakGroupCalled$PeakData$Cont)
InputPeakGroupCalled$UniquePeaks <- list(PD = InputPeakGroupCalled$PeakData$PD[-queryHits(GroupPeakOverlap)],
                                         Cont = InputPeakGroupCalled$PeakData$Cont[-subjectHits(GroupPeakOverlap)])
InputPeakGroupCalled$UniquePeaks <- lapply(InputPeakGroupCalled$UniquePeaks, function(grp){
  PeakAnnoFile <- mergeByOverlaps(annoFileCollapsed, grp, maxgap = 0, type = "any", select = "all") %>% data.frame()
  names(PeakAnnoFile)[1:4] <- c("Region.CHR", "Region.START", "Region.END", "Region.width")
  PeakAnnoFile %<>% select(-matches("strand|_id|^id|^PeakName|annoFileCollapsed|offset", ignore.case = T))
  names(PeakAnnoFile) <- sapply(names(PeakAnnoFile), function(x) gsub("grp", "Peak", x))
  names(PeakAnnoFile)[grepl("PeakName", names(PeakAnnoFile))] <- "PeakName"
  PeakAnnoFile %>% mutate(Peak.Location = paste0(Peak.seqnames, ":", Peak.start, "-", Peak.end)) %>% data.frame
})


InputPeakGroupCalled$CommonPeaks <- list(PD = InputPeakGroupCalled$PeakData$PD[queryHits(GroupPeakOverlap)] %>% data.frame() %>% filter(!duplicated(.$PeakName)) %>% as(., "GRanges"),
                                         Cont = InputPeakGroupCalled$PeakData$Cont[subjectHits(GroupPeakOverlap)] %>% data.frame() %>% filter(!duplicated(.$PeakName))  %>% as(., "GRanges"))
InputPeakGroupCalled$CommonPeaks <- lapply(InputPeakGroupCalled$CommonPeaks, function(grp){
  PeakAnnoFile <- mergeByOverlaps(annoFileCollapsed, grp, maxgap = 0, type = "any", select = "all") %>% data.frame()
  names(PeakAnnoFile)[1:4] <- c("Region.CHR", "Region.START", "Region.END", "Region.width")
  PeakAnnoFile %<>% select(-matches("strand|_id|^id|^PeakName|annoFileCollapsed|offset", ignore.case = T))
  names(PeakAnnoFile) <- sapply(names(PeakAnnoFile), function(x) gsub("grp", "Peak", x))
  names(PeakAnnoFile)[grepl("PeakName", names(PeakAnnoFile))] <- "PeakName"
  PeakAnnoFile %>% mutate(Peak.Location = paste0(Peak.seqnames, ":", Peak.start, "-", Peak.end)) %>% data.frame
})

pPvalueDF <- sapply(c("CommonPeaks", "UniquePeaks"), function(peakType){
  sapply(c("Cont", "PD"), function(Group){
    data.frame(pPvalue = InputPeakGroupCalled[[peakType]][[Group]] %>%
                 filter(!duplicated(PeakName), Peak.pValue > 0) %>% .$pPvalue,
               PeakTypeName = paste0(gsub("Peaks", "", peakType), Group),
               PeakType = peakType,
               Group = Group)
  }, simplify = F) %>% do.call(rbind,.)
}, simplify = F) %>% do.call(rbind,.)   

Plot <- ggplot(pPvalueDF, aes(PeakTypeName, pPvalue)) +
  theme_classic() +
  labs(x = "", y = "-log10(p-value)") +
  geom_violin(width = 0.8, draw_quantiles = c(0.25, 0.5, 0.75), aes(fill = Group, color = PeakType), alpha = 0.6, scale = "width") +
  #geom_boxplot(outlier.shape = NA, aes(fill = Group, color = PeakType), alpha = 0.6) +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("black", "firebrick1"))
ggsave(paste0(ResultsPath,"PeakPvalueDist", Cohort, ".pdf"), plot = Plot,device = "pdf", width = 10, height = 6, dpi = 300)  
closeDev()  

GroupGenomeAnnoDF <- data.frame(GenomeAnno = unique(annoFileCollapsed$type))
CommonCont = table(InputPeakGroupCalled$CommonPeaks$Cont$type)[as.character(GroupGenomeAnnoDF$GenomeAnno)]
UniqueCont = table(InputPeakGroupCalled$UniquePeaks$Cont$type)[as.character(GroupGenomeAnnoDF$GenomeAnno)]                                                   
CommonPD = table(InputPeakGroupCalled$CommonPeaks$PD$type)[as.character(GroupGenomeAnnoDF$GenomeAnno)]
UniquePD = table(InputPeakGroupCalled$UniquePeaks$PD$type)[as.character(GroupGenomeAnnoDF$GenomeAnno)]  
GroupGenomeAnnoDF %<>% mutate(CommonCont = signif(100*CommonCont/nrow(InputPeakGroupCalled$CommonPeaks$Cont), digits = 2) %>% as.vector,
                              UniqueCont = signif(100*UniqueCont/nrow(InputPeakGroupCalled$UniquePeaks$Cont), digits = 2) %>% as.vector,
                              CommonPD = signif(100*CommonPD/nrow(InputPeakGroupCalled$CommonPeaks$PD), digits = 2) %>% as.vector,
                              UniquePD = signif(100*UniquePD/nrow(InputPeakGroupCalled$UniquePeaks$PD), digits = 2) %>% as.vector) %>% arrange(CommonCont)

GroupGenomeAnnoDF$GenomeAnno <- sapply(as.character(GroupGenomeAnnoDF$GenomeAnno), function(x) gsub("hg19_genes_", "", x))
GroupGenomeAnnoMelt <- gather(GroupGenomeAnnoDF, key = "PeakType", value = "Percent", -GenomeAnno)
GroupGenomeAnnoMelt$GenomeAnno <- factor(GroupGenomeAnnoMelt$GenomeAnno, levels = rev(c("intergenic", "introns", "1to5kb",  "3UTRs", "5UTRs", "intronexonboundaries", "exons", "promoters")))                                                           
Plot <- ggplot(GroupGenomeAnnoMelt, aes(PeakType, Percent, fill = GenomeAnno)) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "Annotations (%)") +
  scale_fill_manual(values = c("firebrick1", "darkorange1", "gold", "chartreuse3" ,"chartreuse4", "cyan", "blue4", "magenta3"), name = "Genomic annottation") +
  geom_bar(stat = "identity") 
ggsave(paste0(ResultsPath,"GenomAnnoUniquePeaks", Cohort, ".pdf"), plot = Plot,device = "pdf", width = 10, height = 6, dpi = 300)
closeDev()

Conditions = names(InputPeakGroupCalled$PeakData)

InputPeakGroupCalled$GroupStat <- data.frame(Condition = Conditions,
                                             TotalPeaksInGroup = lapply(InputPeakGroupCalled$PeakData, function(grp){
                                               grp %>% data.frame %>% nrow
                                             }) %>% unlist,
                                             TotalCoverage = lapply(InputPeakGroupCalled$PeakData, function(grp){
                                               grp %>% data.frame %>% .$width %>% sum
                                             }) %>% unlist,
                                             CommonPeaks = lapply(InputPeakGroupCalled$CommonPeaks, function(grp){
                                               grp %>% data.frame %>% filter(!duplicated(.$PeakName)) %>% nrow
                                             }) %>% unlist,
                                             UniquePeaks = lapply(InputPeakGroupCalled$UniquePeaks, function(grp){
                                               grp %>% data.frame %>% filter(!duplicated(.$PeakName)) %>% nrow
                                             }) %>% unlist)

                                               
InputPeakGroupCalled$GroupStat %<>% mutate(UniquePeakGroupPercent = signif(100*UniquePeaks/TotalPeaksInGroup, digits = 3),
                                           TotalCovPercent = signif(100*TotalCoverage/(2.7*10^9), digits = 3))
write.table(InputPeakGroupCalled$GroupStat, paste0(ResultsPath,"PeakStatGroupComparison-PeakData", Cohort, ".tsv"), row.names = F, col.names = T, sep = "\t")
DataToPlot <- InputPeakGroupCalled$GroupStat %>% select(Condition, TotalPeaksInGroup,
                                                        UniquePeakGroupPercent, TotalCovPercent) %>%
  gather(key = "Measure", value = "Value", -Condition)
DataToPlot$Measure <- factor(DataToPlot$Measure, levels = c("TotalPeaksInGroup", "UniquePeakGroupPercent", "TotalCovPercent"))
levels(DataToPlot$Measure) <- c("TotalPeaks", "UniquePeaks (%)", "GenomeCoverage (%)")

Plot <- ggplot(DataToPlot, aes(Condition, Value, fill = Condition)) +
  theme_classic(base_size = 16) +
  labs(x = "", y = "") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~Measure, nrow = 1, scales = "free")
ggsave(paste0("GroupPeaksDifference", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 4, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()
################ Repeat for peaks called on all samples combined ############

InputPeakAllCalled <- list()
InputPeakAllCalled$PeakData <- read.table(CombinedPeakLoc, header = F, sep = "\t")
names(InputPeakAllCalled$PeakData) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")[1:ncol(InputPeakAllCalled$PeakData)]
InputPeakAllCalled$PeakData %<>% mutate(pValue = 10^(-pPvalue)) %>% filter(pValue < 10^(-7))
InputPeakAllCalled$PeakData <- InputPeakAllCalled$PeakData[!grepl("GL|hs", InputPeakAllCalled$PeakData$CHR),] %>% droplevels()
InputPeakAllCalled$PeakData %<>% arrange(CHR, START) %>% as(.,"GRanges")



InputPeakAllCalled$Summary <- data.frame(TotalPeaks = InputPeakAllCalled$PeakData %>% data.frame %>% nrow,
                                         TotalCoverage = InputPeakAllCalled$PeakData %>% data.frame %>% .$width %>% sum)

InputPeakAllCalled$Summary %<>% mutate(TotalCovPercent = signif(100*TotalCoverage/(2.7*10^9), digits = 3))
InputPeakAllCalled$PeakData %>% data.frame() %>% .$width %>% summary()                                         

############################# HTseq counts ######################################################################
HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")

names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  x = gsub(".data.parkome.chipseq.results.bamfiles.0?", "", x)
  strsplit(x,"_")[[1]][1]
})
HTseqCounts %<>% mutate(Peak.Location = paste0("chr", Chr, ":", Start, "-", End))
names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", paste0("X", as.character(Metadata$activemotif_id))))

#Remove peaks in contig regions and peaks with p > 10^-7 (as well as blacklisted peaks)
HTseqCounts <- HTseqCounts %>% filter(Geneid %in% as.character(InputPeakAllCalled$PeakData$PeakName)) %>% droplevels()

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("^X")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  

AllCalledData <- GetCountMatrixHTseq(HTseqCounts, OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29")
closeDev()

##### Get relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
AllCalledData$SampleInfo <- GetCellularProportions(AllCalledData$SampleInfo)

pdf(paste0(ResultsPath, "SampleCorPromoter", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
countMatrixPromotersAllCalled  <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot, collapseBy = "GeneAnnoType",CorMethod = "pearson",
                                                     FilterBy = "promoter", meta = AllCalledData$SampleInfo,
                                                     title = paste0("Sample correlation (promoter), Peaks - called together, ", Cohort))
closeDev()

#Detect outliers
MedianCor <- apply(countMatrixPromotersAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]

Model = as.formula(" ~ condition + sex + age + batch + pm_hours + Oligo_MSP")

DESeqOutAll_promoters <- RunDESeq(data = countMatrixPromotersAllCalled$countMatrix,UseModelMatrix = T, sampleToFilter = paste(sapply(names(Outlier), function(x) gsub("X", "", x)), collapse = "|"),
                                  meta = countMatrixPromotersAllCalled$Metadata, normFactor = "MeanRatioOrg",
                                  FullModel = Model, test = "Wald", FitType = "local")

                                    
DESegResultsSex_promotersAll <- GetDESeqResults(DESeqOutAll_promoters, coef = "sexM") %>% mutate(symbol = sapply(.$PeakName, function(x) {strsplit(x, "_")[[1]][1]}))
DESegResultsAge_promotersAll <- GetDESeqResults(DESeqOutAll_promoters, coef = "age") %>% mutate(symbol = sapply(.$PeakName, function(x) {strsplit(x, "_")[[1]][1]}))
DESegResultsGroup_promotersAll <- GetDESeqResults(DESeqOutAll_promoters, coef = "conditionPD") %>% mutate(symbol = sapply(.$PeakName, function(x) {strsplit(x, "_")[[1]][1]}))

###########################################################################################################
############ RERUN USING PEAKS BASED ON ALL THE SAMPLES, without collapsing peaks #########################
###########################################################################################################
pdf(paste0(ResultsPath, "SampleCorAllPeaks", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",CorMethod = "pearson",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))
closeDev()

#Get the pvalues for associasion of each covariate with the first 3 PCs
PCAsamples <- prcomp(t(countMatrixFullAllCalled$CPMdata), scale. = T)
countMatrixFullAllCalled$Metadata %<>% mutate(PC1 = PCAsamples$x[,1],
                                              PC2 = PCAsamples$x[,2],
                                              PC3 = PCAsamples$x[,3],
                                              PC4 = PCAsamples$x[,4],
                                              PC5 = PCAsamples$x[,5]) 
VarExplained <- PCAsamples %>% summary() %>% .$importance %>%
  .[2, 1:sum(grepl("^PC", names(countMatrixFullAllCalled$Metadata)))]*100 


CovarPvalues <- sapply(grep("^PC", names(countMatrixFullAllCalled$Metadata), value = T), function(PC){
  temp <- lm(as.formula(paste0(PC, "~ condition + sex + age + batch + pm_hours + rin + Oligo_MSP + MeanRatioOrg + RiP_NormMeanRatioOrg")),
             data = countMatrixFullAllCalled$Metadata) %>% summary
  temp$coefficients[-1,4]
}, simplify = F) %>% do.call(cbind, .) %>% data.frame()

names(CovarPvalues) <- paste0(names(CovarPvalues), "(", round(VarExplained, digits = 1), "%)")
CovarPvalues %<>% mutate(Variable = factor(rownames(CovarPvalues), levels = rownames(CovarPvalues)))
levels(CovarPvalues$Variable) <- c("Group", "Sex", "Age", grep("batch", levels(CovarPvalues$Variable), value = T), "PMI", "RIN", "Oligo_MSP", "NormalizingFactor", "NormalizedRiP")
CovarPvaluesMelt <- gather(CovarPvalues, key = "PC", value = "pValue", -Variable)
CovarPvaluesMelt %<>% mutate(pPvalue = -log10(pValue)) 


Plot  <- ggplot(CovarPvaluesMelt, aes(PC, Variable)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0)) +
  labs(x = "", y = "", title = "Association of variables with main PCs") +
  geom_tile(aes(fill = pPvalue), colour = "white") +
  scale_fill_gradient(high = "steelblue", low = "gray94", name = "-log10(p)") +
  geom_text(aes(label = signif(pValue, 2))) 
ggsave(paste0("AssociationWithPCs", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 8, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()


countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|^X"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("^X")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)

################# Check the best normalization method ###############################

## Calculate normalization factors using RLE
DEStemp <- DESeqDataSetFromMatrix(countData = countMatrixDF %>% data.frame %>% select(matches("^X")) %>% as.matrix, colData = countMatrixFullAllCalled$Metadata, design = Model)
DEStemp <-  estimateSizeFactors(DEStemp)
DEStemp$RiP_RLE <- DEStemp$TotalCount/DEStemp$sizeFactor

MeltedDataAll <- DEStemp@colData %>% data.frame %>% gather(key = "MeasureType", value = "Value", TotalCount, library_size, RiP_NormAllCount, RiP_NormBackground, RiP_RLE, RiP_NormMeanRatioOrg, RiP_NormMeanRatioAll)
MeltedDataAll$MeasureType <- factor(MeltedDataAll$MeasureType, levels = c("TotalCount", "library_size",  "RiP_NormAllCount", "RiP_NormBackground", "RiP_RLE", "RiP_NormMeanRatioOrg", "RiP_NormMeanRatioAll"))
levels(MeltedDataAll$MeasureType) <- c("LibrarySize",  "RiP", "RiP/LibrarySize", "RiP/RoP", "RiP/sizeFactor", "RiP/MeanRatio", "RiP/MeanRatio2")

#Add individual vallues
xLabFun <- function(x) signif(as.numeric(as.character(x)), digits = 2)

StatNormMethod <- sapply(levels(MeltedDataAll$MeasureType), function(Type){
  data = MeltedDataAll %>% filter(!is.na(H3K27gapdh_Norm), !activemotif_id %in% c("57","39"), MeasureType == Type) %>% droplevels()
  cor.test(~H3K27gapdh_Norm+Value, data = data)
}, simplify = F)

MeltedDataAll %<>% mutate(MeasureType2 = MeasureType)

levels(MeltedDataAll$MeasureType2) <- sapply(levels(MeltedDataAll$MeasureType2), function(Type){
  temp <- StatNormMethod[[Type]]
  paste0(Type, " (r=", signif(temp$estimate, digits = 2), ", p=",  signif(temp$p.value, digits = 2), ")")
})

  
Plot <- ggplot(MeltedDataAll %>% filter(!is.na(H3K27gapdh_Norm), !activemotif_id %in% c("57", "39")), aes(H3K27gapdh_Norm, Value, color = condition)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Counts", x = "WB, H3K27gapdh_normalized (Final)", title = Cohort) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") + 
  #geom_smooth(method = "lm", color = "black") +
  scale_y_continuous(labels = xLabFun) +
  facet_wrap(~MeasureType2, scales = "free_y")
ggsave(paste0("WB_ChipSeqCor_FinalnormalizedToReplicates", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()


Plot <- ggplot(MeltedDataAll %>% filter(MeasureType == "RiP/MeanRatio", !activemotif_id %in% c("57", "39")), aes(condition, Value, color = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Normalized RiP", title = Cohort) +
  geom_boxplot(outlier.shape = NA, aes(fill = condition), alpha = 0.4) +
  geom_jitter(width = 0.2, aes(color = condition)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_y_continuous(labels = xLabFun)
ggsave(paste0("GroupNormRiP", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()

Plot <- ggplot(MeltedDataAll %>% filter(MeasureType == "RiP/MeanRatio", !activemotif_id %in% c("57", "39")), aes(condition, H3K27gapdh_Norm, color = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "WB H3K27ac/GAPDH ", title = Cohort) +
  geom_boxplot(outlier.shape = NA, aes(fill = condition), alpha = 0.4) +
  geom_jitter(width = 0.2, aes(color = condition)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_y_continuous(labels = xLabFun)
ggsave(paste0("GroupWB", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()


Plot <- ggplot(MeltedDataAll %>% filter(MeasureType == "RiP/MeanRatio", !activemotif_id %in% c("57", "39")), aes(age, Value, fill = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "Age", y = "Normalized RiP", title = Cohort) +
  geom_smooth(method = "lm", aes(fill = condition, color = condition), alpha = 0.3, size = 0.2) +
  geom_point(size = 2, aes(color = condition)) +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_y_continuous(labels = xLabFun) +
  facet_wrap(~condition, scales = "free_x")
ggsave(paste0("AgeRipCorrelation", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()


#Look at the correlation betweem expression-based and ChIP-seq based cell abundance estimates
CellData <- countMatrixFullAllCalled$Metadata %>% select(-matches("__|RiP|All|gapdh|Backgr|library|Total|h3")) %>% gather(matches("Genes"), key = "CellTypeMGP", value = "MGP")
CellData$CellTypeMGP <- sapply(CellData$CellTypeMGP, function(x) gsub("_Genes", "", x))
CellData %<>% gather(matches("MSP"), key = "CellTypeMSP", value = "MSP")
CellData$CellTypeMSP <- sapply(CellData$CellTypeMSP, function(x) gsub("_MSP", "", x))
CellData %<>% filter(CellTypeMGP == CellTypeMSP)
CellData %<>% mutate(MSPnorm = MSP/MeanRatioOrg)

#Add information regarding whether or not samples were analyzed in western blot
WesterAnalized <- c(21, 23, 25, 26, 28, 29, 30, 33, 34, 35, 36, 37, 42, 62)
CellData$Western <- sapply(CellData$activemotif_id, function(Subj){
  if(Subj %in% WesterAnalized){
    "YES"
    } else {
      "NO"
      }
  })

Plot <- ggplot(CellData %>% filter(!activemotif_id %in% c("57", "39")), aes(condition, MSPnorm, color = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  geom_boxplot(outlier.shape = NA, aes(fill = condition), alpha = 0.4) +
  geom_jitter(width = 0.2, aes(color = condition, pch = Western)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_shape_manual(values = c(16, 4)) +
  facet_wrap(~CellTypeMSP, scales = "free")

ggsave(paste0("CellTypeDifferences", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()

#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]

#Run the analysis
DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, 
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg", sampleToFilter = paste(sapply(names(Outlier), function(x) gsub("X", "", x)), collapse = "|"),
                             FullModel = Model, test = "Wald", FitType = "local")
                                  

DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "sexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "conditionPD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")


## Rerun with RLE normalization
DESeqOutAll_Full_RLE <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, 
                             meta = countMatrixFullAllCalled$Metadata, normFactor = NULL, sampleToFilter = paste(sapply(names(Outlier), function(x) gsub("X", "", x)), collapse = "|"),
                             FullModel = Model, test = "Wald", FitType = "local")


DESegResultsSex_FullAll_RLE <- GetDESeqResults(DESeqOutAll_Full_RLE, coef = "sexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll_RLE <- GetDESeqResults(DESeqOutAll_Full_RLE, coef = "age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll_RLE <- GetDESeqResults(DESeqOutAll_Full_RLE, coef = "conditionPD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")



#### Rerun without cell type correction
Model2 = as.formula(" ~ condition + sex + age + batch + pm_hours")

DESeqOutAll_Full_noCorrection <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, 
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg", sampleToFilter = paste(sapply(names(Outlier), function(x) gsub("X", "", x)), collapse = "|"),
                             FullModel = Model2, test = "Wald", FitType = "local")


DESegResultsSex_FullAll_noCorrection <- GetDESeqResults(DESeqOutAll_Full_noCorrection, coef = "sexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll_noCorrection <- GetDESeqResults(DESeqOutAll_Full_noCorrection, coef = "age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll_noCorrection <- GetDESeqResults(DESeqOutAll_Full_noCorrection, coef = "conditionPD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

#RLE normalization
DESeqOutAll_Full_RLE_noCorrection <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, 
                                 meta = countMatrixFullAllCalled$Metadata, normFactor = NULL, sampleToFilter = paste(sapply(names(Outlier), function(x) gsub("X", "", x)), collapse = "|"),
                                 FullModel = Model2, test = "Wald", FitType = "local")


DESegResultsSex_FullAll_RLE_noCorrection <- GetDESeqResults(DESeqOutAll_Full_RLE_noCorrection, coef = "sexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll_RLE_noCorrection <- GetDESeqResults(DESeqOutAll_Full_RLE_noCorrection, coef = "age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll_RLE_noCorrection <- GetDESeqResults(DESeqOutAll_Full_RLE_noCorrection, coef = "conditionPD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")




#Identify top genes
ModifiedDF <- DESegResultsGroup_FullAll
ModifiedDF$type <- sapply(ModifiedDF$type, function(x){
  rev(strsplit(x, "_")[[1]])[1]
}) %>% factor(levels = c("promoters", "exons", "5UTRs", "3UTR", "1to5kb","intronexonboundaries", "introns", "intergenic"))
ModifiedDF %<>% mutate(AnnoOrder = type)
levels(ModifiedDF$AnnoOrder) <- c(1:length(levels(ModifiedDF$AnnoOrder)))
ModifiedDF %<>% arrange(AnnoOrder) 
ModifiedDF %<>% filter(!duplicated(Peak_Gene))

SignifGenes <- ModifiedDF %>% filter(padj < 0.05) %>% .$symbol %>% unique
ModifiedDF %<>% filter(symbol %in% SignifGenes, !is.na(symbol), !is.na(padj))


TopGenes <- sapply(unique(ModifiedDF$symbol), function(gene){
  GenePeaks <- ModifiedDF %>% filter(symbol == gene)
  GenePeaks %<>% mutate(Adj2 = p.adjust(.$padj, method = "BH"))
  TotalPeaks <- nrow(GenePeaks)
  UpPeaks <- sum(GenePeaks$log2FoldChange > 0)
  DownPeaks <- sum(GenePeaks$log2FoldChange < 0)
  TopList <- if(UpPeaks/TotalPeaks == 1 | DownPeaks/TotalPeaks == 1 | "promoters" %in% c(GenePeaks %>%
                                                                                          filter(padj < 0.05) %>%
                                                                                          .$type %>% as.character())){
    "Yes"
  } else {
    "No"
  }
  NomSignif = sum(GenePeaks$pvalue < 0.05)
  AdjSignif = sum(GenePeaks$padj < 0.05)
  PromoterSignif <- if("promoters" %in% c(GenePeaks %>% filter(padj < 0.05) %>% .$type %>% as.character())){
    "Yes"
  } else if("promoters" %in% c(GenePeaks %>% filter(padj > 0.05) %>% .$type %>% as.character())){
    "No"
  } else {
    NA
  }
  Direction <- if(sum(GenePeaks$log2FoldChange) > 0){
    "Up"
  } else {
    "Down"
  }
  data.frame(GeneSymbol = gene,
             TotalPeaks = TotalPeaks,
             TopList = TopList,
             NomSignif = NomSignif,
             AdjSignif = AdjSignif,
             PromoterSignif = PromoterSignif,
             UpPeaks = UpPeaks,
             DownPeaks = DownPeaks,
             Direction = Direction,
             TopP = min(GenePeaks$Adj2)) %>% mutate(NomSignifProp = round(NomSignif/TotalPeaks, digits = 2),
                                               AdjSignifProp = round(AdjSignif/TotalPeaks, digits = 2))
}, simplify = F) %>% do.call(rbind, .) 


TopGenes %<>% filter(TopList == "Yes") %>% .[!grepl("^MIR|^SNOR", .$GeneSymbol),]

save.image(paste0(ResultsPath, Cohort, ".Rdata"))
saveRDS(DESeqOutAll_Full, file = paste0(ResultsPath, Cohort, "DEoutput.Rds"))
saveRDS(DESeqOutAll_Full_noCorrection, file = paste0(ResultsPath, Cohort, "DEoutputNoCorrection.Rds"))
saveRDS(DESeqOutAll_Full_RLE_noCorrection, file = paste0(ResultsPath, Cohort, "DEoutputRLENoCorrection.Rds"))
saveRDS(DESeqOutAll_Full_RLE, file = paste0(ResultsPath, Cohort, "DEoutputRLE.Rds"))
saveRDS(DESegResultsGroup_FullAll, file = paste0(ResultsPath, Cohort, "DEresults.Rds"))
saveRDS(DESegResultsGroup_FullAll_RLE, file = paste0(ResultsPath, Cohort, "DEresultsRLE.Rds"))
saveRDS(DESegResultsGroup_FullAll_noCorrection, file = paste0(ResultsPath, Cohort, "DEresultsNoCorrection.Rds"))
saveRDS(DESegResultsGroup_FullAll_RLE_noCorrection, file = paste0(ResultsPath, Cohort, "DEresultsRLENoCorrection.Rds"))

