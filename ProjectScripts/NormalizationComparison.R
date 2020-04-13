source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
plotMA = DESeq2::plotMA
packageF("org.Hs.eg.db")

ResultsPath = "NormalizationComparison"
SampleFilter =  function(x) x %>% filter(age > 20, cohort == "PV") %>% droplevels()
CombinedPeakLoc = "data/Peaks/PV_ALL_above_20.narrowPeakClean.gz"
CountMatrixLoc = "data/PV_ALL_above_20_counts.tsv.gz"
CellTypePeakCountLoc = "CellTypeH3K27ac/data/PV_ALL_above_20_counts.tsv"
Cohort = "PV_all_adults"


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
MetaChipSeqAllcalled <- MetaChipSeq

MetaChipSeqAllcalled$Peaks <- CombinedPeakLoc


################ Peak analysis ############

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

##### Get relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
AllCalledData$SampleInfo <- AllCalledData$SampleInfo <- GetCellularProportions(AllCalledData$SampleInfo, MetaSamplCol = "SampleID")


Model = as.formula(" ~ condition + sex + age + batch + pm_hours + Oligo_MSP")

countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",CorMethod = "pearson",countSampleRegEx = "^X",MetaSamleCol = "SampleID", MetaSamleIDCol = "SampleID",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))

#Filer peaks with low normalized count
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