ResultsPath = "MPPcomparison/Results"
if(!ResultsPath %in% list.dirs(full.names = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

source("/data/Rprojects/GeneralScripts/generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")

##### Get relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
##### Get expression-based cell type marker genes ###########

if(!"homologene" %in% rownames(installed.packages())){
  install_github("oganm/homologene", force = T)
}
library(homologene)

if(!"markerGeneProfile" %in% rownames(installed.packages())){
  install_github("oganm/markerGeneProfile", force = T)
}
library(markerGeneProfile)
data("mouseMarkerGenesCombined")

CellType_genes <- mouseMarkerGenesCombined$Cortex
for(i in 1:length(CellType_genes)){
  
  CellType_genes[[i]] <- as.vector(mouse2human(CellType_genes[[i]])$humanGene)
  
}


# Get the cell type specific (neuron vs glia) peaks

if("CellTypeH3K27ac.tsv" %in% list.files("CellTypeH3K27ac")){
  CellTypePeaks <- read.table("CellTypeH3K27ac/CellTypeH3K27ac.tsv", header = T, sep = "\t")
} else {
  source(paste0("ProjectScripts/CellProportions.R"))
}

CellTypePeaks %<>% mutate(Peak_Gene = paste(PeakName, symbol, sep = "_"))
CellTypePeaks %<>% filter(!duplicated(Peak_Gene), !is.na(symbol))

#Identify cell type specific peaks by intersecting between the cell type specific genes and cell type specific peaks

CellTypeSpecificPeaks <-lapply(CellType_genes, function(CellType){
  MarkerPeaks <- CellTypePeaks %>% filter(symbol %in% CellType, !duplicated(Peak_Gene)) %>% select(PeakName, symbol, log2FoldChange, padj) %>% arrange(padj) 
  if(nrow(MarkerPeaks) > 2){
    MarkerPeaks
  } else {
    NULL
  }
})

#Add general NeuN enriched peaks
CellTypeSpecificPeaks$NeuNall <- CellTypePeaks %>% filter(log2FoldChange > 3, !duplicated(Peak_Gene)) %>% select(PeakName, symbol, log2FoldChange, padj) %>% arrange(padj)

CellTypeSpecificPeaks <- CellTypeSpecificPeaks[!unlist(lapply(CellTypeSpecificPeaks, function(x) is.null(x)))]

CellTypeSpecificPeakNames <- lapply(CellTypeSpecificPeaks, function(CellType){
  CellType %>% filter(!duplicated(PeakName)) %>% .$PeakName %>% as.character()
})


#### Get the count matrix based on cell type specific peaks
HTseqCountsSamples <- read.table("MPPcomparison/data/sun_et_al_counts.tsv", header = T, sep = "\t")

Metadata <- read.table("MPPcomparison/meta/metadata_download_samples.csv", header = T, sep = ",")
Metadata$BrainRegion <- sapply(as.character(Metadata$BrainRegion), function(x){
  if(x == "C"){
    "Cerebellum"
  } else {
    x
  }
}) %>% factor(levels = c("PFC", "Cerebellum"))

Metadata$Sample_ID <- sapply(as.character(Metadata$Sample_ID), function(x) gsub("-", ".", x))
                       
names(HTseqCountsSamples)[grepl("V|BA", names(HTseqCountsSamples))] <- sapply(names(HTseqCountsSamples)[grepl("final", names(HTseqCountsSamples))], function(x) paste0(strsplit(x, "_")[[1]][c(10,11)], collapse = "_"))

HTseqCountsSamples$NeuronFC <- CellTypePeaks$log2FoldChange[match(HTseqCountsSamples$Geneid, CellTypePeaks$PeakName)]

counMatrixSamples <- HTseqCountsSamples %>% select(as.character(Metadata$Sample_ID)) %>% as.matrix()
rownames(counMatrixSamples) <- as.character(HTseqCountsSamples$Geneid)

samplesCPM <- Count2CPM(counMatrixSamples)

## Run PCA on celltype specific peaks (+ general Neuronal peaks) ##########

PCAcellType <- lapply(CellTypeSpecificPeakNames, function(CellType){
  samplesCPM <- samplesCPM[rownames(samplesCPM) %in% CellType,]
  temp <- prcomp(t(samplesCPM), scale. = T)
  sumPC1 = sum(temp$rotation[,1])
  if(sumPC1 < 0){
    temp$rotation[,1] <- -1*temp$rotation[,1]
    temp$x[,1] <- -1*temp$x[,1]
  }
  while(sum(temp$rotation[,1] > 0) < nrow(temp$rotation)){
    samplesCPM <- samplesCPM[rownames(samplesCPM) %in% rownames(temp$rotation)[temp$rotation[,1] > 0],]
    temp <- prcomp(t(samplesCPM), scale. = T)
    sumPC1 = sum(temp$rotation[,1])
    if(sumPC1 < 0){
      temp$rotation[,1] <- -1*temp$rotation[,1]
      temp$x[,1] <- -1*temp$x[,1]
    }
  }
  temp
})

VarianceExplained <- lapply(PCAcellType, function(cell){
  temp <- cell %>% summary
  if(length(temp) > 3){
    temp$importance %>% .[2,1:3] 
  }
}) %>% do.call(rbind,.) %>% data.frame %>%
  mutate(Peaks = unlist(lapply(PCAcellType, function(x) nrow(x$rotation))))
rownames(VarianceExplained) <- names(PCAcellType)

MPP_out <- lapply(PCAcellType, function(CellType){
  rescale(CellType$x[,1], to = c(0,1))
}) %>% do.call(cbind, .) %>% data.frame %>% mutate(Sample = row.names(.))
  
#Add general calculation based on ratio between bneuronal and glial peaks, this is for comparing different regions
samplesCPM <- samplesCPM[rownames(samplesCPM) %in% (CellTypePeaks %>% filter(abs(log2FoldChange) > 3, !duplicated(PeakName)) %>% .$PeakName),]

CellTypePCA <- prcomp(t(samplesCPM), scale. = T)
#Make sure that the loadings of genes upregulated in neurons is positive
DF <- data.frame(Peak = rownames(CellTypePCA$rotation),
                 PC1 = CellTypePCA$rotation[,1],
                 NeuronFC = HTseqCountsSamples$NeuronFC[match(rownames(CellTypePCA$rotation), as.character(HTseqCountsSamples$Geneid))])
DF$NeuronalPeak <- sapply(DF$NeuronFC, function(x){
  if(x > 0){
    "NeuronalPeak"
  } else {
    "GlialPeak"
  }
})

ggplot(DF, aes(NeuronalPeak,PC1, color = NeuronalPeak)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.5)
MedPC1 <- DF %>% group_by(NeuronalPeak) %>% summarise(Median = median(PC1)) %>% data.frame()

if((MedPC1 %>% filter(MedPC1$NeuronalPeak == "NeuronalPeak") %>% .$Median) < (MedPC1 %>% filter(MedPC1$NeuronalPeak == "GlialPeak") %>% .$Median )){
  CellTypePCA$rotation[,1] <- -1*CellTypePCA$rotation[,1]
  CellTypePCA$x[,1] <- -1*CellTypePCA$x[,1]
}

Metadata$NeuronProp <- rescale(to = c(0,1), x = CellTypePCA$x[,1])

Metadata <- merge(Metadata, MPP_out, by.x = "Sample_ID", by.y = "Sample")

ggplot(Metadata, aes(age, NeuronProp, color = BrainRegion)) +
  theme_classic(base_size = 14) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = BrainRegion), size = 0.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("aquamarine4", "chocolate1")) +
  scale_fill_manual(values = c("aquamarine4", "chocolate1")) 

ggplot(Metadata, aes(age, NeuNall, color = BrainRegion)) +
  theme_classic(base_size = 14) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = BrainRegion), size = 0.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("aquamarine4", "chocolate1")) +
  scale_fill_manual(values = c("aquamarine4", "chocolate1")) +
  facet_wrap(~BrainRegion)

Metadata$SubjectID <- sapply(as.character(Metadata$Sample_ID), function(x){
  strsplit(x, "_")[[1]][1]
})

MetaMelt <- gather(Metadata, key = "CellType", value = "MSP", Astrocyte, Endothelial, Microglia, Oligo, OligoPrecursors)


Plot <- ggplot(MetaMelt, aes(SubjectID, MSP, color = BrainRegion)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.y  = element_blank(), panel.grid.minor.y  = element_blank(), axis.text.x = element_blank()) +
  labs(x = "Subject") +
  geom_point(size = 2) +
  scale_color_manual(values = c("aquamarine4", "chocolate1")) +
  scale_fill_manual(values = c("aquamarine4", "chocolate1")) +
  facet_grid(~CellType)
  
ggsave("MPPRegionComparison.pdf", plot = Plot, device = "pdf", width = 10, height = 4, dpi = 300, useDingbats = F, path = ResultsPath)
ggsave("MPPRegionComparison.png", plot = Plot, device = "png", width = 10, height = 4, dpi = 300, path = ResultsPath)
