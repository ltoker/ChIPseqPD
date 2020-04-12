ResultsPath = "GeneralResults"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")
source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
packageF("ggrepel")
packageF("metap")
packageF("data.table")

resultsPV <- readRDS("ResultsPV_Final/PVDEresults.Rds")
resultsPV_RLE <- readRDS("ResultsPV_Final/PVDEresultsRLE.Rds")

resultsNBB <- readRDS("ResultsNBB_Final/NBBDEresults.Rds")

CompareData <- sapply(c("results", "resultsNoCorrection", "resultsRLE", "resultsRLENoCorrection"), function(ModelType){
  NBBresults <- readRDS(paste0("ResultsNBB_Final/NBBDE", ModelType, ".Rds")) %>%
    arrange(pvalue) %>% mutate(UniqueRegionGeneTypePeak = paste0(UniqueRegionGeneType, PeakName)) %>%
    filter(!duplicated(UniqueRegionGeneTypePeak))
  PVresults <- readRDS(paste0("ResultsPV_Final/PVDE", ModelType, ".Rds")) %>%
    arrange(pvalue) %>% mutate(UniqueRegionGeneTypePeak = paste0(UniqueRegionGeneType, PeakName)) %>%
    filter(!duplicated(UniqueRegionGeneTypePeak))
  
  AllData <- merge(NBBresults %>% select(PeakName, baseMean, log2FoldChange, pvalue, padj, symbol, UniqueRegionGeneType),
                   PVresults %>% select(PeakName, baseMean, log2FoldChange, pvalue, padj, UniqueRegionGeneType), by = "UniqueRegionGeneType",
                   suffixes = c(".NBB", ".PV"), all.x = F, all.y = F)
  AllData$Region <- sapply(as.character(AllData$UniqueRegionGeneType), function(x){
    strsplit(x, "_")[[1]][5]
  }) %>% factor(levels = c("promoters", "exons", "5UTRs", "3UTR", "1to5kb","intronexonboundaries", "introns", "intergenic"))
  
  AllData %<>% mutate(Region2 = Region)                                 
  levels(AllData$Region2) <- c(1:8)
  
  AllData %<>% arrange(Region2)
  
  SignifBoth <- AllData %>% filter(padj.NBB < 0.05, padj.PV < 0.05) %>%
    mutate(DuplCol = paste0(baseMean.NBB, baseMean.PV, UniqueRegionGeneType), CohorSignif = "Both") %>%
    filter(!duplicated(DuplCol)) %>% select(-matches("^Regio|Dupl")) %>% filter(!duplicated(UniqueRegionGeneType))
  SignifNBB <- AllData %>% filter(padj.NBB < 0.05, padj.PV > 0.05) %>% 
    mutate(DuplCol = paste0(baseMean.NBB, baseMean.PV, UniqueRegionGeneType), CohorSignif = "NBB only") %>%
    filter(!duplicated(DuplCol)) %>% select(-matches("^Regio|Dupl")) %>% filter(!duplicated(UniqueRegionGeneType))
  SignifPV <- AllData %>% filter(padj.NBB > 0.05, padj.PV < 0.05) %>% 
    mutate(DuplCol = paste0(baseMean.NBB, baseMean.PV, UniqueRegionGeneType), CohorSignif = "PV only") %>%
    filter(!duplicated(DuplCol)) %>% select(-matches("^Regio|Dupl")) %>% filter(!duplicated(UniqueRegionGeneType))
  NonSignif <- AllData %>% filter(padj.NBB > 0.05, padj.PV > 0.05) %>% 
    mutate(DuplCol = paste(baseMean.NBB, baseMean.PV, UniqueRegionGeneType), CohorSignif = "NS") %>%
    filter(!duplicated(DuplCol)) %>% select(-matches("^Regio|Dupl")) %>% filter(!duplicated(UniqueRegionGeneType))
  Filtered <- AllData %>% filter(is.na(padj.NBB) | is.na(padj.PV)) %>% 
    mutate(DuplCol = paste(baseMean.NBB, baseMean.PV, UniqueRegionGeneType), CohorSignif = "Filtered") %>%
    filter(!duplicated(DuplCol)) %>% select(-matches("^Regio|Dupl")) %>% filter(!duplicated(UniqueRegionGeneType))
  AllData2 <- rbind(NonSignif, SignifNBB, SignifPV, SignifBoth, Filtered)
  
  #remove intergenic regions - these are too big to consider common and genes annotated as NA since they are uncomparable
  AllData2 %<>% .[!grepl("inter", .$UniqueRegionGeneType),] %>% filter(!is.na(symbol))
  
  AllData2$UniqueRegionGeneType <- sapply(as.character(AllData2$UniqueRegionGeneType), function(x){
    gsub("_hg19_genes", "", x)
  })
  

  AllData2$CohorSignif <- factor(AllData2$CohorSignif, levels = c("Filtered", "NS", "NBB only", "PV only", "Both"))
  AllData2
}, simplify = F)

GenesToMark = c("DLG2","PTPRH")

#GenesToMark = c("DLG2", "CTNND2", "RNF216", "PDE4B", "SIK3", "SLC46A3", "FRMD5", "SV2B", "PTPRH", "MKLN1")
CompareData <- lapply(CompareData, function(ModelData){
  ModelData$Mark <- sapply(as.character(ModelData$symbol), function(x){
    if(is.na(x)){
      ""
    } else if(x %in% GenesToMark){
      x
    } else {
      ""
    }
  })
  ModelData
})

#Print the plots
sapply(names(CompareData), function(ModelType){
  
  #The filtering is to remove genes filtered by DESeq IndependFiltering 
  
  AllData2 = CompareData[[ModelType]] %>% filter(CohorSignif != "Filtered") %>% droplevels()
  if(grepl("Correct", ModelType)){
    Alpha = NA
  } else {
    Alpha = 0.3
  }
  Plot <- ggplot(AllData2, aes(log2FoldChange.PV, log2FoldChange.NBB)) +
    theme_classic() +
    xlim(-1.2, 1.4) +
    ylim(-2, 2) +
    geom_point(color = "grey", size = 0.4) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(data = AllData2 %>% filter(CohorSignif %in%  c("NBB only", "PV only")),
               aes(log2FoldChange.PV, log2FoldChange.NBB, color = CohorSignif), alpha = Alpha, size = 1) +
    geom_point(data = AllData2 %>% filter(CohorSignif == "Both"),
               aes(log2FoldChange.PV, log2FoldChange.NBB, color = CohorSignif), size = 1) +
    scale_color_manual(values = c( "grey", "cornflowerblue", "darkolivegreen","coral2"),
                       name = "Significance", drop = F)
  if(grepl("Correct", ModelType)){
    Plot <- Plot +
      geom_label_repel(arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
                     nudge_y = 0.2, nudge_x = 0.2, aes(label = Mark, color = CohorSignif), show.legend = F)
  } else {
    Plot <- Plot +
      geom_density2d(data = AllData2 %>% filter(CohorSignif != "NS"),
                     aes(color = CohorSignif)) +
      geom_label_repel(arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
                       nudge_y = 0.2, nudge_x = 0.2, aes(label = Mark, color = CohorSignif), show.legend = F)
  }

  #plot(Plot) 
  ggsave(paste0(ResultsPath,"NBB_PVsignif", gsub("results", "", ModelType), ".png"), Plot, device = "png", width = 10, height = 6, dpi = 300)
  ggsave(paste0(ResultsPath,"NBB_PVsignif", gsub("results", "", ModelType), ".pdf"), Plot, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F)
  closeDev()
})


#Add information whether the genes are implicated in PD
PDgenes <- read.table("data/PD_implicated_genes.txt", header = T, sep = "\t")

#Get gene information
packageF("EnsDb.Hsapiens.v75")
packageF("ensembldb")
edb <- EnsDb.Hsapiens.v75
txdb <- as.data.frame(transcriptsBy(edb, by="gene",
                                    columns=c("gene_id", "tx_id", "gene_name", "gene_biotype")))

#Restrict to nuclear genes and genes on canonical chromosomes
txdb <- subset(txdb, seqnames %in% c(1:22,"X","Y", "MT") & startsWith(gene_id, 'ENSG'))
txdb$seqnames <- as.character(txdb$seqnames)

geneNames <- txdb %>% select(gene_name, gene_id, gene_biotype, width)
names(geneNames) <- c("hgnc_symbol", "ensembl_gene_id", "gene_biotype", "GeneLength")


#Filter miRNA and sn/snoRNAs
symbolToFilter <- geneNames %>% filter(gene_biotype %in% c("miRNA", "snRNA", "snoRNA")) %>%
  .$hgnc_symbol %>% as.character()

CommonRegions <- CompareData$results
CommonRegions %<>% filter(!symbol %in% symbolToFilter, CohorSignif != "Filtered")


#Add whether the gene is implicated in PD
CommonRegions$PDgene <- sapply(CommonRegions$symbol, function(gene){
  if(as.character(gene) %in% as.character(PDgenes$Gene)){
    "Yes"
  } else {
    "No"
  }
})


CommonRegions$metaP <- apply(CommonRegions %>% select(matches("pvalue")), 1,  function(region){
  metap::sumlog(region)$p
})

CommonRegions$GWAS <- "No"
CommonRegions$GWAS[CommonRegions$symbol %in% as.character(PDgenes %>%
                                                          filter(GWAS_Nalls2019 == "YES") %>% .$Gene)] <- "Yes"

CommonRegions %<>% mutate(CommonPair = paste0("NBB.", PeakName.NBB, "PV.", PeakName.PV))

#Adjusting metaP for the common peaks. The adjsutment is for the specific CommonPeaks pair 
AdjPvalues <- data.frame(metaP = CommonRegions %>%filter(!duplicated(CommonPair)) %>% .$metaP)
AdjPvalues$AdjmetaP <- p.adjust(AdjPvalues$metaP, method = "BH")
AdjPvalues %<>% filter(!duplicated(metaP))

CommonRegions <- merge(CommonRegions, AdjPvalues, by = "metaP")
CommonRegions$metaP <- sapply(CommonRegions$metaP, function(x) signif(x, digits = 2))


PeakNum <- CommonRegions %>%
  group_by(symbol) %>% summarise(n = n()) %>% data.frame

CommonRegions$PeakNum <- PeakNum$n[match(CommonRegions$symbol, PeakNum$symbol)]

CommonRegions$RegionLength <- sapply(CommonRegions$UniqueRegionGeneType, function(x){
  x = strsplit(x, "_")[[1]][1]
  x =  strsplit(x,":")[[1]][2]
  x = strsplit(x,"-")[[1]] %>% as.numeric()
  x[2]-x[1]
})


CommonEffectiveLength <- CommonRegions %>%
  group_by(symbol) %>% summarise(Length = sum(as.numeric(RegionLength))) %>% data.frame()

CommonRegions$EffectiveLength <- CommonEffectiveLength$Length[match(CommonRegions$symbol, CommonEffectiveLength$symbol)]

MetaPThresh = CommonRegions %>% filter(AdjmetaP < 0.05) %>% .$metaP %>% max

CommonRegions$DirectionChange <- apply(CommonRegions %>% select(matches("log2")), 1, function(x){
  if(x[1] > 0 & x[2] > 0){
    "Hyperacetylated"
  } else if(x[1] < 0 & x[2] < 0){
    "Hypoacetylated"
  } else {
    "Mixed"
  }
}) %>% factor(levels = c("Hypoacetylated", "Hyperacetylated", "Mixed"))



Plot <- ggplot(CommonRegions %>% filter(!duplicated(CommonPair)), aes(log10(EffectiveLength), -log10(metaP))) +
  theme_classic(14) +
  geom_point(alpha = 0.2, aes(color = PDgene), shape = 16) +
  geom_point(data = CommonRegions %>% filter(PDgene == "Yes") %>% filter(!duplicated(CommonPair)),
             aes(log10(EffectiveLength), -log10(metaP)), color = "goldenrod1") +
  geom_point(data = CommonRegions %>% filter(CohorSignif != "NS", PDgene == "Yes") %>% filter(!duplicated(CommonPair)),
             aes(log10(EffectiveLength), -log10(metaP)), shape = 1, size = 2 , color = "red") +
  scale_color_manual(values = c("black", "goldenrod1"), name = "PDgene") +
  geom_hline(yintercept = -log10(MetaPThresh),
             color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~DirectionChange, nrow = 1)


ggsave(paste0(ResultsPath,"EffectiveGeneLengt.pdf"), 
       device = "pdf", Plot, width = 10, height = 4, dpi = 300, useDingbats = F)


lmMetahypo <- lm(-log10(metaP)~PDgene+log10(EffectiveLength) + PeakNum,
                 data = CommonRegions %>% filter(DirectionChange == "Hypoacetylated"))

lmMetahyper <- lm(-log10(metaP)~PDgene+log10(EffectiveLength) + PeakNum,
                  data = CommonRegions %>% filter(DirectionChange == "Hyperacetylated"))

ConfIntMeta <- rbind(confint.lm(lmMetahypo) %>% data.frame() %>% .[-1,] %>%
                       mutate(beta = summary(lmMetahypo)$coef[,1][-1], Covar = rownames(.), Direction = "Hypoacetylated"),
                     confint.lm(lmMetahyper) %>% data.frame() %>% .[-1,] %>%
                       mutate(beta = summary(lmMetahyper)$coef[,1][-1], Covar = rownames(.), Direction = "Hypercetylated"))
names(ConfIntMeta)[1:2] <- c("low2.5", "high97.5")
ConfIntMeta$Cohort  = "MetaP"


#Repeat for PV only
resultsPV %<>% filter(!symbol %in% symbolToFilter)
resultsPV$GeneLength <- geneNames$GeneLength[match(resultsPV$symbol, geneNames$hgnc_symbol)]
resultsPV$PDgene <- sapply(resultsPV$symbol, function(gene){
  if(gene %in% as.character(PDgenes$Gene)){
    "Yes"
  } else {
    "No"
  }
})

PeakNumPV <- resultsPV %>% filter(!is.na(symbol)) %>%
  group_by(symbol) %>% summarise(n = n(), num = mean(GeneLength)) %>% data.frame

GeneEffectiveLengthPV <- resultsPV %>% filter(!is.na(symbol)) %>%
  group_by(symbol) %>% summarise(Length = sum(as.numeric(Region.width))) %>% data.frame()

resultsPV$PeakNum <- PeakNumPV$n[match(resultsPV$symbol, PeakNumPV$symbol)]
resultsPV$EffectiveLength <- GeneEffectiveLengthPV$Length[match(resultsPV$symbol, GeneEffectiveLengthPV$symbol)]
resultsPV$DirectionChange <- sapply(resultsPV$log2FoldChange, function(x){
  if(x > 0){
    "Hyperacetylated"
  } else {
    "Hypoacetylated"
  }
}) %>% factor(levels = c("Hypoacetylated", "Hyperacetylated"))


lmPVhypo <- lm(-log10(pvalue)~PDgene+log10(EffectiveLength) + PeakNum,
               data = resultsPV %>% filter(DirectionChange == "Hypoacetylated"))

lmPVhyper <- lm(-log10(pvalue)~PDgene+log10(EffectiveLength) + PeakNum,
                data = resultsPV %>% filter(DirectionChange == "Hyperacetylated"))

ConfIntPV <- rbind(confint.lm(lmPVhypo) %>% data.frame() %>% .[-1,] %>%
                     mutate(beta = summary(lmPVhypo)$coef[,1][-1], Covar = rownames(.), Direction = "Hypoacetylated"),
                   confint.lm(lmPVhyper) %>% data.frame() %>% .[-1,] %>%
                     mutate(beta = summary(lmPVhyper)$coef[,1][-1], Covar = rownames(.), Direction = "Hypercetylated"))
names(ConfIntPV)[1:2] <- c("low2.5", "high97.5")
ConfIntPV$Cohort  = "PW"


#Repeat for NBB
resultsNBB %<>% filter(!symbol %in% symbolToFilter)
resultsNBB$GeneLength <- geneNames$GeneLength[match(resultsNBB$symbol, geneNames$hgnc_symbol)]
resultsNBB$PDgene <- sapply(resultsNBB$symbol, function(gene){
  if(gene %in% as.character(PDgenes$Gene)){
    "Yes"
  } else {
    "No"
  }
})

PeakNumNBB <- resultsNBB %>% filter(!is.na(symbol)) %>%
  group_by(symbol) %>% summarise(n = n(), num = mean(GeneLength)) %>% data.frame

GeneEffectiveLengthNBB <- resultsNBB %>% filter(!is.na(symbol)) %>%
  group_by(symbol) %>% summarise(Length = sum(as.numeric(Region.width))) %>% data.frame()

resultsNBB$PeakNum <- PeakNumNBB$n[match(resultsNBB$symbol, PeakNumNBB$symbol)]
resultsNBB$EffectiveLength <- GeneEffectiveLengthNBB$Length[match(resultsNBB$symbol, GeneEffectiveLengthNBB$symbol)]
resultsNBB$DirectionChange <- sapply(resultsNBB$log2FoldChange, function(x){
  if(x > 0){
    "Hyperacetylated"
  } else {
    "Hypoacetylated"
  }
}) %>% factor(levels = c("Hypoacetylated", "Hyperacetylated"))


lmNBBhypo <- lm(-log10(pvalue)~PDgene+log10(EffectiveLength) + PeakNum,
                data = resultsNBB %>% filter(DirectionChange == "Hypoacetylated"))

lmNBBhyper <- lm(-log10(pvalue)~PDgene+log10(EffectiveLength) + PeakNum,
                 data = resultsNBB %>% filter(DirectionChange == "Hyperacetylated"))

ConfIntNBB <- rbind(confint.lm(lmNBBhypo) %>% data.frame() %>% .[-1,] %>%
                      mutate(beta = summary(lmNBBhypo)$coef[,1][-1], Covar = rownames(.), Direction = "Hypoacetylated"),
                    confint.lm(lmNBBhyper) %>% data.frame() %>% .[-1,] %>%
                      mutate(beta = summary(lmNBBhyper)$coef[,1][-1], Covar = rownames(.), Direction = "Hypercetylated"))
names(ConfIntNBB)[1:2] <- c("low2.5", "high97.5")
ConfIntNBB$Cohort  = "NBB"

ConfIntAll <- rbind(ConfIntMeta, ConfIntPV, ConfIntNBB)
ConfIntAll %<>% mutate(Covar = factor(Covar, levels = c("PDgeneYes","log10(EffectiveLength)","PeakNum")),
                       Cohort = as.factor(Cohort),
                       Direction = as.factor(Direction))

levels(ConfIntAll$Covar) <- c("PDgene", "log10(EGL)", "TotalPeaks")
ConfIntAll$Covar <- relevel(ConfIntAll$Covar, ref = "PDgene")
ConfIntAll$Direction <- relevel(ConfIntAll$Direction, ref = "Hypoacetylated")
ConfIntAll$Cohort <- factor(ConfIntAll$Cohort, levels = c("MetaP", "PW", "NBB"))

Plot <-   ggplot(ConfIntAll %>% filter(Cohort == "MetaP"), aes(Covar, beta, color = Covar)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "Coefficient") +
  geom_hline(yintercept = 0, col  = "darkgrey", linetype = "dashed", size = 0.5) +
  geom_point() +
  geom_errorbar(aes(ymin = low2.5, ymax = high97.5, color = Covar), size = 1) +
  scale_color_manual(values = c("brown4" ,"darkcyan", "blue4"), name = "Covariate") +
  facet_wrap(Cohort~Direction, nrow = 1)

ggsave(paste0(ResultsPath,"EffectiveGeneLengtMetaCI.pdf"), 
       device = "pdf", Plot, width = 5, height = 2.5, dpi = 300, useDingbats = F)

Plot <-   ggplot(ConfIntAll %>% filter(Cohort == "PW"), aes(Covar, beta, color = Covar)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "Coefficient") +
  geom_hline(yintercept = 0, col  = "darkgrey", linetype = "dashed", size = 0.5) +
  geom_point() +
  geom_errorbar(aes(ymin = low2.5, ymax = high97.5, color = Covar), size = 1) +
  scale_color_manual(values = c("brown4" ,"darkcyan", "blue4"), name = "Covariate") +
  facet_wrap(Cohort~Direction, nrow = 1)

ggsave(paste0(ResultsPath,"EffectiveGeneLengtPVCI.pdf"), 
       device = "pdf", Plot, width = 5, height = 2.5, dpi = 300, useDingbats = F)


Plot <-   ggplot(ConfIntAll %>% filter(Cohort == "NBB"), aes(Covar, beta, color = Covar)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "", y = "Coefficient") +
  geom_hline(yintercept = 0, col  = "darkgrey", linetype = "dashed", size = 0.5) +
  geom_point() +
  geom_errorbar(aes(ymin = low2.5, ymax = high97.5, color = Covar), size = 1) +
  scale_color_manual(values = c("brown4" ,"darkcyan", "blue4"), name = "Covariate") +
  facet_wrap(Cohort~Direction, nrow = 1)

ggsave(paste0(ResultsPath,"EffectiveGeneLengtNBBCI.pdf"), 
       device = "pdf", Plot, width = 5, height = 2.5, dpi = 300, useDingbats = F)


SummaryLM <- sapply(ls(pat = "lm.*hyp"), function(x){
  eval(as.name(x)) %>% summary
})



pValLengthDF <- data.frame(Type = c("Hypoacetylated", "Hyperacetylated"),
                           CommonRegionsPDgene = c(paste0("beta = ", round(SummaryLM$lmMetahypo$coef[2,1],digits = 2),
                                                        ", p = ", signif(SummaryLM$lmMetahypo$coef[2,4], digits = 1)),
                                                 paste0("beta = ", round(SummaryLM$lmMetahyper$coef[2,1],digits = 2),
                                                        ", p = ", signif(SummaryLM$lmMetahyper$coef[2,4], digits = 1))),
                           CommonRegionsEGL = c(paste0("beta = ", round(SummaryLM$lmMetahypo$coef[3,1],digits = 2),
                                                       ", p = ", signif(SummaryLM$lmMetahypo$coef[3,4], digits = 1)),
                                                paste0("beta = ", round(SummaryLM$lmMetahyper$coef[3,1],digits = 2),
                                                       ", p = ", signif(SummaryLM$lmMetahyper$coe[3,4], digits = 1))),
                           PVPDgene = c(paste0("beta = ", round(SummaryLM$lmPVhypo$coef[2,1],digits = 2),
                                             ", p = ", signif(SummaryLM$lmPVhypo$coef[2,4], digits = 1)),
                                      paste0("beta = ", round(SummaryLM$lmPVhyper$coef[2,1],digits = 2),
                                             ", p = ", signif(SummaryLM$lmPVhyper$coef[2,4], digits = 1))),
                           PVEGL = c(paste0("beta = ", round(SummaryLM$lmPVhypo$coefficients[3,1],digits = 2),
                                            ", p = ", signif(SummaryLM$lmPVhypo$coef[3,4], digits = 1)),
                                     paste0("beta = ", round(SummaryLM$lmPVhyper$coef[3,1],digits = 2),
                                            ", p = ", signif(SummaryLM$lmPVhyper$coef[3,4], digits = 1))),
                           NBBPDgene = c(paste0("beta = ", round(SummaryLM$lmNBBhypo$coef[2,1],digits = 2),
                                              ", p = ", signif(SummaryLM$lmNBBhypo$coef[2,4], digits = 1)),
                                       paste0("beta = ", round(SummaryLM$lmNBBhyper$coef[2,1],digits = 2),
                                              ", p = ", signif(SummaryLM$lmNBBhyper$coef[2,4], digits = 1))),
                           NBBEGL = c(paste0("beta = ", round(SummaryLM$lmNBBhypo$coef[3,1],digits = 2),
                                             ", p = ", signif(SummaryLM$lmNBBhypo$coef[3,4], digits = 1)),
                                      paste0("beta = ", round(SummaryLM$lmNBBhyper$coef[3,1],digits = 2),
                                             ", p = ", signif(SummaryLM$lmNBBhyper$coef[3,4], digits = 1))))

write.table(pValLengthDF, file = paste0(ResultsPath, "pValLengthDF.tsv"), sep = "\t", row.names = F, col.names = T)

write.table(resultsPV %>%
              select(symbol, PeakName, baseMean, Peak.width,
                     log2FoldChange, pvalue, padj, type, UniqueRegion, PDgene, PeakNum, EffectiveLength),
            file = paste0(ResultsPath, "ChIPresultsPW.tsv"), sep = "\t", row.names = F, col.names = T)

write.table(resultsNBB %>%
              select(symbol, PeakName, baseMean, Peak.width,
                     log2FoldChange, pvalue, padj, type, UniqueRegion, PDgene, PeakNum, EffectiveLength),
            file = paste0(ResultsPath, "ChIPresultsNBB.tsv"), sep = "\t", row.names = F, col.names = T)

write.table(CommonRegions %>%
              select(symbol, PeakName.PV, PeakName.NBB,
                     log2FoldChange.PV, log2FoldChange.NBB,
                     pvalue.PV, pvalue.NBB,
                     padj.PV, padj.NBB, metaP, AdjmetaP,CohorSignif, CommonPair,
                     UniqueRegionGeneType, PDgene, PeakNum, EffectiveLength),
            file = paste0(ResultsPath, "CommonRegions.tsv"), sep = "\t", row.names = F, col.names = T)


#Hypergeometric test for gene overlap
#Genes are excluded if they were filtered out based on DESeq2 independent filtering

AllGenesPV <- resultsPV %>% arrange(padj) %>% filter(!is.na(symbol), !is.na(padj), PDgene == "No", !duplicated(symbol))
AllGenesNBB <- resultsNBB %>% arrange(padj) %>% filter(!is.na(symbol), !is.na(padj), PDgene == "No", !duplicated(symbol))
AllGenesOverlap <- intersect(AllGenesPV$symbol, AllGenesNBB$symbol)

SignifGenesPVdown <- resultsPV %>% arrange(padj) %>% filter(!is.na(symbol), padj < 0.05, log2FoldChange < 0, !duplicated(symbol))
SignifGenesPVup <- resultsPV %>% arrange(padj) %>% filter(!is.na(symbol), padj < 0.05, log2FoldChange > 0, !duplicated(symbol))

SignifGenesNBBdown <- resultsNBB %>% arrange(padj) %>% filter(!is.na(symbol), padj < 0.05, log2FoldChange < 0) %>%  filter(!duplicated(symbol))
SignifGenesNBBup <- resultsNBB %>% arrange(padj) %>% filter(!is.na(symbol), padj < 0.05, log2FoldChange > 0) %>%  filter(!duplicated(symbol))

SignifGenesOverlapDown <- intersect(SignifGenesPVdown$symbol, SignifGenesNBBdown$symbol)
SignifGenesOverlapUp <- intersect(SignifGenesPVup$symbol, SignifGenesNBBup$symbol)

SignifGenesOverlap <- c(SignifGenesOverlapDown, SignifGenesOverlapUp)

SignifGenesOverlapDFPV <- resultsPV %>% filter(symbol %in% SignifGenesOverlap) %>%
  arrange(padj) %>% filter(!duplicated(symbol)) %>% select(symbol, log2FoldChange, padj) %>% arrange(log2FoldChange)
SignifGenesOverlapDFNBB <- resultsNBB %>% filter(symbol %in% SignifGenesOverlap) %>%
  arrange(padj) %>% filter(!duplicated(symbol)) %>% select(symbol, log2FoldChange, padj) %>% arrange(log2FoldChange)
SignifGenesOverlapDF <- merge(SignifGenesOverlapDFPV, SignifGenesOverlapDFNBB, by = "symbol", suffixes = c("_PW", "_NBB")) %>% arrange(padj_PW)
SignifGenesOverlapDF %<>% mutate_if(is.numeric, function(x) signif(x, digits = 2))

write.table(SignifGenesOverlapDF, paste0(ResultsPath, "ReplicatedGenes.tsv"), sep = "\t", row.names = F, col.names = T)

#Hypoacetylated genes
dhyper(x = length(SignifGenesOverlapDown),  #number of genes hypoacetylated DARs in both cohorts
       m = nrow(SignifGenesPVdown), #number of genes with hypoacetylated DARs in PV cohort
       n = length(AllGenesOverlap) - length(SignifGenesPVdown), #number of genes represented by ChIP data in both cohorts exluding genes with hypoacetylated DARs in PV cohort, which are not NAs, miRNAs and snRNAs
       k = nrow(SignifGenesNBBdown)) #number of genes with hypoacetylated DARs in NBB cohort

#Hyperacetylated genes
dhyper(x = length(SignifGenesOverlapUp),  #number of genes hypoacetylated DARs in both cohorts
       m = nrow(SignifGenesPVup), #number of genes with hypoacetylated DARs in PV cohort
       n = length(AllGenesOverlap) - length(SignifGenesPVdown), #number of genes represented by ChIP data in both cohorts exluding genes with hypoacetylated DARs in PV cohort, which are not NAs, miRNAs and snRNAs
       k = nrow(SignifGenesNBBup)) #number of genes with hypoacetylated DARs in NBB cohort

##Hypergeometric test for DAR overlap.
#Genes are excluded if they were filtered out based on DESeq2 independent filtering
dhyper(x = nrow(CommonRegions %>% filter(CohorSignif == "Both") %>% filter(!(duplicated(CommonPair)))),
       m = nrow(CommonRegions %>% filter(CohorSignif %in%  c("PV only", "Both"))%>% filter(!(duplicated(CommonPair)))),
       n = nrow(CommonRegions %>% filter(!CohorSignif %in%  c("PV only", "Both")) %>% filter(!(duplicated(CommonPair)))),
       k = nrow(CommonRegions %>% filter(CohorSignif %in%  c("NBB only", "Both")) %>% filter(!(duplicated(CommonPair)))))


#Hypergeomettric test for enrichment of PD implicated genes:
PDallGenesPV <- resultsPV %>% arrange(padj) %>% filter(!is.na(symbol), !is.na(padj), PDgene == "Yes", !duplicated(symbol))
PDallGenesNBB <- resultsNBB %>% arrange(padj) %>% filter(!is.na(symbol), !is.na(padj), PDgene == "Yes", !duplicated(symbol))
PDallOverlap <- intersect(PDallGenesPV$symbol, PDallGenesNBB$symbol)


PDsignifGenesPV <- PDallGenesPV %>% filter(padj < 0.05)
PDsignifGenesNBB <- PDallGenesNBB %>% filter(padj < 0.05)
PDsignifOverlap <- intersect(PDsignifGenesPV$symbol, PDsignifGenesNBB$symbol)

dhyper(x = length(PDsignifOverlap),  #number of PD hits with DARs in both cohorts
       m = length(PDallOverlap), #number of PD hits appearing in both cohorts
       n = length(AllGenesOverlap), #number of genes represented by ChIP data in both cohorts which are not PDgene, NAs, miRNAs and snRNAs
       k = length(SignifGenesOverlap)) #number of genes with DARs in the same direction in both cohorts


#Hypergeomettric test for enrichment of GWAS in based on metaP:
CommonPDgeneAll <- CommonRegions %>% arrange(metaP) %>% filter(PDgene == "Yes", CohorSignif != "Filtered") %>%  filter(!duplicated(symbol))
CommonPDgene_DARs <- CommonRegions %>% arrange(metaP) %>% filter(PDgene == "Yes", AdjmetaP < 0.05) %>%  filter(!duplicated(symbol))
AllCommonPeakGenes <- CommonRegions %>% arrange(metaP) %>% filter(PDgene == "No") %>%  filter(!duplicated(symbol))
SignifCommonPeakGenes <- CommonRegions %>% arrange(metaP) %>% filter(AdjmetaP < 0.05) %>%  filter(!duplicated(symbol))

dhyper(x = nrow(CommonPDgene_DARs),  #number of PD genes with common DARs in both cohorts
       m = nrow(CommonPDgeneAll), #number of PD genes common regions in both cohortsn
       n = nrow(AllCommonPeakGenes), #number of genes represented by ChIP data in both cohorts which are not PD genes, NAs, miRNAs and snRNAs
       k =  nrow(SignifCommonPeakGenes)) #number of genes with common regions with adjusted metaP < 0.05 

#Venn diagrams
packageF("RVenn")
packageF("VennDiagram")
packageF("RAM")

DFsignif <- list()
DFsignif$AllGenes_PW <- resultsPV %>% filter(!is.na(padj)) %>% filter(!duplicated(symbol), !is.na(symbol)) %>% .$symbol
DFsignif$HypoAc_PW <- resultsPV %>% filter(padj < 0.05, !is.na(symbol), log2FoldChange < 0) %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$HyperAc_PW <- resultsPV %>% filter(padj < 0.05, !is.na(symbol), log2FoldChange > 0) %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$AllPDGenes_PW <- resultsPV %>% filter(PDgene == "Yes", !is.na(padj)) %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$SignifPDGenes_PW <- resultsPV %>% filter(padj < 0.05, PDgene == "Yes") %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$AllGenes_NBB <- resultsNBB %>% filter(!is.na(padj)) %>% filter(!duplicated(symbol), !is.na(symbol)) %>% .$symbol
DFsignif$HypoAc_NBB <- resultsNBB %>% filter(padj < 0.05, !is.na(symbol), log2FoldChange < 0) %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$HyperAc_NBB <- resultsNBB %>% filter(padj < 0.05, !is.na(symbol), log2FoldChange > 0) %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$SignifPDGenes_NBB <- resultsNBB %>% filter(padj < 0.05, PDgene == "Yes") %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$AllPDGenes_NBB <- resultsNBB %>% filter(PDgene == "Yes", !is.na(padj)) %>% filter(!duplicated(symbol)) %>% .$symbol
DFsignif$AllCommonGenes <- intersect(DFsignif$AllGenes_PW, DFsignif$AllGenes_NBB)
DFsignif$AllCommonPDGenes <- intersect(DFsignif$AllPDGenes_PW, DFsignif$AllPDGenes_PW)

group.venn(vectors=DFsignif[c("SignifPDGenes_PW", "SignifPDGenes_NBB")], label= T,
           fill = c("blue", "red"),
           cat.pos = c(-20, 15),
           cat.dist = c(0.1,0.1),cat.cex = 1.2,
           lab.cex=1,
           file = paste0(ResultsPath, "VennDiagramSignifPDgenes"), ext = "pdf", height = 3, width = 4)

#Venn digram overlap replicated hypoacetylated genes
group.venn(vectors=DFsignif[c("HypoAc_PW", "HypoAc_NBB",  "AllCommonGenes")], label= F,
           fill = c("darkgreen", "firebrick4", "white"),
           cat.pos = c(-160, 0,  160),
           cat.dist = c(0.01,0.01, 0.01),cat.cex = 1,
           lab.cex=1,
           file = paste0(ResultsPath, "VennDiagramReplicatedGenesDown"), ext = "pdf", height = 3, width = 4)

#Venn digram overlap replicated hyperacetylated genes
group.venn(vectors=DFsignif[c("HyperAc_PW", "HyperAc_NBB",  "AllCommonGenes")], label= F,
           fill = c("darkgreen", "firebrick4", "white"),
           cat.pos = c(-160, 160, 0),
           cat.dist = c(0.01,0.01, 0.01),cat.cex = 1,
           lab.cex=1,
           file = paste0(ResultsPath, "VennDiagramReplicatedGenesUp"), ext = "pdf", height = 3, width = 4)

group.venn(vectors=DFsignif[c("HyperAc_PW", "HyperAc_NBB",  "HypoAc_PW", "HypoAc_NBB")], label= F,
           #fill = c("darkgreen", "firebrick4", "white"),
           #cat.pos = c(-160, 160, 0),
           #cat.dist = c(0.01,0.01, 0.01),cat.cex = 1,
           lab.cex=1)

group.venn(vectors=DFsignif[c("AllGenes_PW", "AllGenes_NBB")], label= F,
           fill = c("blue", "red"),
           cat.pos = c(-15, 10),
           cat.dist = c(0.1,0.1),cat.cex = 1.2,
           lab.cex=1)

group.venn(vectors=DFsignif[c("AllPDGenes_PW", "AllPDGenes_NBB")], label= F,
           fill = c("blue", "red"),
           cat.pos = c(-20, 16),
           cat.dist = c(0.1,0.1),cat.cex = 1.2,
           lab.cex=1)

DFcommonRegions <- list()
DFcommonRegions$AllGenes_CommonPeaks <- CommonRegions %>% filter(!duplicated(symbol)) %>% .$symbol
DFcommonRegions$MetaSignif_Genes <- CommonRegions %>% filter(AdjmetaP < 0.05) %>% filter(!duplicated(symbol)) %>% .$symbol
DFcommonRegions$PD_Genes <- CommonRegions %>% filter(PDgene == "Yes") %>% filter(!duplicated(symbol)) %>% .$symbol
DFcommonRegions$CommonPairs <- CommonRegions %>% filter(!duplicated(CommonPair)) %>% .$CommonPair
DFcommonRegions$DAR_PW <- CommonRegions %>% filter(padj.PV < 0.05) %>% filter(!duplicated(CommonPair)) %>% .$CommonPair
DFcommonRegions$DAR_NBB <- CommonRegions %>% filter(padj.NBB < 0.05) %>% filter(!duplicated(CommonPair)) %>% .$CommonPair  
DFcommonRegions$Gene_PW <- CommonRegions %>% filter(CohorSignif < 0.05) %>% filter(!duplicated(symbol)) %>% .$CommonPair
DFcommonRegions$Gene_NBB <- CommonRegions %>% filter(padj.NBB < 0.05) %>% filter(!duplicated(symbol)) %>% .$symbol  


group.venn(vectors=DFcommonRegions[c("MetaSignif_Genes", "PD_Genes", "AllGenes_CommonPeaks")], label= F,
           fill = c("darkgreen", "firebrick4", "white"),
           cat.pos = c(-160, 0, 160),
           cat.dist = c(0.01,0.01, 0.01),cat.cex = 1,
           lab.cex=1,
           file = paste0(ResultsPath, "VennDiagramMetaP_PDgenes"), ext = "pdf", height = 3, width = 4)

group.venn(vectors=DFcommonRegions[c("DAR_PW", "DAR_NBB", "CommonPairs")], label= F,
           fill = c("blue", "red", "white"),
           cat.pos = c(-170, 0, 170),
           cat.dist = c(0.01,0.01, 0.01),cat.cex = 1,
           lab.cex=1,
           file = paste0(ResultsPath, "VennDiagramReplicatedDARs"), ext = "pdf", height = 3, width = 4)

#Add information which of the PD genes are significan
PDgenes$AdjMetaP <- sapply(PDgenes$Gene, function(gene){
  temp <- CommonRegions %>% filter(symbol == gene) %>% arrange(metaP) %>%
    filter(!duplicated(symbol)) %>% .$AdjmetaP
  if(length(temp) == 0){
    NA
  } else if(temp < 0.0099){
    sprintf("%0.1e", temp)
  } else {
    sprintf("%0.2f",temp)
  }
}) %>% unlist

PDgenes$PW_AdjPval <- sapply(PDgenes$Gene, function(gene){
  temp <- resultsPV %>% filter(symbol == gene) %>% arrange(padj) %>%
    filter(!duplicated(symbol)) %>% .$padj
  if(length(temp) == 0){
    NA
  } else if(is.na(temp)){
    NA
  } else if(temp < 0.0099){
    sprintf("%0.1e", temp)
  } else {
    sprintf("%0.2f",temp)
  }
}) %>% unlist

PDgenes$NBB_AdjPval <- sapply(PDgenes$Gene, function(gene){
  temp <- resultsNBB %>% filter(symbol == gene) %>% arrange(padj) %>%
    filter(!duplicated(symbol)) %>% .$padj
  if(length(temp) == 0){
    NA
  } else if(is.na(temp)){
    NA
  } else if(temp < 0.0099){
    sprintf("%0.1e", temp)
  } else {
    sprintf("%0.2f",temp)
  }
}) %>% unlist

write.table(PDgenes, paste0(ResultsPath, "PDgeneStat.tsv"), sep = "\t", row.names = F, col.names = T)

CommonPeaksPDgenes <- CompareData$results %>% filter(CohorSignif == "Both", log2FoldChange.NBB*log2FoldChange.PV > 0) %>%
  filter(symbol %in% PDgenes$Gene) %>%
  .$symbol %>% unique()



CommonRegionsUP <- CompareData$results %>% filter(CohorSignif == "Both", log2FoldChange.NBB > 0, log2FoldChange.PV > 0)
CommonRegionsUP$metaP <- apply(CommonRegionsUP %>% select(matches("pvalue")), 1,  function(region){
  metap::sumlog(region)$p %>% signif(digits = 2)
})

TopCommonRegionsUP <- CommonRegionsUP %>% arrange(metaP) %>%
  select(symbol, PeakName.NBB, PeakName.PV, log2FoldChange.NBB, log2FoldChange.PV, UniqueRegionGeneType, metaP) %>%
  filter(!duplicated(symbol)) %>% head(15)

PDgenesTop <- CommonRegionsUP %>% arrange(metaP) %>% filter(!duplicated(symbol)) %>%
  filter(symbol %in% PDgenes$Gene) %>%
  select(symbol, PeakName.NBB, PeakName.PV, log2FoldChange.NBB, log2FoldChange.PV, UniqueRegionGeneType, metaP)

CommonRegionsDown <- CompareData$results %>% filter(CohorSignif == "Both", log2FoldChange.NBB < 0, log2FoldChange.PV < 0)
CommonRegionsDown$metaP <- apply(CommonRegionsDown %>% select(matches("pvalue")), 1,  function(region){
  metap::sumlog(region)$p %>% signif(digits = 2)
})

TopCommonRegionsDown <- CommonRegionsDown %>% arrange(metaP) %>%
  select(symbol, PeakName.NBB, PeakName.PV, log2FoldChange.NBB, log2FoldChange.PV, UniqueRegionGeneType, metaP) %>%
  filter(!duplicated(symbol)) %>% head(15)

TopRegions <- rbind(TopCommonRegionsDown, TopCommonRegionsUP, PDgenesTop) %>% filter(!duplicated(symbol)) 


TopRegionsInfo <- sapply(TopRegions$symbol, function(Gene){
  dataPV <- resultsPV %>% filter(symbol == Gene)
  dataPV$UniqueRegionGeneType <- sapply(dataPV$UniqueRegionGeneType, function(x){
    gsub("_hg19_genes", "", x)
  })
  dataNBB <- resultsNBB %>% filter(symbol == Gene)
  dataNBB$UniqueRegionGeneType <- sapply(dataNBB$UniqueRegionGeneType, function(x){
    gsub("_hg19_genes", "", x)
  })
  PeakNameNBB = TopRegions %>% filter(symbol == Gene) %>% .$PeakName.NBB
  PeakNamePV = TopRegions %>% filter(symbol == Gene) %>% .$PeakName.PV
  RegionGeneType = TopRegions %>% filter(symbol == Gene) %>% .$UniqueRegionGeneType
  Data = CompareData$results %>% filter(symbol == Gene)
  TotalRegions = nrow(Data)
  SignifRegions = Data %>% filter(CohorSignif != "NS") %>% nrow
  
  if(sum(grepl("promoter", Data %>% filter(CohorSignif == "NS") %>% .$UniqueRegionGeneType)) > 0){
    PromNS = "p"
  } else {
    PromNS = ""
  }
  
  if(sum(grepl("promoter", Data %>% filter(CohorSignif != "NS") %>% .$UniqueRegionGeneType)) > 0){
    PromSig = "p"
  } else {
    PromSig  = ""
  }
  
  PeakLocationNBB = dataNBB %>% filter(PeakName == PeakNameNBB, UniqueRegionGeneType == RegionGeneType) %>%
    .$Peak.Location 
  PeakLocationPW = dataPV %>% filter(PeakName == PeakNamePV, UniqueRegionGeneType == RegionGeneType) %>%
    .$Peak.Location
  data.frame(Gene = Gene,
             PropSignif = paste0(SignifRegions, PromSig, "/", TotalRegions, PromNS, "(", signif(SignifRegions/TotalRegions, digits = 1), ")"),
             Chr = strsplit(PeakLocationPW, ":")[[1]][1] %>% gsub("chr", "", .),
             PeakLocationPW = strsplit(PeakLocationPW, ":")[[1]][2],
             PeakLocationNBB = strsplit(PeakLocationNBB, ":")[[1]][2],
             RegionAnnotation = rev(strsplit(RegionGeneType, "_")[[1]])[1])
}, simplify = F) %>% do.call(rbind, .)

TopRegionsFull <- merge(TopRegions, TopRegionsInfo, by.x = "symbol", by.y = "Gene", sort = FALSE)
TopRegionsFull$RegionLocation <- sapply(TopRegionsFull$UniqueRegionGeneType, function(x){
  x = gsub("chr[0-9]+:", "", x)
  strsplit(x, "_")[[1]][1]
})

TopRegionsFull %<>% select(-UniqueRegionGeneType)
names(TopRegionsFull) <- c("Symbol", "Peak_NBB", "Peak_PW", "LFC_NBB", "LFC_PW", "MetaP", "PropSignif", "CHR", "PeakLocation_PW", "PeakLocation_NBB", "RegionAnnot", "RegionLocation")
TopRegionsFull %<>% mutate(LFC_NBB = round(LFC_NBB, digits = 2),
                           LFC_PW = round(LFC_PW, digits = 2))

TopRegionsFull %<>% select(Symbol, RegionAnnot, CHR, RegionLocation, PeakLocation_PW, PeakLocation_NBB, LFC_PW, LFC_NBB, MetaP, Peak_PW, Peak_NBB, PropSignif)
TopRegionsFull$Symbol <- sapply(TopRegionsFull$Symbol, function(x){
  if(x %in% PDgenes$Nearest.Gene){
    paste0(x, "*")
  } else if(x == "PTPRH") {
    paste0(x, "**")
  } else {
    x
  }
})

TopRegionsFull %<>% mutate("LFC(PW,NBB)" = paste0(LFC_PW, ",", LFC_NBB),
                           "PeakName(PW, NBB)" = paste0(Peak_PW, ",", Peak_NBB),
                           "PeakLocation(PW, NBB)" = paste0(PeakLocation_PW, ",", PeakLocation_NBB))

write.table(TopRegionsFull, file = "GeneralResults/TopRegionsGWASBoth.tsv", row.names = F, col.names = T, sep = "\t")


#Identify genes which are significant with either RLE or house-keeping normalization
SingnifBothHKandRLE <- merge(CompareData$results %>% filter(CohorSignif == "Both") %>%
                               select(symbol, PeakName.NBB, PeakName.PV, log2FoldChange.NBB, log2FoldChange.PV, padj.NBB, padj.PV, UniqueRegionGeneType),
                             CompareData$resultsRLE %>% filter(CohorSignif == "Both") %>%
                               select(log2FoldChange.NBB, log2FoldChange.PV, padj.NBB, padj.PV, UniqueRegionGeneType),
                             by = "UniqueRegionGeneType", suffixes = c("_HK", "RLE"), all.x = F, ally = F) %>% arrange(log2FoldChange.NBB_HK)

SingnifBothHKandRLEInfo <- sapply(SingnifBothHKandRLE$UniqueRegionGeneType, function(DAR){
  Gene = rev(strsplit(DAR, "_")[[1]])[2]
  dataPV <- resultsPV %>% filter(symbol == Gene)
  dataPV$UniqueRegionGeneType <- sapply(dataPV$UniqueRegionGeneType, function(x){
    gsub("_hg19_genes", "", x)
  })
  dataPV %<>% filter(UniqueRegionGeneType ==  DAR)
  
  dataNBB <- resultsNBB %>% filter(symbol == Gene)
  dataNBB$UniqueRegionGeneType <- sapply(dataNBB$UniqueRegionGeneType, function(x){
    gsub("_hg19_genes", "", x)
  })
  dataNBB %<>% filter(UniqueRegionGeneType ==  DAR)
  
  PeakNameNBB = SingnifBothHKandRLE %>% filter(UniqueRegionGeneType == DAR) %>% .$PeakName.NBB
  PeakNamePV = SingnifBothHKandRLE %>% filter(UniqueRegionGeneType == DAR) %>% .$PeakName.PV

  PeakLocationNBB = dataNBB %>% filter(PeakName == PeakNameNBB, UniqueRegionGeneType == DAR) %>%
    .$Peak.Location 
  PeakLocationPW = dataPV %>% filter(PeakName == PeakNamePV, UniqueRegionGeneType == DAR) %>%
    .$Peak.Location
  data.frame(UniqueRegionGeneType = DAR,
             Chr = strsplit(PeakLocationPW, ":")[[1]][1] %>% gsub("chr", "", .),
             PeakLocationPW = strsplit(PeakLocationPW, ":")[[1]][2],
             PeakLocationNBB = strsplit(PeakLocationNBB, ":")[[1]][2],
             RegionAnnotation = rev(strsplit(DAR, "_")[[1]])[1])
}, simplify = F) %>% do.call(rbind, .)

SingnifBothHKandRLEFull <- merge(SingnifBothHKandRLE, SingnifBothHKandRLEInfo, by = "UniqueRegionGeneType", sort = FALSE)
SingnifBothHKandRLEFull$RegionLocation <- sapply(SingnifBothHKandRLEFull$UniqueRegionGeneType, function(x){
  x = gsub("chr[0-9]+:", "", x)
  strsplit(x, "_")[[1]][1]
})

SingnifBothHKandRLEFull %<>% select(-UniqueRegionGeneType)
names(SingnifBothHKandRLEFull) <- c("Symbol", "Peak_NBB", "Peak_PW", "LFC.HK_NBB", "LFC.HK_PW", "padj.HK_NBB", "padj.HK_PW",
                                    "LFC.RLE_NBB", "LFC.RLE_PW", "padj.RLE_NBB", "padj.RLE_PW", "CHR", "PeakLocation_PW", "PeakLocation_NBB", "RegionAnnot", "RegionLocation")
SingnifBothHKandRLEFull %<>% mutate_at(vars(matches("padj")), signif, digits = 1)
SingnifBothHKandRLEFull %<>% mutate_at(vars(matches("LFC")), round, digits = 2)

SingnifBothHKandRLEFull %<>% select(Symbol,CHR, RegionAnnot,  RegionLocation, PeakLocation_PW, PeakLocation_NBB, LFC.HK_PW, LFC.HK_NBB, LFC.RLE_PW, LFC.RLE_NBB,
                                    padj.HK_PW, padj.HK_NBB, padj.RLE_PW, padj.RLE_NBB, Peak_PW, Peak_NBB) %>% arrange(Symbol)

SingnifBothHKandRLEFull$Symbol <- sapply(SingnifBothHKandRLEFull$Symbol, function(x){
  if(x %in% PDgenes$Nearest.Gene){
    paste0(x, "*")
  } else if(x == "PTPRH") {
    paste0(x, "**")
  } else {
    x
  }
})

SingnifBothHKandRLEFull %<>% gather(Cohort, value = Info, matches("PW|NBB"))
SingnifBothHKandRLEFull$Var <- sapply(SingnifBothHKandRLEFull$Cohort, function(x){
  x = strsplit(x, "_")[[1]][1]
})

SingnifBothHKandRLEFull$Cohort <- sapply(SingnifBothHKandRLEFull$Cohort, function(x){
  x = strsplit(x, "_")[[1]][2]
})


SingnifBothHKandRLEFull %<>% spread(Var, value = Info) %>% mutate(RegionAnnot = paste0(RegionLocation, "_", RegionAnnot)) %>% select(-RegionLocation) %>% arrange(RegionAnnot, desc(Cohort))




write.table(SingnifBothHKandRLEFull, file = "GeneralResults/SingnifBothHKandRLEFull.tsv", row.names = F, col.names = T, sep = "\t")

######################### Count correlation ##############
RNAcounts <- readRDS("data/AdjustedCountsGon.Rds")
NBB_DESeq2 <-readRDS("ResultsNBB_Final/NBBDEoutput.Rds")
PV_DESeq2 <- readRDS("ResultsPV_Final/PVDEoutput.Rds")
PV_DESeq2RLE <- readRDS("ResultsPV_Final/PVDEoutputRLE.Rds")
NBB_DESeq2RLE <- readRDS("ResultsNBB_Final/NBBDEoutputRLE.Rds")

########## Getting the metadata #############
load("meta/Metadata.Rda")
Metadata %<>% mutate(subjectID = paste0("X", activemotif_id))

###################### Setting the adjustment covariates #######################
AdjCovarPV <- data.frame(Cov = c("sexM","batchC", "batchD",  "age","pm_hours",
                                 "Oligo_MSP"),
                         adjType = c(rep("base", 3), rep("mean", 3)))

AdjCovarNBB <- data.frame(Cov = c("sexM","batchB", "age","pm_hours",
                                  "Oligo_MSP"),
                          adjType = c(rep("base", 2), rep("mean", 3)))
###############################################################################
promoterPV <- data.table(resultsPV)[type == "hg19_genes_promoters",.(symbol, PeakName, Peak_Gene, UniqueRegion, type)]
RNApeaksPV <- promoterPV[symbol %in% RNAcounts$AdjustedPV$GeneSymbol,]


#Repeat for NBB
promoterNBB <- data.table(resultsNBB)[type == "hg19_genes_promoters",.(symbol, PeakName, Peak_Gene, UniqueRegion, type)]
RNApeaksNBB <- promoterNBB[symbol %in% RNAcounts$AdjustedNBB$GeneSymbol,]

GetChIP_RNAcor <- function(RNApeaks, DESseqOut, AdjCovar, Cohort = "PV", Name){
  RNAcohort = paste0("Adjusted", Cohort)
  GeneList <- as.list(RNApeaks$PeakName)
  names(GeneList) <- RNApeaks$PeakName
  AdjustedPromoterPeak <- lapply(GeneList, function(PeakName){
    GetAdjCountDESeq(dds = DESseqOut, Gene = PeakName, adjCov = AdjCovar) %>% t %>% data.frame
  }) %>% rbindlist()#, mc.cores = detectCores()) %>% rbindlist()
  
  names(AdjustedPromoterPeak) <- attr(DESseqOut, "colData")$SampleID
  
  AdjustedPromoterPeak$PeakName <- sapply(rownames(AdjustedPromoterPeak), function(x){
    strsplit(x, "\\.")[[1]][1]
  })
  
  AdjustedPromoterPeak$PeakGene <- RNApeaks$Peak_Gene
  AdjustedPromoterPeak$GeneSymbol  <- sapply(AdjustedPromoterPeak$PeakGene, function(x){
    strsplit(x, "_")[[1]][3]
  })
  
  #Remove duplicates and genes annotated to miRNAs/snRNAs
  AdjustedPromoterPeak %<>% filter(!duplicated(PeakGene)) %>% .[!.$GeneSymbol %in% symbolToFilter,]
  
  AdjustedPromoterPeak_Melt <- gather(AdjustedPromoterPeak, key = "SubjectID", value = AdjPromoter, -PeakName, -PeakGene, -GeneSymbol)
  AdjustedPromoterPeak_Melt$RNAid <- Metadata$sample_id_rna[match(AdjustedPromoterPeak_Melt$SubjectID, Metadata$subjectID)]
  
  CommonSamples <- intersect(as.character(AdjustedPromoterPeak_Melt$RNAid), names(RNAcounts[[RNAcohort]]))
  CombinedData <- sapply(CommonSamples, function(subj){
    SubjData <- RNAcounts[[RNAcohort]] %>% select(subj, GeneSymbol)
    temp <- merge(AdjustedPromoterPeak_Melt %>% filter(RNAid == subj), SubjData, by = "GeneSymbol", sort = F)
    names(temp)[ncol(temp)] <- "AdjExpr"
    temp
  }, simplify = F) %>% do.call(rbind, .) #%>% filter(is.finite(AdjPromoter), is.finite(AdjExpr))
  
  CombinedData$Group <- Metadata$condition[match(CombinedData$RNAid, Metadata$sample_id_rna)] %>%
    droplevels() %>% factor(levels = c("Cont", "PD"))
  
  ChIP_expCor <- sapply(unique(CombinedData$GeneSymbol), function(Gene){
    subData <-  CombinedData %>% filter(GeneSymbol == Gene)
    CorCont <- sapply(unique(subData$PeakGene), function(Peak){
      cor.test(~AdjPromoter+AdjExpr, data = subData %>% filter(Group == "Cont", PeakGene == Peak))$estimate %>% signif(digits = 2)
    }) %>% max
    
    CorPD <- sapply(unique(subData$PeakGene), function(Peak){
      cor.test(~AdjPromoter+AdjExpr, data = subData %>% filter(Group == "PD", PeakGene == Peak))$estimate %>% signif(digits = 2)
    }) %>% max
    data.frame(GeneSymbol = Gene, CorCont = CorCont, CorPD = CorPD )
  }, simplify = F) %>% do.call(rbind,.)
  
  #Look at the delta shift between the correlations
  ThreshData <- sapply(seq(0,1, 0.1), function(Thresh){
    subData = ChIP_expCor %>% filter(!(abs(CorCont) < Thresh & abs(CorPD) < Thresh))
    subData %<>% mutate(Threshold = Thresh, GeneNum = nrow(subData))
    gather(subData, key = "Group", value = "Cor", -GeneSymbol, -Threshold, -GeneNum)
  }, simplify = F) %>% do.call(rbind, .) 
  
  ThreshData$Group <- factor(ThreshData$Group)
  ThreshData$Threshold2 <- paste0("|Cor|>", ThreshData$Threshold, " , (",  ThreshData$GeneNum ,")")
  
  Plot <- ggplot(ThreshData, aes(Group, Cor, color = Group)) +
    theme_classic() +
    geom_boxplot(aes(fill = Group), alpha = 0.6) +
    labs(x = "", y = "Pearson's correlation") +
    scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
    scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
    geom_hline(yintercept = 0, col = "red", linetype = 2) +
    facet_wrap(~Threshold2, nrow = 2)
  
  ggsave(paste0(ResultsPath,"ChiP_ExprCorByGroup", Name, ".pdf"), Plot, device = "pdf", width = 8, height = 2.5, dpi = 300, useDingbats = F)
  
  WilCoxOut <- sapply(seq(0,0.9, 0.1), function(Thresh){
    subData = ChIP_expCor %>% select(-GeneSymbol) %>% filter(!(abs(CorCont) < Thresh & abs(CorPD) < Thresh))
    WilcoxDelta <- if(nrow(subData) > 0){
      MedianCont = subData$CorCont %>% median(na.rm = T)
      MedianPD = subData$CorPD %>% median(na.rm = T)
      round(wilcox.test(subData$CorCont, subData$CorPD, conf.int = T)$estimate, digits = 2)
    } else {
      WilcoxDelta = NA
      MedianCont = NA
      MedianPD = NA
    }
    data.frame(Threshold = Thresh, MedianCont = MedianCont,
               MedianPD = MedianPD, WilcoxDelta = WilcoxDelta)
  }, simplify = F) %>% do.call(rbind, .)
  
  
  Plot2 <- ggplot(gather(WilCoxOut,
                         key = "Group", value = "MedianCor", - Threshold, -WilcoxDelta),
                  aes(Threshold, MedianCor)) +
    theme_classic() +
    labs(y = "Pearson's correlation") +
    geom_col(data = WilCoxOut, aes(x = Threshold, y = WilcoxDelta), fill = "grey")  +
    geom_point(aes(color = Group)) +
    scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "") +
    geom_text(aes(x = Threshold, y = c(WilcoxDelta + 0.06) , label = WilcoxDelta),
              size = 2, color = "grey")
  ggsave(paste0(ResultsPath,"ChiP_ExprCorMedian", Name, ".pdf"), Plot2, device = "pdf", width = 4, height = 2.5, dpi = 300, useDingbats = F)
  
    return(list(AdjustedPromoterPeak_Melt = AdjustedPromoterPeak_Melt,
                CombinedData = CombinedData,
                ChIP_expCor = ChIP_expCor,
                ThreshData = ThreshData,
                BoxPlot = Plot,
                MedianPlot = Plot2,
                WilCoxOut = WilCoxOut))
}

ChIPrnaPV <- GetChIP_RNAcor(RNApeaks = RNApeaksPV, AdjCovar = AdjCovarPV,
                            DESseqOut = PV_DESeq2, Cohort = "PV", Name = "PV")



ChIPrnaPV$ChIP_expCor$FC <- resultsPV$log2FoldChange[match(ChIPrnaPV$ChIP_expCor$GeneSymbol, resultsPV$symbol)]

ChIPrnaPV$ChIP_expCor %<>% mutate(Delta = CorPD - CorCont)


ChIPrnaPV_RLE <- GetChIP_RNAcor(RNApeaks = RNApeaksPV, AdjCovar = AdjCovarPV,
                                DESseqOut = PV_DESeq2RLE, Cohort = "PV", Name = "PV_RLE")
ChIPrnaPV_RLE$ChIP_expCor %<>% mutate(Delta = CorPD - CorCont)


ChIPrnaNBB <- GetChIP_RNAcor(RNApeaks = RNApeaksNBB, AdjCovar = AdjCovarNBB,
                             DESseqOut = NBB_DESeq2, Cohort = "NBB", Name = "NBB")

ChIPrnaNBB$ChIP_expCor$FC <- resultsNBB$log2FoldChange[match(ChIPrnaNBB$ChIP_expCor$GeneSymbol, resultsNBB$symbol)]
ChIPrnaNBB$ChIP_expCor %<>% mutate(Delta = CorPD - CorCont)

MergedChIPrna <- merge(ChIPrnaPV$ChIP_expCor, ChIPrnaNBB$ChIP_expCor, by = "GeneSymbol", suffixes = c(".PV", ".NBB")) %>% select(GeneSymbol, CorCont.PV, CorCont.NBB, CorPD.PV, CorPD.NBB, Delta.PV, Delta.NBB, FC.PV, FC.NBB)
rownames(MergedChIPrna) <- MergedChIPrna$GeneSymbol
MergedChIPrna <- merge(MergedChIPrna, PDgenes, by.x = "GeneSymbol", by.y = "Gene", all.x = T, sort =F)
MergedChIPrna %<>% mutate(AdjmetaP = as.numeric(AdjMetaP))


ggplot(MergedChIPrna %>% filter(CorCont.PV > 0.5 & CorCont.NBB > 0.5), aes(Delta.PV, Delta.NBB)) +
  geom_point(alpha = 0.2) +
  geom_label_repel(data = MergedChIPrna %>% filter(CorCont.PV > 0.5 & CorCont.NBB > 0.5) %>% filter(Delta.PV < -0.8 & Delta.NBB < -0.8),
                   arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
                   nudge_y = 0.3, nudge_x = 0.3, aes(label = GeneSymbol), show.legend = F)


ChIPrnaNBB_RLE <- GetChIP_RNAcor(RNApeaks = RNApeaksNBB, AdjCovar = AdjCovarNBB,
                                DESseqOut = NBB_DESeq2RLE, Cohort = "NBB", Name = "NBB_RLE")
ChIPrnaNBB_RLE$ChIP_expCor %<>% mutate(Delta = CorPD - CorCont)

MergedChIPrnaRLE <- merge(ChIPrnaPV_RLE$ChIP_expCor, ChIPrnaNBB_RLE$ChIP_expCor, by = "GeneSymbol", suffixes = c(".PV", ".NBB")) %>% select(GeneSymbol, CorCont.PV, CorCont.NBB, CorPD.PV, CorPD.NBB, Delta.PV, Delta.NBB)
rownames(MergedChIPrnaRLE) <- MergedChIPrna$GeneSymbol



GetSpecificGeneCor <- function(Gene){
  chipPV <- ChIPrnaPV$AdjustedPromoterPeak_Melt %>% filter(GeneSymbol == Gene)
  rnaPV <- RNAcounts$AdjustedPV %>% filter(GeneSymbol == Gene) %>% select(-GeneSymbol) %>% t %>% data.frame()
  names(rnaPV) <- "adjCounts"
  
  chipPV$adjCounts <-rnaPV$adjCounts[match(chipPV$RNAid,rownames(rnaPV))]   

  chipNBB <- ChIPrnaNBB$AdjustedPromoterPeak_Melt %>% filter(GeneSymbol == Gene)
  rnaNBB <- RNAcounts$AdjustedNBB %>% filter(GeneSymbol == Gene) %>% select(-GeneSymbol) %>% t %>% data.frame()
  names(rnaNBB) <- "adjCounts"
  
  chipNBB$adjCounts <-rnaNBB$adjCounts[match(chipNBB$RNAid,rownames(rnaNBB))]   
  
  combinedData <- rbind(chipPV, chipNBB)
  combinedData$Group <- Metadata$condition[match(combinedData$SubjectID,Metadata$subjectID)] %>% factor(levels = c("Cont", "PD"))
  combinedData$Cohort <- Metadata$cohort[match(combinedData$SubjectID,Metadata$subjectID)]  %>% factor(levels = c("PV", "NBB"))
  levels(combinedData$Cohort) <- c("PW", "NBB")
  combinedData$PeakName <- factor(combinedData$PeakName)

  Stat <- sapply(levels(combinedData$Cohort), function(cohortName){
    CohortData <- combinedData %>% filter(Cohort == cohortName) %>% droplevels
    sapply(levels(CohortData$PeakName), function(peakName){
      PeakData <- CohortData %>% filter(PeakName == peakName) %>% droplevels
      sapply(levels(PeakData$Group), function(groupName){
        groupData = PeakData %>% filter(Group == groupName) %>% droplevels
        cor.test(~adjCounts+AdjPromoter, data = groupData, method = "pearson")$est
      }, simplify = F) %>% do.call(cbind, .) %>% data.frame %>% mutate(PeakName = peakName, Cohort = cohortName)
    }, simplify = F) %>% do.call(rbind, .) %>% data.frame
  }, simplify = F)
  
  MaxConPW <- Stat$PW %>% arrange(desc(Cont)) %>% .[1,] 
  MaxPDPW <- Stat$PW %>% arrange(desc(PD)) %>% .[1,] 
  MaxConNBB <- Stat$NBB %>% arrange(desc(Cont)) %>% .[1,] 
  MaxPDNBB <- Stat$NBB %>% arrange(desc(PD)) %>% .[1,]
  
  DataToPlot <- rbind(combinedData %>% filter(Cohort == "PW", PeakName == MaxConPW$PeakName),
                      combinedData %>% filter(Cohort == "NBB", PeakName == MaxConNBB$PeakName)) %>%
    data.frame() %>% filter(!is.na(adjCounts)) %>% droplevels()
  
  StatPlot <- rbind(MaxConPW, MaxConNBB) %>% data.frame() %>%
    gather(key = "Group", value = "Cor", -PeakName, -Cohort) %>%
    mutate(x = rep(c(4.5,8.3), 2),
           y = c(12.1, 12.28, 12.05, 12.25),
           Text = paste0("Cor(", Group, ") = ", round(Cor, digits = 2))) 
  StatPlot$Cohort <- factor(StatPlot$Cohort, levels = c("PW", "NBB"))
  
  levels(DataToPlot$Cohort) <- c(paste0("PW Cor(Ctrl) = ",  round(MaxConPW$Cont, digits = 2), " Cor(PD) = ",  round(MaxConPW$PD, digits = 2)),
                                 paste0("NBB Cor(Ctrl) = ",  round(MaxConNBB$Cont, digits = 2), " Cor(PD) = ",  round(MaxConNBB$PD, digits = 2)))
  
  ggplot(DataToPlot, aes(AdjPromoter, adjCounts, color = Group)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = paste0("Adjusted log2(", Gene, " ChIPseq promoter counts)"), y = paste0("Adjusted log2(", Gene, " RNAseq counts)")) +
    geom_smooth(method = "lm", formula = y~x, aes(fill = Group),  alpha = 0.2) +
    geom_point() +
    scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
    scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
    #geom_text(data = StatPlot, aes(x, y, label = Text, color = Group), inherit.aes = F, show.legend = F) +
    facet_wrap(~Cohort, scales = "free")

}

Plot <- GetSpecificGeneCor("MED13")
ggsave(paste0(ResultsPath,"CorMED13.pdf"), 
       device = "pdf", Plot, width = 5, height = 2.5, dpi = 300, useDingbats = F)

Plot <- GetSpecificGeneCor("PTPRH")
ggsave(paste0(ResultsPath,"CorPTPRH.pdf"), 
       device = "pdf", Plot, width = 5, height = 2.5, dpi = 300, useDingbats = F)

Plot <- GetSpecificGeneCor("DLG2")
ggsave(paste0(ResultsPath,"CorDLG2.pdf"), 
       device = "pdf", Plot, width = 5, height = 2.5, dpi = 300, useDingbats = F)

save.image(file = "GeneralResults/PVNBBcomparison.RData")
save(CommonRegions, MergedChIPrna, file = paste0(ResultsPath, "CombinedData.Rda"))
