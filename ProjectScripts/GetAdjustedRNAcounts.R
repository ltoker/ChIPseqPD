GonDESeqOut <- readRDS("data/DESeqOut_CT.Rds")

GetDESeq2Results <- function(DESeqOut, coef, alpha = 0.05, indepFilter = TRUE){
  DEresults <- results(DESeqOut, name = coef, alpha = alpha, format = "DataFrame", independentFiltering = indepFilter)
  DEresults$GeneSymbol <- geneNames$hgnc_symbol[match(rownames(DEresults), geneNames$ensembl_gene_id)]
  DEresults$EnsemblID <- rownames(DEresults)
  DEresults %<>% data.frame %>% filter(GeneSymbol != "")
  return(DEresults)
}

GetAdjCountDESeq <- function(dds, Gene,  adjCov){
  GeneRow <- which(rownames(dds) == Gene)
  Mu <- log2(t(t(assays(dds)$mu[GeneRow,])/sizeFactors(dds)))
  Counts <-  log2(t(t(assays(dds)$counts[GeneRow,])/sizeFactors(dds)))
  Resid <- Counts - Mu
  corrFactor <- sizeFactors(dds)
  Coef <- coef(dds)[GeneRow,]
  mod <- attr(dds, "modelMatrix")
  modAdj <-mod
  for(cov in as.character(adjCov$Cov)){
    adjType = adjCov %>% filter(Cov == cov) %>% .$adjType %>% as.character()
    if(adjType == "mean"){
      modAdj[,cov] <- mean(modAdj[,cov], na.rm=T)
      
    } else if (adjType == "base"){
      modAdj[,cov] <- 0
    }
  }
  AdjValue <- (modAdj %*% Coef) + Resid
  return(AdjValue)
}




packageF("EnsDb.Hsapiens.v75")
packageF("ensembldb")
edb <- EnsDb.Hsapiens.v75
txdb <- as.data.frame(transcriptsBy(edb, by="gene",
                                    columns=c("gene_id", "tx_id", "gene_name", "gene_biotype")))

# RESTRICT TO NUCLEAR GENES IN CANONICAL CHROMOSOMES
txdb <- subset(txdb, seqnames %in% c(1:22,"X","Y", "MT") & startsWith(gene_id, 'ENSG'))
txdb$seqnames <- as.character(txdb$seqnames)

geneNames <- txdb %>% select(gene_name, gene_id, gene_biotype)
names(geneNames) <- c("hgnc_symbol", "ensembl_gene_id", "gene_biotype")


GonDESeqResultsDF <- lapply(GonDESeqOut, function(cohort){
  GetDESeq2Results(cohort)
})


AdjCovarPV <- data.frame(Cov = c("sex_M_vs_F", 
                                 "Batch_2_vs_1", "Batch_3_vs_1", "Batch_4_vs_1",
                                 "age_years","PMI_hours", "rin",
                                 "Oligo_Genes"),
                         adjType = c(rep("base", 4), rep("mean", 4)))

AdjCovarNBB <- data.frame(Cov = c("sex_M_vs_F", "Batch_1_vs_0",
                                  "Batch_2_vs_0", "Batch_3_vs_0", "Batch_4_vs_0",
                                  "age_years","PMI_hours", "rin",
                                  "Oligo_Genes", "Microglia_Genes"),
                          adjType = c(rep("base", 5), rep("mean", 5)))



#Get the corrected counts for all genes
AdjustedPV <- as.list(GonDESeqResultsDF$PW$GeneSymbol)
names(AdjustedPV) <- GonDESeqResultsDF$PW$GeneSymbol

AdjustedPV <- mclapply(AdjustedPV, function(Gene){
  EnsemblID = GonDESeqResultsDF$PW %>% filter(GeneSymbol == Gene) %>% .$EnsemblID
  temp <- GetAdjCountDESeq(dds = GonDESeqOut$PW, Gene = EnsemblID, adjCov = AdjCovarPV) %>% t %>% data.frame()
}, mc.cores = detectCores()) %>% rbindlist()

AdjustedPV %<>% mutate(GeneSymbol = as.character(GonDESeqResultsDF$PW$GeneSymbol))


AdjustedNBB <- as.list(GonDESeqResultsDF$NBB$GeneSymbol)
names(AdjustedNBB) <- GonDESeqResultsDF$NBB$GeneSymbol

AdjustedNBB <- mclapply(AdjustedNBB, function(Gene){
  EnsemblID = GonDESeqResultsDF$NBB %>% filter(GeneSymbol == Gene) %>% .$EnsemblID
  temp <- GetAdjCountDESeq(dds = GonDESeqOut$NBB, Gene = EnsemblID, adjCov = AdjCovarNBB) %>% t %>% data.frame()
}, mc.cores = detectCores()) %>% rbindlist()

AdjustedNBB %<>% mutate(GeneSymbol = as.character(GonDESeqResultsDF$NBB$GeneSymbol))


saveRDS(list(AdjustedPV = AdjustedPV, AdjustedNBB = AdjustedNBB), "data/AdjustedCountsGon.Rds")
