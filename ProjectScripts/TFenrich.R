source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
packageF("ggpubr")

ResultsPath = "GeneralResults/"


DAresults <- list(NBB = readRDS("ResultsNBB_Final/NBBDEresults.Rds") %>%
                    select(PeakName, log2FoldChange, pvalue, padj, type, ),
                  PV = readRDS("ResultsPV_Final/PVDEresults.Rds") %>%
                    select(PeakName, log2FoldChange, pvalue, padj, type))


TF_PeakData <- sapply(names(DAresults), function(Cohort){
  if(Cohort == "NBB"){
    narrowPeak <- read.table(paste0("data/Peaks/", Cohort, "_ALL.narrowPeakClean.gz"), header = F, sep = "\t")
  } else {
    narrowPeak <- read.table(paste0("data/Peaks/", Cohort, "_ALL_RNA.narrowPeakClean.gz"), header = F, sep = "\t")
  }
  names(narrowPeak)[1:5] <- c("CHR", "START", "END", "PeakName", "Score")
  narrowPeak <- narrowPeak[!grepl("GL|hs",narrowPeak$CHR),]
  narrowPeak %<>%  mutate(CHR = paste0("chr", CHR))
  
  Data <- merge(DAresults[[Cohort]] %>% filter(!duplicated(PeakName)), narrowPeak, by = "PeakName", all.x = T, all.y = T, sort = F)
  Data %<>% mutate(Width = END - START) %>%
    select(-Score, -matches("^V.+$"))
  
  Data %<>% filter(!is.na(padj))
  
  Data$Signif <- "NS"
  Data$Signif[Data$padj < 0.05 & Data$log2FoldChange < 0] <- "Hypo"
  Data$Signif[Data$padj < 0.05 & Data$log2FoldChange > 0] <- "Hyper"
  
  BackgroundRand <-   sapply(c("Hyper", "Hypo"), function(Type){
    DataSignif <- Data %>% filter(Signif == Type)
    SignifTot = nrow(DataSignif %>% filter(!duplicated(PeakName)))
    SignifWidth <- quantile(DataSignif$Width, seq(0,1, 0.1))
    
    MatchedNSPeaks <- list()
    for(i in 1:10){
      min = SignifWidth[i]
      max = SignifWidth[i+1]
      temp <- Data %>% filter(Signif != Type, Width >= min, Width <= max, !duplicated(PeakName))
      MatchedNSPeaks[[paste0("Width_", min, "_", max)]] <- temp$PeakName
    }
    
    lapply(MatchedNSPeaks, function(WidthGroup){
      sample(WidthGroup, 0.5*SignifTot, replace = F)
    }) %>% unlist
  }, simplify = F)
  

  list(Data = Data,
       BackgroundRand = BackgroundRand)
}, simplify = F)


#Background for DARs: peak-length matched. Bed files are used to retrieve the peak sequences using http://rsat.sb-roscoff.fr/
GetBed <- function(Cohort, Type){
  CohortData = TF_PeakData[[Cohort]]
  
  DARs <- CohortData$Data  %>% filter(Signif == Type) %>% 
    select(CHR, START,END, PeakName)
  
  WidthMatched <-  CohortData$Data %>%
    filter(PeakName %in% CohortData$BackgroundRand[[Type]]) %>%
    select(CHR, START,END, PeakName)
  

  write.table(DARs, paste0(ResultsPath, Type, "Peak", Cohort, ".bed"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  write.table(WidthMatched, paste0(ResultsPath, "Non", Type, "Peak", Cohort, ".bed"),
              sep = "\t", row.names = F, col.names = F, quote = F)
}


GetBed("NBB", "Hyper")
GetBed("NBB", "Hypo")
GetBed("PV", "Hyper")
GetBed("PV", "Hypo")



#Fix the Peak naming because of the MEME-suite limitations
FixID <- function(Cohort, Type){
  FileIn = paste0("data/", Type, Cohort, "fasta.txt")
  FileOut = paste0("data/", Type, Cohort, "fastaFixed.txt")
  IDform = paste0(">", Type, "Peak_")
  
  temp <- readLines(FileIn)
  temp[grep(">", temp)] <- paste0(IDform, 1:(length(temp)/2))
  write.table(temp, FileOut, quote = F, row.names = F, col.names = F)
}

FixID("NBB", "NonHyper")
FixID("NBB", "NonHypo")
FixID("NBB", "Hyper")
FixID("NBB", "Hypo")
FixID("PV", "NonHyper")
FixID("PV", "NonHypo")
FixID("PV", "Hyper")
FixID("PV", "Hypo")


GetBaseFreq <- function(SeqObject){
  SeqOnly <- SeqObject[!grepl(">", SeqObject)]
  
  SeqFreqPeak <- sapply(SeqOnly, function(x){
    Length = stringr::str_count(x, ".")
    data.frame(A = stringr::str_count(x, "A")/Length,
               "T" = stringr::str_count(x, "T")/Length,
               C = stringr::str_count(x, "C")/Length,
               G= stringr::str_count(x, "G")/Length)
  }, simplify = F) %>% rbindlist()
  
  
  SummaryFreq <- sapply(names(SeqFreqPeak), function(x){
    temp <- summary(SeqFreqPeak[[x]])
    DF <- temp %>% array %>% data.frame() %>% t
    names(DF) <- attr(temp, "names")
    DF$Base <- x
    DF 
  }, simplify = F) %>% rbindlist %>%
    data.frame() %>% select(c(ncol(.), 1:ncol(.)-1))
  
  return(list(PerPeak = SeqFreqPeak,
              SummFreq = SummaryFreq))
}


AllSeq <- sapply(list.files(path = "data/", pattern = "Fixed"), function(Type){
  readLines(paste0("data/", Type))
},  simplify = F)


AllBaseFreq <- lapply(AllSeq, function(x){
  GetBaseFreq(x)
}) 

AllBaseFreqDF <- sapply(names(AllBaseFreq), function(x){
  Cohort = ifelse(grepl("NBB", x), "NBB", "PW")
  
  Type = if(grepl("^Hyper", x)){
    "Hyper"
  } else if (grepl("^Hypo", x)){
    "Hypo"
  } else {
    "NotDAR"
  }
  
  Data = AllBaseFreq[[x]]$PerPeak
  Data %>% mutate(GC = 100*(G+C)/(G+C+A+T), Cohort = Cohort, Type = Type) %>% data.frame() 
}, simplify = F) %>% rbindlist()

AllBaseFreqDF$Cohort <- factor(AllBaseFreqDF$Cohort, levels = c("PW", "NBB"))

GCplot <- ggplot(AllBaseFreqDF) + geom_density(aes(GC, fill = Type, color = Type), alpha = 0.3) +
  theme_minimal() +
  labs(x = "GC%") +
  facet_wrap(~Cohort)

ggsave("GCcontent.pdf", plot = GCplot, device = "pdf", path = ResultsPath, width = 8, height = 5, useDingbats = F)



#Repeat for exons and introns separately

lapply(names(TF_PeakData), function(Cohort){
  
  CohortData1 = TF_PeakData[[Cohort]]$Data %>% filter(!is.na(padj))
  CohortData2 = DAresults[[Cohort]] %>% filter(!is.na(padj))
  
  sapply(c("Hypo", "Hyper"), function(SignifType){
    sapply(c("exons", "introns"), function(RegionType){
      RegionType = paste0("hg19_genes_", RegionType)
      SignifName = paste0(RegionType, "_", SignifType, Cohort)
      BackgroundName = paste0("Non", RegionType, "_", SignifType, Cohort)
      
      SignifPeaksAll = CohortData1 %>% filter(Signif == SignifType) %>% .$PeakName
      
      
      SignifPeaks = CohortData2 %>% filter(PeakName %in% SignifPeaksAll, type == RegionType) %>%
        .$PeakName %>% unique
      
      SignifPeakData <- CohortData1 %>% filter(PeakName %in% SignifPeaks)
      SignifTot = length(SignifPeaks)
      
      SignifWidth <- quantile(SignifPeakData$Width, seq(0,1, 0.1))
      
      MatchedNSPeaks <- list()
      for(i in 1:10){
        min = SignifWidth[i]
        max = SignifWidth[i+1]
        temp <- CohortData2 %>% filter(!PeakName %in% SignifPeaks, type == RegionType) %>%
          filter(!duplicated(PeakName))
        temp2 <-  CohortData1 %>% filter(Width >= min, Width <= max)
        MatchedNSPeaks[[paste0("Width_", min, "_", max)]] <- temp$PeakName
      }
      
      BackgroundPeaks <- lapply(MatchedNSPeaks, function(WidthGroup){
        sample(WidthGroup, 0.5*SignifTot, replace = F)
      }) %>% unlist
      
      
      BackgroundData <- CohortData1 %>% filter(PeakName %in% BackgroundPeaks, !duplicated(PeakName))
      

      SignifPeakData %>% select(Peak.CHR, Peak.START,
                                Peak.END, PeakName) %>% write.table(paste0(ResultsPath, SignifName, ".bed"),
                                                                    sep = "\t", row.names = F, col.names = F, quote = F)
      BackgroundData %>% select(Peak.CHR, Peak.START,
                                Peak.END, PeakName) %>% write.table(paste0(ResultsPath, BackgroundName, ".bed"),
                                                                    sep = "\t", row.names = F, col.names = F, quote = F)
             
    })
  })
})

DAresults$PV %>% filter(!is.na(padj), type == "hg19_genes_exons") %>% 
  filter(!duplicated(PeakName)) %>%
  select(Peak.CHR, Peak.START,
         Peak.END, PeakName) %>% write.table("ResultsPV_Final/ExonAnnoAllPeak.bed",
                                             sep = "\t", row.names = F, col.names = F, quote = F)

DAresults$PV %>% filter(!is.na(padj), type == "hg19_genes_introns") %>% 
  filter(!duplicated(PeakName)) %>%
  select(Peak.CHR, Peak.START,
         Peak.END, PeakName) %>% write.table("ResultsPV_Final/IntronAnnoAllPeak.bed",
                                             sep = "\t", row.names = F, col.names = F, quote = F)
DAresults$PV %>% filter(padj < 0.05, type == "hg19_genes_exons") %>% 
  filter(!duplicated(PeakName)) %>%
  select(Peak.CHR, Peak.START,
         Peak.END, PeakName) %>% write.table("ResultsPV_Final/ExonAnnoSignifPeak.bed",
                                             sep = "\t", row.names = F, col.names = F, quote = F)

DAresults$PV %>% filter(padj < 0.05, type == "hg19_genes_introns") %>% 
  filter(!duplicated(PeakName)) %>%
  select(Peak.CHR, Peak.START,
         Peak.END, PeakName) %>% write.table("ResultsPV_Final/IntronAnnoSignifPeak.bed",
                                             sep = "\t", row.names = F, col.names = F, quote = F)

All_Seq <- list(Seq_ExonsAll = readLines("ResultsPV_Final/All_Exon.txt"),
                Seq_IntonAll = readLines("ResultsPV_Final/All_Intron.txt"),
                Seq_Introns_Signif = readLines("ResultsPV_Final/Signif_Intron.txt"),
                Seq_Exons_Signif = readLines("ResultsPV_Final/Signif_Exon.txt"))



AllBaseFreq <- lapply(All_Seq, function(x){
  GetBaseFreq(x)
}) 


SummaryFreqDF <- sapply(names(AllBaseFreq), function(x){
  AllBaseFreq[[x]]$SummFreq %>% data.frame %>%
    select(Base, Median) %>% mutate(Median = signif(Median, digits = 3),
                                    Type = x)
}, simplify = F) %>% rbindlist() %>% data.frame %>% 
  select(Type, Base, Median)


SummaryFreqDF$regionType <- sapply(SummaryFreqDF$Type, function(x){
  if(grepl("Exon", x)){
    "Exons"
  } else {
    "Introns"
  }
})

SummaryFreqDF$PeakType <- sapply(SummaryFreqDF$Type, function(x){
  if(grepl("All", x)){
    "AllPeaks"
  } else {
    "SignifPeaks"
  }
})

SummaryFreqDF %>% arrange(PeakType, regionType) %>%
  select(regionType, PeakType, Base, Median)
