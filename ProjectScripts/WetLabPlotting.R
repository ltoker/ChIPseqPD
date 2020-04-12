ResultsPath = "GeneralResults/"
source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")

StatData <- read.table("data/Supplementary_Table_2.txt", sep = "\t", header = T, comment.char = "!", quote = "'")

Data1 <- read.table("data/PFC_acetylation_data_RESCALED_0-1.txt", header = T, sep = "\t", comment.char = "")
Data2 <- read.table("data/PanAcetyl.txt", header = T, sep = "\t", comment.char = "") %>% filter(!is.na(Status)) 
Data3 <- read.table("data/Striatum_Cerebellum_acetylation_mixed_models_rescaled_0-1_ALLFINAL.txt", header = T, sep = "\t", comment.char = "") %>% filter(!is.na(Status))

Data <- merge(Data1 %>% select(-Status), Data2 %>% select(matches("ID$|Acet|Status")), by = "ID", all.y = T)

Data$Group <- sapply(Data$Status, function(x){
  if(x == 0){
    "Cont"
  } else {
    "PD"
  }
}) %>% factor


DataAcetylation <- Data %>% gather(key = "Measure", value = "Value", matches("ac$|Acet", ignore.case = F)) %>% filter(!is.na(Value))

DataAcetylationSum <- DataAcetylation %>% group_by(Measure, ID) %>% summarise(n = n(), Mean = mean(Value), SD = sd(Value)) %>% data.frame()

DataAcetylationSum$Measure <- sapply(DataAcetylationSum$Measure, function(x){
  gsub("Acetyl_lysin_17KDa", "Lys-Ac(17kDa)", x)
})

DataAcetylationSum <- merge(DataAcetylationSum, DataAcetylation %>% filter(!duplicated(ID)) %>%
                              select(ID, Group, Age, Gender,
                                     PM), by = "ID", sort = F, all.y = F)

#This is to fix the cases with only one replicate - for Lys-Ac multiple entrences were created because of the merging
DataAcetylationSum[DataAcetylationSum$Measure == "Lys-Ac(17kDa)",]$n <- 1
DataAcetylationSum[is.na(DataAcetylationSum$SD),]$SD <- 0
DataAcetylationSum$Measure <- factor(DataAcetylationSum$Measure, levels = c("Lys-Ac(17kDa)", "H2AK5ac", "H2BK15ac",
                                                                            "H3K9K14ac", "H3K27ac","H3K56ac",
                                                                            "H4K5ac", "H4K12ac", "H4K16ac"))
levels(DataAcetylationSum$Measure) <- c("Lys-Ac(17kDa):GAPDH", "H2AK5ac:H2A", "H2BK15ac:H2B",
                                        "H3K9/K14ac:H3", "H3K27ac:H3","H3K56ac:H3",
                                        "H4K5ac:H4", "H4K12ac:H4", "H4K16ac:H4")
StatAcetylCortex <- StatData %>% filter(Region == "Prefrontal cortex") %>% select(Marker, p.value)

AcetylStat <- data.frame(Group = 1.5,Mean = 1.1,
                         label = paste0 ("p = ", StatAcetylCortex$p.value[match(levels(DataAcetylationSum$Measure),
                                                              StatAcetylCortex$Marker)]),
                         Measure = factor(levels(DataAcetylationSum$Measure), levels = levels(DataAcetylationSum$Measure)))

Plot <- ggplot(DataAcetylationSum %>% filter(Measure != "Lys-Ac(17kDa):GAPDH"), aes(Group, Mean)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 4) +
  geom_text(data = AcetylStat %>% filter(Measure != "Lys-Ac(17kDa):GAPDH") , aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WBacetylResults.pdf"), Plot, device = "pdf", width = 4.5, height = 8, dpi = 300, useDingbats = F)
#ggsave(paste0(ResultsPath,"WBacetylResults.pdf"), Plot, device = "pdf", width = 16, height = 3, dpi = 300, useDingbats = F)

Plot <- ggplot(DataAcetylationSum %>% filter(Measure == "Lys-Ac(17kDa):GAPDH"), aes(Group, Mean)) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 1) +
  geom_text(data = AcetylStat %>% filter(Measure == "Lys-Ac(17kDa):GAPDH") , aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WBPanAcetylResults.pdf"), Plot, device = "pdf", width = 3.5, height = 3, dpi = 300, useDingbats = F)



#Histones
DataHistones <- Data %>% gather(key = "Measure", value = "Value", matches("^H..?$")) %>% filter(!is.na(Value))

DataHistoneSum <- DataHistones %>% group_by(Measure, ID) %>% summarise(n = n(), Mean = mean(Value), SD = sd(Value)) %>% data.frame()

DataHistoneSum <- merge(DataHistoneSum, DataHistones %>% filter(!duplicated(ID)) %>%
                          select(ID, Group, Age, Gender,
                                 PM), by = "ID", sort = F, all.y = F)
                              

DataHistoneSum$Measure <- factor(DataHistoneSum$Measure, unique(DataHistoneSum$Measure))
levels(DataHistoneSum$Measure) <- c("H2A:GAPDH", "H2B:GAPDH", "H3:GAPDH","H4:GAPDH")

StatHistoneCortex <- StatData %>% filter(Marker %in% levels(DataHistoneSum$Measure),Region == "Prefrontal cortex") %>% select(Marker, p.value)

HistoneStat <- data.frame(Group = 1.5,Mean = 1.1,
                         label = paste0 ("p = ", StatHistoneCortex$p.value[match(levels(DataHistoneSum$Measure),
                                                                                 StatHistoneCortex$Marker)]),
                         Measure = factor(levels(DataHistoneSum$Measure), levels = levels(DataHistoneSum$Measure)))

Plot <- ggplot(DataHistoneSum, aes(Group, Mean)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 1) +
  geom_text(data = HistoneStat, aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WB_HistoneResults.pdf"), Plot, device = "pdf", width = 10, height = 3.3, dpi = 300, useDingbats = F)

#a-tubulin
DataAlphaTub <- read.table("data/a-tubulin_NEW.txt", header = T, sep = "\t") %>%
  gather(key = "Measure", value = "Value", matches("A_")) %>% filter(!is.na(Value))

DataAlphaTub$Group <- sapply(DataAlphaTub$Status, function(x){
  if(x == 0){
    "Cont"
  } else {
    "PD"
  }
}) %>% factor

DataAlphaTubSum <- DataAlphaTub %>% group_by(Measure, ID) %>% summarise(n = n(), Mean = mean(Value), SD = sd(Value)) %>% data.frame()

DataAlphaTubSum <- merge(DataAlphaTubSum, DataAlphaTub %>% filter(!duplicated(ID)) %>%
                           select(ID, Group, Age, Gender,
                                  PM), by = "ID", sort = F, all.y = F)
                          


DataAlphaTubSum$Measure <- factor(DataAlphaTubSum$Measure, unique(DataAlphaTubSum$Measure))
levels(DataAlphaTubSum$Measure) <- c("a-tubulinK40c:a-tubulin")

StatAlphaTub <- StatData %>% filter(Marker == "a-tubulinK40c:a-tubulin") %>% select(Marker, p.value)

AlphaTubStat <- data.frame(Group = 1.5,Mean = 0.9,
                           label = paste0 ("p = ", StatAlphaTub$p.value[match(levels(DataAlphaTubSum$Measure),
                                                                              StatAlphaTub$Marker)]),
                           Measure = factor(levels(DataAlphaTubSum$Measure), levels = levels(DataAlphaTubSum$Measure)))
                                         
                         

Plot <- ggplot(DataAlphaTubSum, aes(Group, Mean)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 1) +
  geom_text(data = AlphaTubStat, aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WB_AlphaTub.pdf"), Plot, device = "pdf", width = 4.2, height = 3.3, dpi = 300, useDingbats = F)

#Sirtuins
DataSirtuin <- Data %>% gather(key = "Measure", value = "Value", matches("Sirt")) %>% filter(!is.na(Value))

DataSirtuinSum <- DataSirtuin %>% group_by(Measure, ID) %>% summarise(n = n(), Mean = mean(Value), SD = sd(Value)) %>% data.frame()

DataSirtuinSum <- merge(DataSirtuinSum, DataAlphaTub %>% filter(!duplicated(ID)) %>%
                          select(ID, Group, Age, Gender,
                                 PM), by = "ID", sort = F, all.y = F)
                           



DataSirtuinSum$Measure <- factor(DataSirtuinSum$Measure, unique(DataSirtuinSum$Measure))
levels(DataSirtuinSum$Measure) <- c("SIRT1:b-tubulin", "SIRT2:b-tubulin", "SIRT3:b-tubulin")

StatSirtuin <- StatData %>% filter(Marker %in% levels(DataSirtuinSum$Measure)) %>% select(Marker, p.value)

SirtStat <- data.frame(Group = 1.5,Mean = 1.1,
                       label = paste0 ("p = ", StatSirtuin$p.value[match(levels(DataSirtuinSum$Measure),
                                                                         StatSirtuin$Marker)]),
                       Measure = factor(levels(DataSirtuinSum$Measure), levels = levels(DataSirtuinSum$Measure)))
                           
Plot <- ggplot(DataSirtuinSum, aes(Group, Mean)) +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 1) +
  geom_text(data = SirtStat, aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WB_Sirtuin.pdf"), Plot, device = "pdf", width = 8, height = 3, dpi = 300, useDingbats = F)

#Cell markers
DataMarkers <- Data %>% gather(key = "Measure", value = "Value", NUP188, NeuN, GFAP, CNP1, CX3CR1) %>% filter(!is.na(Value))

DataMarkersSum <- DataMarkers %>% group_by(Measure, ID) %>% summarise(n = n(), Mean = mean(Value), SD = sd(Value)) %>% data.frame()

DataMarkersSum <- merge(DataMarkersSum, DataMarkers %>% filter(!duplicated(ID)) %>%
                          select(ID, Group, Age, Gender,
                                 PM), by = "ID", sort = F, all.y = F)


DataMarkersSum$Measure <- factor(DataMarkersSum$Measure, levels = c("NUP188", "NeuN", "GFAP", "CNP1", "CX3CR1"))
levels(DataMarkersSum$Measure) <- c("NUP188:GAPDH", "NeuN:GAPDH", "GFAP:GAPDH", "CNP1:GAPDH", "CX3CR1:GAPDH")

StatMarker <- StatData %>% filter(Marker %in% levels(DataMarkersSum$Measure), Region == "Prefrontal cortex") %>% select(Marker, p.value)

MarkerStat <- data.frame(Group = 1.5,Mean = 1.1,
                         label = paste0 ("p = ", StatMarker$p.value),
                         Measure = factor(StatMarker$Marker, levels = StatMarker$Marker))

Plot <- ggplot(DataMarkersSum, aes(Group, Mean)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 1) +
  geom_text(data = MarkerStat, aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WBMarkersResults.pdf"), Plot, device = "pdf", width = 10, height = 3, dpi = 300, useDingbats = F)                              


#Repeat for other brain regions
Data3$Group <- sapply(Data3$Status, function(x){
  if(x == 0){
    "Cont"
  } else {
    "PD"
  }
}) %>% factor


DataRegionsAcetylation <- Data3 %>% gather(key = "Measure", value = "Value", matches("H3")) %>% filter(!is.na(Value))
DataRegionsAcetylation$Region <- sapply(DataRegionsAcetylation$Measure, function(x){
  if(grepl("Striat", x)){
    "Striatum"
  } else {
    "Cerebellum"
  }
})

DataRegionsAcetylation$MeasureType <- sapply(DataRegionsAcetylation$Measure, function(x){
  if(grepl("_to_", x)){
    "Acetylation"
  } else {
    "Histone"
  }
})

DataRegionsAcetylationSum <- DataRegionsAcetylation %>% group_by(Region, MeasureType, Measure, ID1) %>% summarise(n = n(), Mean = mean(Value), SD = sd(Value)) %>% data.frame()

DataRegionsAcetylationSum$Measure <- sapply(DataRegionsAcetylationSum$Measure, function(x){
  strsplit(x, "_")[[1]][1]
})

DataRegionsAcetylationSum <- merge(DataRegionsAcetylationSum, DataRegionsAcetylation %>% filter(!duplicated(ID1)) %>%
                                     select(ID1, Group, Age, Gender,
                                            PM), by = "ID1", sort = F, all.y = F)
                              


DataRegionsAcetylationSum$Measure <- factor(DataRegionsAcetylationSum$Measure, levels = c("H3K9K14ac", "H3K27ac","H3"))

levels(DataRegionsAcetylationSum$Measure) <- c("H3K9/K14ac:H3", "H3K27ac:H3","H3:GAPDH")

StatAcetylRegions <- StatData %>% filter(Region != "Prefrontal cortex", Marker %in% levels(DataRegionsAcetylationSum$Measure)) %>% select(Marker, p.value, Region)

AcetylRegionStat <- StatAcetylRegions %>% filter(!is.na(p.value)) %>%
                                                   mutate(Group = 1.5,Mean = 1.1,
                                                          label = paste0 ("p = ", p.value),
                                                          Measure = Marker)
                              
                         

Plot <- ggplot(DataRegionsAcetylationSum %>% filter(Region == "Striatum"), aes(Group, Mean)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 1) +
  geom_text(data = AcetylRegionStat %>% filter(Region == "Striatum") , aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WBacetylResultsStriatum.pdf"), Plot, device = "pdf", width = 8, height = 3, dpi = 300, useDingbats = F)

Plot <- ggplot(DataRegionsAcetylationSum %>% filter(Region == "Cerebellum"), aes(Group, Mean)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Mean ratio") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group, color = Group), alpha = 0.4) +
  geom_jitter(width = 0.3, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group")  +
  facet_wrap(~Measure, nrow = 1) +
  geom_text(data = AcetylRegionStat %>% filter(Region == "Cerebellum") , aes(label = label), size = 5)
ggsave(paste0(ResultsPath,"WBacetylResultsCerebellum.pdf"), Plot, device = "pdf", width = 8, height = 3, dpi = 300, useDingbats = F)


#Drug cell reatment
DataDrugs <- read.table("data/PD_Drugs_cells_mixed_models_2.txt", header = T, sep = "\t") %>%
  gather(key = "Measure", value = "Value", matches("H3"))

DataDrugs$Sample <- factor(DataDrugs$Sample, levels =  DataDrugs$Sample %>% unique())
DataDrugs$Measure <- factor(DataDrugs$Measure, levels = c("H3K9K14ac", "H3K27ac","H3"))
levels(DataDrugs$Measure) <- c("H3K9/K14ac:H3", "H3K27ac:H3","H3:GAPDH")

Plot <- ggplot(DataDrugs) +
    theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid= element_blank()) +
  labs(x = "", y = "Ratio") +
  geom_bar(aes(Sample, Value), stat = "summary", fun.y = "mean", alpha = 0.2) +
  geom_jitter(aes(Sample, Value, color = Sample), height = 0, width = 0.2, size = 1, show.legend = F) +
  scale_color_manual(values = c("black", "red", "chartreuse4", "chartreuse3",
                                "cadetblue4", "cadetblue3", "aquamarine4","darkmagenta" ), name = "") +
  # scale_fill_manual(values = c("black", "coral4", "chartreuse4", "chartreuse3",
  #                               "cadetblue4", "cadetblue3", "aquamarine4","darkmagenta" ))
  facet_wrap(~Measure, nrow = 3)
ggsave(paste0(ResultsPath,"DrugTreatment.pdf"), Plot, device = "pdf", width = 6, height = 4, dpi = 300, useDingbats = F)
