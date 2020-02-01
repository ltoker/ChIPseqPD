ComplexIdef <- read.table("data/complex_I_counts_NDUFS4_frontalcortex.csv", sep = "\t", header = T, quote = "'", comment.char = "#")

PWmeta <- PV_DESeq2@colData %>% data.frame()

SampleMapping <- read.table("meta/RawMetata/SampleID_mapping.csv", sep = "\t", header = T, quote = "'", comment.char = "#")
ComplexIdef$activemotif_id <- SampleMapping$ChIPseq_id[match(ComplexIdef$ID, SampleMapping$Autopsy_id)]

PWmeta$ComplexIintact <- ComplexIdef$X..Positive[match(PWmeta$activemotif_id, ComplexIdef$activemotif_id)]
PWmeta$Autopsy_ID <- ComplexIdef$ID[match(PWmeta$activemotif_id, ComplexIdef$activemotif_id)]



Data1 <- read.table("data/PFC_acetylation_data_RESCALED_0-1.csv", header = T, sep = "\t")
DataWestern <- Data1 %>% select(ID, ID2, H3K27ac, SIRT1,   SIRT2,   SIRT3)
DataWestern <- merge(DataWestern, ComplexIdef, by.x = "ID2", by.y = "ID")
DataWestern$Group <- PWmeta$condition[match(DataWestern$ID2, PWmeta$Autopsy_ID)]

DataWesternMean <- DataWestern %>% group_by(ID2, Group) %>% summarise_if(is.numeric, mean, na.rm = T)
DataWesternMeanMelt <- gather(DataWesternMean, key = "Protein", value = "NormalizedIntensity", H3K27ac, SIRT1, SIRT2, SIRT3)

ggplot(DataWesternMeanMelt, aes(X..Positive, NormalizedIntensity)) +
  theme_classic() +
  geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", color = "black", size = 0.3) + 
  facet_wrap(~Protein, scales = "free_y")


ggplot(PWmeta, aes(ComplexIdef, c(RiP_NormMeanRatioOrg/Oligo_MSP))) + geom_point(aes(color = condition)) + geom_smooth(method = "lm")
ggplot(PWmeta, aes(ComplexIintact, RiP_NormMeanRatioOrg/Oligo_MSP)) + geom_point(aes(color = condition)) + geom_smooth(method = "lm")
write.table(PWmeta, "PWmetaComplexI.tsv", sep = "\t", row.names = F, col.names = T)


write.table(PWmeta %>%
              select(activemotif_id, Autopsy_ID, rin, condition, sex, age, batch, pm_hours, matches("MSP"), ComplexIdef),
            "PWmetaComplexI_2.tsv", sep = "\t", row.names = F, col.names = T)
