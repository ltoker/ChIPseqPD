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

ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersPV$`respiratory chain`), aes(Group, Cor)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "respiratory chain genes", x = "") +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4)+
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(~Cohort, scales = "free")

MergedChIPrna %>% filter(GeneSymbol %in% GeneMembersPV$`mitochondrial translation`) %>% arrange(desc(CorCont.PV))



PVenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)





PVenrich2 <- gsr(scores = MergedChIPrna[MergedChIPrna$CorCont.PV > 0.5,], scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrich2 <- gsr(scores = MergedChIPrna %>% filter(CorCont.NBB > 0.6), scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)

PVenrichDown <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)
NBBenrichDown <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)


