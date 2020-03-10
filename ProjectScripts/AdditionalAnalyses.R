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

ggplot(MergedChIPrnaMelt %>% filter(GeneSymbol %in% GeneMembersPV$`mitochondrial respiratory chain complex I assembly`), aes(Group, Cor)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "mitochondrial respiratory chain complex I assembly", x = "") +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4)+
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(~Cohort, scales = "free")


PVenrichCont <- ora(threshold = 0.8, scores = MergedChIPrna, scoreColumn = "CorCont.PV", aspects = c("M", "B"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrichCont <- ora(threshold = 0.8, scores = MergedChIPrna, scoreColumn = "CorCont.NBB", aspects = c("M", "B"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)



PVenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)
NBBenrichRLE <- gsr(scores = MergedChIPrnaRLE, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = T, logTrans = F)





PVenrichDown <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.PV", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)
NBBenrichDown <- gsr(scores = MergedChIPrna, scoreColumn = "Delta.NBB", aspects = c("M", "B", "C"), annotation = HumanAnno, bigIsBetter = F, logTrans = F)


