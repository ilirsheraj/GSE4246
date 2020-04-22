#GSE4246 Data Analysis
#Note: Raw data was not available, so series matrix ones are used for this study
exp_matrix <- read.csv("series_matrix_all.csv")
ago2 <- c("ID_REF", "GSM96818", "GSM96819","GSM96816", "GSM96817")

ago2_exp <- exp_matrix %>% select(all_of(ago2))
summary(ago2_exp)
rownames(ago2_exp) <- ago2_exp$ID_REF
ago2_exp <- ago2_exp[,-c(1)]
ago2_exp_nonzero <- ago2_exp[rowSums(ago2_exp[1:3] > 0) > 0, ]
ago2_exp_nonzero_log <- log2(ago2_exp_nonzero + 1)

summary(ago2_exp_nonzero_log)
#    GSM96818         GSM96819         GSM96816         GSM96817     
# Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 3.463   1st Qu.: 3.407   1st Qu.: 3.428   1st Qu.: 3.201  
# Median : 5.817   Median : 5.805   Median : 5.804   Median : 5.724  
# Mean   : 5.719   Mean   : 5.651   Mean   : 5.673   Mean   : 5.447  
# 3rd Qu.: 8.136   3rd Qu.: 8.133   3rd Qu.: 8.142   3rd Qu.: 8.092  
# Max.   :14.062   Max.   :14.202   Max.   :14.075   Max.   :14.087  

boxplot(ago2_exp_nonzero_log, las = 2)
write.csv(ago2_exp_nonzero_log, file = "ago2_expression_log2.csv")

#Differential expression analysis for AGO2
library(limma)

design <- model.matrix(~0+as.factor(c(rep(1,2), rep(0,2))))
colnames(design) <- c("Control", "AGO_KD")
contrast.matrix <- makeContrasts(AGO_KD_v_Control = AGO_KD -Control, levels=design)
fit <- lmFit(ago2_exp_nonzero_log,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2, coef="AGO_KD_v_Control", number=nrow(fit2), adjust.method="fdr")
summary(decideTests(fit2))
#Down                637
#NotSig            25908
#Up                  784

write.csv(top, file = "Ago2_Control_diff_expression.csv")

#Check the expression of AGO2 Probes just in case
#4 AGO2 Probes
#1561042_at
#213310_at
#225569_at
#225827_at
top["1561042_at",]
#               logFC  AveExpr         t      P.Value  adj.P.Val         B
#1561042_at -4.130131 2.065066 -14.81944 0.0008008792 0.02260162 0.2941895
top["213310_at",]
#              logFC  AveExpr         t     P.Value  adj.P.Val          B
#213310_at -3.399559 5.888901 -12.94617 0.001179407 0.02911654 -0.1898804
top["225569_at",]
#              logFC  AveExpr         t     P.Value  adj.P.Val        B
#225569_at -2.147532 6.715038 -20.10881 0.000332919 0.01195577 1.353261
top["225827_at",]
#              logFC  AveExpr         t      P.Value  adj.P.Val        B
#225827_at -3.166302 6.521015 -24.14223 0.0001964764 0.00838985 1.953054

#Volcano Plot for all
with(top, plot(logFC, -log10(adj.P.Val), pch=20, cex = 0.8, 
               main= "AGO2 - Cont"))
with(subset(top, adj.P.Val <.05 ), points(logFC, -log10(adj.P.Val), 
                                          pch=20, cex = 0.8, col="red"))
with(subset(top, adj.P.Val <.05 & abs(logFC)>1), 
     points(logFC, -log10(adj.P.Val), pch=20, cex = 0.8, col="blue"))
#Volcano had to be truncated to on y-axis, so p-values up to 300 are not visible
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1, col="black", lty=3, lwd=1.0)
abline(v=1, col="black", lty=3, lwd=1.0)
abline(h = -log10(0.05))
legend("bottomleft", legend = "Control", text.font = 2, cex = 1.0, bty = "n")
legend("bottomright", legend = "AGO2_KD", text.font = 2, cex = 1.0, bty = "n")

###########################################################
#GO Enrichment to be done from here
ago2 <- read.csv("Ago2_Control_diff_expression.csv")
ago2_significant <- ago2 %>% filter(adj.P.Val < 0.05)
#1421 significant
hist(ago2_significant$logFC, breaks = 100)
abline(v=-1, col="black", lty=3, lwd=1.0)
abline(v=1, col="black", lty=3, lwd=1.0)
ago2_upregulated <- ago2_significant %>% filter(logFC > 0)
range(ago2_upregulated$logFC)
# 0.6891435 6.8579296
rownames(ago2_upregulated) <- ago2_upregulated$X
ago2_upregulated <- ago2_upregulated[,-c(1)]

ago2_downregulated <- ago2_significant %>% filter(logFC < 0)
range(ago2_downregulated$logFC)
#-7.7219150 -0.6657667
rownames(ago2_downregulated) <- ago2_downregulated$X
ago2_downregulated <- ago2_downregulated[,-c(1)]

detach(package:dplyr)
ago2_upregulated_ensemble <- select(hgu133plus2.db, keys=rownames(ago2_upregulated),
                                    columns=c("ENSEMBL", "SYMBOL","GENENAME"),
                                    keytype="PROBEID")
write.csv(ago2_upregulated, file = "ago2_upregulated_significant.csv")
ago2_upregulated_na <- ago2_upregulated_ensemble[!is.na(ago2_upregulated_ensemble$ENSEMBL), ]
ago2_upregulated_unique <- ago2_upregulated_na %>% distinct(PROBEID, .keep_all = TRUE)
ago2_up_ensemble_unique <- ago2_upregulated_unique %>% distinct(ENSEMBL, .keep_all = TRUE)
write.csv(ago2_up_ensemble_unique, file = "ago2_up_ensemble_unique.csv")


ago2_downregulated_ensemble <- select(hgu133plus2.db, keys=rownames(ago2_downregulated),
                                    columns=c("ENSEMBL", "SYMBOL","GENENAME"),
                                    keytype="PROBEID")
ago2_downregulated_na <- ago2_downregulated_ensemble[!is.na(ago2_downregulated_ensemble$ENSEMBL), ]
ago2_downregulated_unique <- ago2_downregulated_na %>% distinct(PROBEID, .keep_all = TRUE)
ago2_down_ensemble_unique <- ago2_downregulated_unique %>% distinct(ENSEMBL, .keep_all = TRUE)
write.csv(ago2_down_ensemble_unique, file = "ago2_down_ensemble_unique.csv")

#Enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)
detach(package:dplyr)
library(dplyr)

ago_up_bp <- enrichGO(gene         = ago2_up_ensemble_unique$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
#0 GO Terms
barplot(cms3_bp, drop=TRUE, showCategory=50)
write.csv(cms3_bp, file = "cms3_bp_goterm.csv")

ago_up_cc <- enrichGO(gene         = ago2_up_ensemble_unique$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
#0 GO Terms
barplot(cms3_cc, drop=TRUE, showCategory=50)
write.csv(cms3_cc, file = "cms3_cc_goterm.csv")

ago_down_bp <- enrichGO(gene         = ago2_down_ensemble_unique$ENSEMBL,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
#15 GO Terms
barplot(ago_down_bp, drop=TRUE, showCategory=50)
write.csv(ago_down_bp, file = "ago_down_bp_goterm.csv")

ago_down_cc <- enrichGO(gene         = ago2_down_ensemble_unique$ENSEMBL,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
#16 GO Terms
barplot(ago_down_cc, drop=TRUE, showCategory=50)
write.csv(ago_down_cc, file = "ago_down_cc_goterm.csv")

########################################################################
#Dicer Differential Expression
dicer <- c("ID_REF", "GSM96842", "GSM96843", "GSM96844", "GSM96845", "GSM96840", "GSM96841")
dicer_exp <- exp_matrix %>% select(all_of(dicer))
summary(dicer_exp)
rownames(dicer_exp) <- dicer_exp$ID_REF
dicer_exp <- dicer_exp[,-c(1)]
dicer_exp_nonzero <- dicer_exp[rowSums(dicer_exp[1:6] > 0) > 0, ]
dicer_exp_nonzero_log <- log2(dicer_exp_nonzero + 1)

summary(dicer_exp_nonzero_log)
# GSM96842         GSM96843         GSM96844         GSM96845         GSM96840         GSM96841     
# Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000  
# 1st Qu.: 2.266   1st Qu.: 2.472   1st Qu.: 2.164   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 2.152  
# Median : 4.422   Median : 4.444   Median : 4.356   Median : 4.300   Median : 4.344   Median : 4.391  
# Mean   : 4.621   Mean   : 4.718   Mean   : 4.594   Mean   : 4.527   Mean   : 4.542   Mean   : 4.593  
# 3rd Qu.: 7.434   3rd Qu.: 7.441   3rd Qu.: 7.435   3rd Qu.: 7.433   3rd Qu.: 7.416   3rd Qu.: 7.420  
# Max.   :14.142   Max.   :14.195   Max.   :14.279   Max.   :14.174   Max.   :14.221   Max.   :14.171  

boxplot(dicer_exp_nonzero_log, las = 2)
write.csv(dicer_exp_nonzero_log, file = "dicer_expression_log2.csv")

design <- model.matrix(~0+as.factor(c(rep(1,4), rep(0,2))))
colnames(design) <- c("Control", "Dicer_KD")
contrast.matrix <- makeContrasts(Dicer_KD_v_Control = Dicer_KD -Control, levels=design)
fit <- lmFit(dicer_exp_nonzero_log,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2, coef="Dicer_KD_v_Control", number=nrow(fit2), adjust.method="fdr")
summary(decideTests(fit2))
#       Dicer_KD_v_Control
#Down                  181
#NotSig              34981
#Up                    247

write.csv(top, file = "Dicer_KD_Control_differential_expression.csv")

#Dicer probes
#206061_s_at
#212888_at
#213229_at
#216260_at
#216280_s_at
#216281_at
top["206061_s_at",]
#                 logFC  AveExpr         t     P.Value adj.P.Val         B
# 206061_s_at -2.741087 6.383193 -4.684894 0.005087655 0.2283254 -2.278999
top["212888_at",]
#               logFC  AveExpr         t      P.Value  adj.P.Val         B
# 212888_at -2.655847 7.950367 -8.365138 0.0003552331 0.03090528 0.7253557
top["213229_at",]
#              logFC  AveExpr         t      P.Value   adj.P.Val        B
#213229_at -1.925008 8.842109 -11.20806 8.493913e-05 0.008195122 2.290237
top["216260_at",]
# 0 signal, filtered out
top["216280_s_at",]
#                logFC   AveExpr         t   P.Value adj.P.Val         B
#216280_s_at -0.9912107 0.9128815 -0.858092 0.4291584 0.7227735 -6.910033
top["216281_at",]
# 0 Signal, filtered out
#Volcano Plot for all
with(top, plot(logFC, -log10(adj.P.Val), pch=20, cex = 0.8, 
               main= "Dicer - Cont"))
with(subset(top, adj.P.Val <.05 ), points(logFC, -log10(adj.P.Val), 
                                          pch=20, cex = 0.8, col="red"))
with(subset(top, adj.P.Val <.05 & abs(logFC)>1), 
     points(logFC, -log10(adj.P.Val), pch=20, cex = 0.8, col="blue"))
#Volcano had to be truncated to on y-axis, so p-values up to 300 are not visible
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1, col="black", lty=3, lwd=1.0)
abline(v=1, col="black", lty=3, lwd=1.0)
abline(h = -log10(0.05))
legend("bottomleft", legend = "Control", text.font = 2, cex = 1.0, bty = "n")
legend("bottomright", legend = "Dicer_KD", text.font = 2, cex = 1.0, bty = "n")

##################################################
#GO Enrichment
dicer <- read.csv("Dicer_KD_Control_differential_expression.csv")
dicer_significant <- dicer %>% filter(adj.P.Val < 0.05)
#428
hist(dicer_significant$logFC, breaks = 100)
abline(v=-1, col="black", lty=3, lwd=1.0)
abline(v=1, col="black", lty=3, lwd=1.0)

dicer_upregulated <- dicer_significant %>% filter(logFC > 0)
#247
rownames(dicer_upregulated) <- dicer_upregulated$X
dicer_upregulated <- dicer_upregulated[,-c(1)]
write.csv(dicer_upregulated, file = "dicer_upregulated_significant.csv")

dicer_downregulated <- dicer_significant %>% filter(logFC < 0)
#181
rownames(dicer_downregulated) <- dicer_downregulated$X
dicer_downregulated <- dicer_downregulated[,-c(1)]
write.csv(dicer_downregulated, file = "dicer_downregulated_significant.csv")

detach(package:dplyr)
dicer_upregulated_ensemble <- select(hgu133plus2.db, keys=rownames(dicer_upregulated),
                                    columns=c("ENSEMBL", "SYMBOL","GENENAME"),
                                    keytype="PROBEID")

dicer_upregulated_na <- dicer_upregulated_ensemble[!is.na(dicer_upregulated_ensemble$ENSEMBL), ]
dicer_upregulated_unique <- dicer_upregulated_na %>% distinct(PROBEID, .keep_all = TRUE)
dicer_up_ensemble_unique <- dicer_upregulated_unique %>% distinct(ENSEMBL, .keep_all = TRUE)
write.csv(dicer_up_ensemble_unique, file = "dicer_up_ensemble_unique.csv")


dicer_downregulated_ensemble <- select(hgu133plus2.db, keys=rownames(dicer_downregulated),
                                      columns=c("ENSEMBL", "SYMBOL","GENENAME"),
                                      keytype="PROBEID")
dicer_downregulated_na <- dicer_downregulated_ensemble[!is.na(dicer_downregulated_ensemble$ENSEMBL), ]
dicer_downregulated_unique <- dicer_downregulated_na %>% distinct(PROBEID, .keep_all = TRUE)
dicer_down_ensemble_unique <- dicer_downregulated_unique %>% distinct(ENSEMBL, .keep_all = TRUE)
write.csv(dicer_down_ensemble_unique, file = "dicer_down_ensemble_unique.csv")

dicer_up_bp <- enrichGO(gene         = dicer_up_ensemble_unique$ENSEMBL,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
#0 GO Terms

dicer_up_cc <- enrichGO(gene         = dicer_up_ensemble_unique$ENSEMBL,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
#0 GO Terms

dicer_down_bp <- enrichGO(gene         = dicer_down_ensemble_unique$ENSEMBL,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
#0 GO Terms

dicer_down_cc <- enrichGO(gene         = dicer_down_ensemble_unique$ENSEMBL,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
#0 GO Terms
