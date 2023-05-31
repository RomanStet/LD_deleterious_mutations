setwd("/shared/projects/selfrecomb/Lise/variant_calling_Cg/VariantFiltration")

library("gridExtra")
library("ggplot2")

VCFsnps <- read.csv("cr145_allScaffolds_snps_raw.table", header = T, na.strings=c("","NA"), sep = "\t") 
dim(VCFsnps)
VCF <- rbind(VCFsnps)
VCF$Variant <- factor(rep("SNPs", dim(VCFsnps)[1]))

snps <- "#A9E2E4"

DP <- ggplot(VCF, aes(x=DP, fill=Variant)) + geom_density(alpha=.3) + xlim(0,30000) + geom_vline(xintercept=c(0,30000), size=1) + theme(legend.position="none")

QD <- ggplot(VCF, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) + xlim(0,45) + geom_vline(xintercept=2, size=1) + theme(legend.position="none")

FS <- ggplot(VCF, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) + xlim(0,100) + geom_vline(xintercept=c(60, 200), size=1) + ylim(0,0.08) + theme(legend.position="none")

MQ <- ggplot(VCF, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) + xlim(15,65) + geom_vline(xintercept=50, size=1) + theme(legend.position="none")

MQRankSum <- ggplot(VCF, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3) + xlim(-8,8) + geom_vline(xintercept=-5, size=1) + theme(legend.position="none")

SOR <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) + xlim(0,10) + geom_vline(xintercept=3, size=1) + theme(legend.position="none")

ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) + geom_vline(xintercept=c(-8,8), size=1) + xlim(-15, 15) + theme(legend.position="none")

svg("Variant_Scores_Cr145.svg", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow=4)
dev.off()
