# ==============================================================================
# SCRIPT 4: FIGURE GENERATION
# ==============================================================================
pacman::p_load(ggplot2, gridExtra, grid, pheatmap, svglite, cowplot, DESeq2, tidyverse)
load("02_Results/analysis_workspace.RData") # Load stats
dds <- readRDS("02_Results/dds_object.rds") # Load counts

# Colors
col_sex <- c("Male" = "#2E86AB", "Female" = "#A23B72")
col_treat <- c("Seawater" = "#009E73", "DMSO" = "#999999", "Phenanthrene" = "#E69F00")
col_heatmap <- colorRampPalette(c("#2D004B", "#542788", "#F7F7F7", "#B35806", "#7F3B08"))(100)

# ------------------------------------------------------------------------------
# FIGURE 6B: DUAL VOLCANO
# ------------------------------------------------------------------------------
prepare_volcano <- function(res, sex) {
  as.data.frame(res) %>% rownames_to_column("GeneID") %>%
    mutate(Sex = sex, padj = ifelse(is.na(padj), 1, padj), logP = -log10(padj),
           Significance = case_when(padj < 0.05 & log2FoldChange > 1 ~ "Up",
                                    padj < 0.05 & log2FoldChange < -1 ~ "Down", TRUE ~ "NS"))
}
df_all <- bind_rows(prepare_volcano(res_male_phen, "Male"), prepare_volcano(res_female_phen, "Female")) %>%
  mutate(Sex = factor(Sex, levels = c("Male", "Female")))

max_y <- max(df_all$logP, na.rm=T) * 1.05
max_x <- max(abs(df_all$log2FoldChange), na.rm=T) * 1.05

p_volc <- ggplot(df_all, aes(x=log2FoldChange, y=logP, color=Significance)) +
  geom_point(alpha=0.7, size=1.5) +
  scale_color_manual(values=c("Down"="blue", "NS"="grey90", "Up"="firebrick")) +
  facet_grid(.~Sex) + theme_bw(base_size=14) +
  coord_cartesian(ylim=c(0, max_y), xlim=c(-max_x, max_x)) +
  labs(title="Figure 6B: Dual Volcano")

ggsave("03_Figures/Figure6B_Dual_Volcano.pdf", p_volc, width=10, height=6)
ggsave("03_Figures/Figure6B_Dual_Volcano.svg", p_volc, width=10, height=6)

# ------------------------------------------------------------------------------
# FIGURE 7: INTERACTION HEATMAP
# ------------------------------------------------------------------------------
vsd <- vst(dds, blind=FALSE)
sig_genes <- rownames(subset(res_interaction, padj < 0.05))
mat <- assay(vsd)[sig_genes, ]
mat <- mat - rowMeans(mat)
meta_ord <- as.data.frame(colData(dds)) %>% arrange(Sex, Treatment)
mat_ord <- mat[, rownames(meta_ord)]
df_anno <- meta_ord[, c("Sex", "Treatment")]
ann_col <- list(Sex=col_sex, Treatment=col_treat)

pdf("03_Figures/Figure7_Interaction_Heatmap.pdf", width=8, height=10)
pheatmap(mat_ord, color=col_heatmap, annotation_col=df_anno, annotation_colors=ann_col,
         scale="row", cluster_cols=FALSE, show_rownames=FALSE,
         gaps_col=which(diff(as.numeric(as.factor(meta_ord$Sex))) != 0))
dev.off()

svg("03_Figures/Figure7_Interaction_Heatmap.svg", width=8, height=10)
pheatmap(mat_ord, color=col_heatmap, annotation_col=df_anno, annotation_colors=ann_col,
         scale="row", cluster_cols=FALSE, show_rownames=FALSE,
         gaps_col=which(diff(as.numeric(as.factor(meta_ord$Sex))) != 0))
dev.off()

# ------------------------------------------------------------------------------
# FIGURE 8: BOXPLOTS
# ------------------------------------------------------------------------------
ids <- c("phaw_30_tra_m.009099", "phaw_30_tra_m.012099", "phaw_30_tra_m.006274", "phaw_30_tra_m.004007")
titles <- c("A. Chitin", "B. Ovarian", "C. NPF", "D. Trypsin")
plots <- list()

for(i in 1:4) {
  d <- plotCounts(dds, gene=ids[i], intgroup=c("Treatment", "Sex"), returnData=TRUE)
  plots[[i]] <- ggplot(d, aes(x=Treatment, y=count, fill=Sex)) +
    geom_boxplot(alpha=0.8, outlier.shape=NA) + geom_point(position=position_jitterdodge()) +
    scale_y_log10() + scale_fill_manual(values=col_sex) + theme_bw() + 
    labs(title=titles[i], x=NULL) + theme(legend.position="none")
}

p_final <- grid.arrange(arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=2))
ggsave("03_Figures/Figure8_Boxplots.pdf", p_final, width=9, height=8)
ggsave("03_Figures/Figure8_Boxplots.svg", p_final, width=9, height=8)

print("Script 4 Complete. All Figures Generated.")
