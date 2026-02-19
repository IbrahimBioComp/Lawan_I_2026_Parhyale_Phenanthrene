# ==============================================================================
# SCRIPT 3: FUNCTIONAL ENRICHMENT
# ==============================================================================
pacman::p_load(topGO, Rgraphviz, dplyr, stringr, ggplot2)
load("02_Results/analysis_workspace.RData") # Load stats

# 1. IMPORT ANNOTATIONS (Note "00_Raw_Data" path)
emapper_raw <- read.table("00_Raw_Data/out.emapper.annotations", sep="\t", header=FALSE, 
                          comment.char="#", quote="", fill=TRUE, stringsAsFactors=FALSE)

emapper_clean <- emapper_raw %>%
  dplyr::select(V1, V8, V10) %>%
  dplyr::rename(GeneID = V1, Description = V8, GO_terms = V10) %>%
  dplyr::filter(GO_terms != "-") %>%
  dplyr::distinct(GeneID, .keep_all = TRUE)

geneID2GO <- strsplit(as.character(emapper_clean$GO_terms), ",")
names(geneID2GO) <- emapper_clean$GeneID

# 2. RUN TOPGO
all_genes <- rownames(dds)
sig_genes_interaction <- rownames(subset(res_interaction, padj < 0.05))
geneList_int <- factor(as.integer(all_genes %in% sig_genes_interaction))
names(geneList_int) <- all_genes

GOdata_int <- new("topGOdata", ontology = "BP", allGenes = geneList_int, 
                  nodeSize = 10, annot = annFUN.gene2GO, gene2GO = geneID2GO)
result_weight <- runTest(GOdata_int, algorithm = "weight01", statistic = "fisher")
allRes_int <- GenTable(GOdata_int, raw_p = result_weight, topNodes = 20)
allRes_int$raw_p <- as.numeric(gsub("<", "", allRes_int$raw_p))

# 3. PLOT BARCHART (FIGURE S2)
allRes_int$Term <- factor(allRes_int$Term, levels = rev(allRes_int$Term)) 
p_go <- ggplot(head(allRes_int, 15), aes(x = Term, y = -log10(raw_p))) +
  geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
  coord_flip() + theme_bw() + labs(title="Figure S2: TopGO Enrichment")

ggsave("03_Figures/FigureS2_TopGO_Barplot.pdf", p_go, width = 8, height = 6)
ggsave("03_Figures/FigureS2_TopGO_Barplot.svg", p_go, width = 8, height = 6)

# 4. MECHANISM MINING (TABLE EXPORT)
# Helper function
get_genes_for_go <- function(go_id) {
  all_go <- genesInTerm(GOdata_int, go_id)[[1]]
  sig_go <- intersect(all_go, sig_genes_interaction)
  as.data.frame(res_female_phen[sig_go, ]) %>%
    dplyr::select(log2FoldChange, padj) %>%
    dplyr::mutate(GeneID = rownames(.)) %>%
    dplyr::left_join(emapper_clean, by="GeneID")
}

write.csv(get_genes_for_go("GO:0040003"), "02_Results/Table_S1_Chitin.csv")
write.csv(get_genes_for_go("GO:0008069"), "02_Results/Table_S2_Repro.csv")

print("Script 3 Complete.")
