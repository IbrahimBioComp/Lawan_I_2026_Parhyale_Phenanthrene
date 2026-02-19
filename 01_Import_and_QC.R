# ==============================================================================
# PROJECT: Parhyale Ecotoxicology (Sex x Treatment Interaction)
# AUTHOR: Ibrahim Lawan
# DATE: 25-01-2026
# GOAL: Robust Import, QC, and Interaction Modelling for T24
# ==============================================================================
##

##
# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
# We use 'pacman' to ensure libraries are installed and loaded efficiently.
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,    # Data manipulation
  DESeq2,       # Differential expression
  pheatmap,     # Heatmaps
  RColorBrewer, # Palettes
  viridis,      # Accessibility palettes
  ggplot2,      # Plotting
  gridExtra,    # Arranging plots
  grid,         # Graphics control
  svglite       # Vector graphics export
)
# ==============================================================================
# SCRIPT 1: IMPORT & QC
# ==============================================================================

# Global Colors (Used in PCA)
col_treat <- c("Seawater" = "#009E73", "DMSO" = "#999999", "Phenanthrene" = "#E69F00")

# 2. IMPORT METADATA (Note the "00_Raw_Data" path)
raw_meta <- read.table("00_Raw_Data/sample_metadata.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE, stringsAsFactors = FALSE)

colnames(raw_meta) <- gsub("Time_Point \\(h\\)", "TimePoint", colnames(raw_meta))
colnames(raw_meta) <- gsub(" ", "_", colnames(raw_meta))
meta_t24 <- raw_meta %>% filter(TimePoint == 24)

# 3. CLEAN FACTORS
meta_t24_clean <- meta_t24 %>%
  mutate(
    Treatment_Clean = case_when(
      grepl("Phen", Treatment) ~ "Phenanthrene",
      grepl("DMSO", Treatment) ~ "DMSO",
      grepl("Seawater|Control", Treatment) ~ "Seawater",
      TRUE ~ "Unknown" 
    ),
    Batch_Status = ifelse(Replacement == 2, "Backup", "Original")
  )

if(any(meta_t24_clean$Treatment_Clean == "Unknown")) stop("Error: Unknown treatment.")

meta_t24_clean$Treatment <- factor(meta_t24_clean$Treatment_Clean, 
                                   levels = c("DMSO", "Seawater", "Phenanthrene"))
meta_t24_clean$Sex <- factor(meta_t24_clean$Sex, levels = c("Male", "Female"))
rownames(meta_t24_clean) <- meta_t24_clean$Sample_ID

# 4. MATCH FILES (Note the "00_Raw_Data/htseq" path)
files_list <- list.files("00_Raw_Data/htseq", pattern = "\\.htseq$", full.names = TRUE)

file_matcher <- data.frame(Full_Path = files_list) %>%
  mutate(Extracted_Base = str_extract(basename(Full_Path), "IL-?[0-9]+"),
         Clean_ID = str_replace(Extracted_Base, "-", ""))

meta_final <- meta_t24_clean %>%
  inner_join(file_matcher, by = c("Sample_ID" = "Clean_ID"))

# 5. DESEQ2 OBJECT
sampleTable_full <- data.frame(
  sampleName   = meta_final$Sample_ID,
  fileName     = basename(meta_final$Full_Path),
  Sex          = meta_final$Sex,
  Treatment    = meta_final$Treatment,
  Batch_Status = meta_final$Batch_Status
)

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_full,
  directory   = "00_Raw_Data/htseq", # CRITICAL PATH UPDATE
  design      = ~ Sex + Treatment + Sex:Treatment
)
colData(dds) <- cbind(colData(dds), meta_final)

# Filter (>=10 reads in >=3 samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# 6. GENERATE PCA (FIGURE 6A)
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("Treatment", "Sex"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = Treatment, shape = Sex)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = col_treat) +
  scale_shape_manual(values = c(19, 17)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  coord_fixed(ratio = percentVar[2] / percentVar[1]) +
  ggtitle("Figure 6A: PCA")

# Save Figure
ggsave("03_Figures/Figure6A_PCA.pdf", p_pca, width = 7, height = 5)
ggsave("03_Figures/Figure6A_PCA.svg", p_pca, width = 7, height = 5)

# 7. SAVE DATA FOR NEXT SCRIPT
saveRDS(dds, "02_Results/dds_object.rds")
print("Script 1 Complete. Data saved to 02_Results.")
