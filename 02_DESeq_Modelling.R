# ==============================================================================
# SCRIPT 2: DIFFERENTIAL EXPRESSION
# ==============================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(DESeq2, tidyverse)

# 1. LOAD DATA
dds <- readRDS("02_Results/dds_object.rds")

# 2. RUN DESEQ
dds <- DESeq(dds)

# 3. EXTRACT RESULTS
# A. Interaction (Sex Diff)
res_interaction <- results(dds, name="SexFemale.TreatmentPhenanthrene", alpha=0.05)
sig_interaction <- subset(res_interaction, padj < 0.05)

# B. Female Response (Phen vs DMSO)
res_female_phen <- results(dds, contrast=list(c("Treatment_Phenanthrene_vs_DMSO", "SexFemale.TreatmentPhenanthrene")), alpha=0.05)
sig_female <- subset(res_female_phen, padj < 0.05)

# C. Male Response (Phen vs DMSO)
res_male_phen <- results(dds, name="Treatment_Phenanthrene_vs_DMSO", alpha=0.05)

# D. Solvent Check (DMSO vs Seawater)
res_solvent <- results(dds, contrast=c("Treatment", "DMSO", "Seawater"))

# 4. EXPORT TABLES (Note "02_Results" path)
write.csv(as.data.frame(sig_interaction), "02_Results/Results_Interaction_Candidates.csv")
write.csv(as.data.frame(sig_female), "02_Results/Results_Female_Phen_Response.csv")

# 5. AUDIT (The Supervisor Check)
audit_genes <- c("phaw_30_tra_m.009099", "phaw_30_tra_m.012099", "phaw_30_tra_m.006274", "phaw_30_tra_m.004007", "phaw_30_tra_m.011201")
print("--- AUDIT: PHEN vs DMSO ---")
print(as.data.frame(res_female_phen[audit_genes, c("log2FoldChange", "padj")]))

# 6. SAVE WORKSPACE FOR PLOTTING
save(dds, res_interaction, res_female_phen, res_male_phen, res_solvent, file = "02_Results/analysis_workspace.RData")
print("Script 2 Complete. Results saved.")
