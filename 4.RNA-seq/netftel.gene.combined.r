# Load required library
library(ggplot2)
library(dplyr)
library(ggpubr)

# Define the file names and their corresponding subtypes
files <- c("single.mes_new1", "single.ac-1", "single.mes_new2", "single.npc_new1", "single.npc_new2", "single.opc_new")
subtypes <- c("MES_NEW", "AC_NEW", "MES_NEW", "NPC_NEW", "NPC_NEW", "OPC_NEW")

# Define the path to the base directory for the input files
base_dir <- "/proj/nobackup/sens2023512/wharf/lucycy/lucycy-sens2023512/database/list/single/"

# Read the gene count matrix
data <- read.table("/proj/nobackup/sens2023512/wharf/lucycy/lucycy-sens2023512/glioblastoma_single_cell/10.osmr_rna_seq2/gene_count_matrix.clean.cpm.log.csv", header = T, sep = "\t")

# Define the group names
new_names <- c(
  "P13ACTRL", "P13AOSM", "P13ECTRL", 
  "P15ACTRL", "P15AOSM", "P15ECTRL", 
  "P15EOSM", "P17ACTRL", "P17AOSM", 
  "P17ECTRL", "P17EOSM", "P1E3OSM"
)

# Initialize an empty data frame for combining data
my_data <- data.frame(
  group = character(),
  rpkm = numeric(),
  subtype = character()
)

# Loop over each file to process and append data to `my_data`
for (i in seq_along(files)) {
  # Read the specific subtype file
  file_path <- paste0(base_dir, files[i])
  subtype_data <- read.table(file_path, header = F, sep = "\t")
  
  # Filter data based on the current subtype
  filtered_data <- data[as.character(subtype_data[,1]),]
  
  # Reshape the filtered data into long format
  filtered_data_long <- data.frame(
    group = rep(new_names, each = nrow(filtered_data)),
    rpkm = as.vector(as.matrix(filtered_data)),
    subtype = subtypes[i]
  )
  
  # Append the data to `my_data`
  my_data <- rbind(my_data, filtered_data_long)
}

# Convert `subtype` column to a factor with specified levels
my_data$subtype <- factor(my_data$subtype, levels = c("NPC_NEW", "OPC_NEW", "AC_NEW", "MES_NEW"))
my_data <- my_data[my_data$group %in% c("P13ECTRL", "P15ECTRL", "P15EOSM", "P17ECTRL", "P17EOSM", "P1E3OSM"),]
my_data$group <- factor(my_data$group, levels =c(  "P13ECTRL","P1E3OSM",
  "P15ECTRL","P15EOSM",
  "P17ECTRL", "P17EOSM"
))

custom_colors <- c(
  "NPC_NEW" = "#3b6a7f",  # Corrected hex code without trailing comma
  "OPC_NEW" = "#4683a0",   # Corrected hex code without trailing comma
  "AC_NEW" = "#c65250",    # Corrected hex code without trailing comma
  "MES_NEW" = "#aa2e36"  # Corrected hex code without trailing comma
)

# Plot the combined data using ggplot2
p <- ggplot(my_data, aes(x = group, y = rpkm, fill = subtype)) +
  geom_boxplot(outlier.shape = NA) + # Avoid plotting outliers if not needed
  scale_fill_manual(values = custom_colors) + # Use manual colors for subtypes
  theme_classic() +
  labs(
    title = "Combined Barplot of Subtypes",
    x = "Group",
    y = "CPM",
    fill = "Subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for readability
    text = element_text(size = 14)
  )

# Save the plot to a PDF
pdf("combined_barplot.combined.pdf", width = 10, height = 6)
print(p)
dev.off()

###compute the significance of MES new ####
mes_only <- my_data[my_data$subtype == "MES_NEW",]

# 2) Define the per-sample CTRL vs OSM comparisons (match your factor levels)
cmp <- list(
  c("P13ECTRL", "P1E3OSM"),
  c("P15ECTRL", "P15EOSM"),
  c("P17ECTRL", "P17EOSM")
)

# 3) Plot + Wilcoxon (rank-sum) comparisons
p_mes <- ggplot(mes_only, aes(x = group, y = rpkm, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  theme_classic() +
  labs(
    title = "MES_NEW: CTRL vs OSM per Sample (Wilcoxon rank-sum)",
    x = "Group",
    y = "log2(CPM)",
    fill = "Subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    legend.position = "none"  # only MES_NEW is shown anyway
  ) +
  stat_compare_means(
    comparisons = cmp,
    method = "wilcox.test",
    label = "p.signif",       # shows significance stars; use "p.format" for exact p
    step.increase = 0.12,     # stacks brackets nicely
    tip.length = 0.01
  )

pdf("mes_new.samples.pdf")
print (p_mes)
dev.off()

###compute the significance of MES new ####
npc_only <- my_data[my_data$subtype == "NPC_NEW",]

# 2) Define the per-sample CTRL vs OSM comparisons (match your factor levels)
cmp <- list(
  c("P13ECTRL", "P1E3OSM"),
  c("P15ECTRL", "P15EOSM"),
  c("P17ECTRL", "P17EOSM")
)

# 3) Plot + Wilcoxon (rank-sum) comparisons
p_npc <- ggplot(npc_only, aes(x = group, y = rpkm, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  theme_classic() +
  labs(
    title = "NPC_NEW: CTRL vs OSM per Sample (Wilcoxon rank-sum)",
    x = "Group",
    y = "log2(CPM)",
    fill = "Subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    legend.position = "none"  # only MES_NEW is shown anyway
  ) +
  stat_compare_means(
    comparisons = cmp,
    method = "wilcox.test",
    label = "p.signif",       # shows significance stars; use "p.format" for exact p
    step.increase = 0.12,     # stacks brackets nicely
    tip.length = 0.01
  )

pdf("NPC_NEW.samples.pdf")
print (p_npc)
dev.off()

###compute the significance of OPC_NEW ####
opc_only <- my_data[my_data$subtype == "OPC_NEW",]

# 2) Define the per-sample CTRL vs OSM comparisons (match your factor levels)
cmp <- list(
  c("P13ECTRL", "P1E3OSM"),
  c("P15ECTRL", "P15EOSM"),
  c("P17ECTRL", "P17EOSM")
)

# 3) Plot + Wilcoxon (rank-sum) comparisons
p_opc <- ggplot(opc_only, aes(x = group, y = rpkm, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  theme_classic() +
  labs(
    title = "OPC_NEW: CTRL vs OSM per Sample (Wilcoxon rank-sum)",
    x = "Group",
    y = "log2(CPM)",
    fill = "Subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    legend.position = "none"  # only MES_NEW is shown anyway
  ) +
  stat_compare_means(
    comparisons = cmp,
    method = "wilcox.test",
    label = "p.signif",       # shows significance stars; use "p.format" for exact p
    step.increase = 0.12,     # stacks brackets nicely
    tip.length = 0.01
  )

pdf("OPC_NEW.samples.pdf")
print (p_opc)
dev.off()

###compute the significance of AC_NEW ####
ac_only <- my_data[my_data$subtype == "AC_NEW",]

# 2) Define the per-sample CTRL vs OSM comparisons (match your factor levels)
cmp <- list(
  c("P13ECTRL", "P1E3OSM"),
  c("P15ECTRL", "P15EOSM"),
  c("P17ECTRL", "P17EOSM")
)

# 3) Plot + Wilcoxon (rank-sum) comparisons
p_ac <- ggplot(ac_only, aes(x = group, y = rpkm, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  theme_classic() +
  labs(
    title = "ac_NEW: CTRL vs OSM per Sample (Wilcoxon rank-sum)",
    x = "Group",
    y = "log2(CPM)",
    fill = "Subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    legend.position = "none"  # only MES_NEW is shown anyway
  ) +
  stat_compare_means(
    comparisons = cmp,
    method = "wilcox.test",
    label = "p.signif",       # shows significance stars; use "p.format" for exact p
    step.increase = 0.12,     # stacks brackets nicely
    tip.length = 0.01
  )

pdf("ac_NEW.samples.pdf")
print (p_ac)
dev.off()

