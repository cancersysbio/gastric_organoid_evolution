#############################################################################################################################
##                                                                                                                      
##  SCRIPT TO GENERATE VOLCANO PLOTS FOR ECB SAMPLES
##                                                                                                                      
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


library("reshape2")
library("stringr")
library("readr")
library("tidyverse")
library("ggpubr")
library("ggplot2")
library("RColorBrewer")
library("ggsci")
library("ggrepel")
library("dplyr")



############### D2C2R2 - FIGURE 6F #####################

############## READ IN FILE WITH DIFFERENTIAL EXPRESSION DATA
setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Data_objects/Seq8_D2C2R2T2"))
results <- read.csv("DEGs_allGenes_againstRG0/markers_RG_0a_vs_0b.txt",sep="\t")


############## SET PARAMETERS FOR VOLCANO PLOT
colnames(results) <- c("p_value","log2FoldChange","pct.1","pct.2","p_val_adj")
results$color <- "NS"
results[abs(results$log2FoldChange) > 1.5, "color"] <- "Log2FC"
results[results$p_value < 1e-2, "color"] <- "p-value"
results[abs(results$log2FoldChange) > 1.5 & results$p_value < 1e-2, "color"] <- "Log2FC & p-value"
results$color <- factor(results$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))
col.values <- c("lightgrey", "blue","red","green")
options(ggrepel.max.overlaps = Inf)
results$hgnc_symbol <- rownames(results)


############## PLOT FIGURE
ggplot(results, aes(x=log2FoldChange, y = -log10(p_value))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-1.5, 1.5), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(1e-2), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  #xlim(-2.5, 2.5) + ## THIS NEEDS TO BE ADAPTED TOO
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = results %>% filter(p_value < 1e-2 , log2FoldChange <= -1.5), # HERE YOU CAN SPECIFY THE GENE NAMES THAT SHOULD BE LABELLED
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 4) + 
  geom_text_repel(data = results %>% filter(p_value < 1e-2 , log2FoldChange >= 1.5),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 4)

ggsave("Volcano_0a_vs_0b.pdf", width = 6, height = 6, dpi = 600)




############### D2C1R1 - FIGURE ED7F #####################


setwd("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_6/Volcano/DEGs_special_D2C1")
results <- read.csv("markers_R1T7_0a_vs_0b.txt",sep="\t")
head(results)

colnames(results) <- c("p_value","log2FoldChange","pct.1","pct.2","p_val_adj")

results$color <- "NS"
results[abs(results$log2FoldChange) > 1, "color"] <- "Log2FC"
results[results$p_value < 1e-5, "color"] <- "p-value"
results[abs(results$log2FoldChange) > 1 & results$p_value < 1e-5, "color"] <- "Log2FC & p-value"
results$color <- factor(results$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey","red","green") 

options(ggrepel.max.overlaps = Inf)

#results %>% filter(p_value < 0.01 , log2FoldChange <= -1.5)
results$hgnc_symbol <- rownames(results)

ggplot(results, aes(x=log2FoldChange, y = -log10(p_value))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-1, 1), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(1e-5), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  xlim(-2.5, 2.5) + ## THIS NEEDS TO BE ADAPTED TOO
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = results %>% filter(p_value < 1e-5 , log2FoldChange <= -1), # HERE YOU CAN SPECIFY THE GENE NAMES THAT SHOULD BE LABELLED
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 4) + 
  geom_text_repel(data = results %>% filter(p_value < 1e-5 , log2FoldChange >= 1),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 4)

ggsave("Volcano_R1T7_0a_vs_0b.pdf", width = 6, height = 6, dpi = 600)


############### D2C1R2 - FIGURE ED8E #####################


setwd("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_6/Volcano/DEGs_special_D2C1")
results <- read.csv("markers_R2T7_1a_vs_1b.txt",sep="\t")
head(results)

colnames(results) <- c("p_value","log2FoldChange","pct.1","pct.2","p_val_adj")

results$color <- "NS"
results[abs(results$log2FoldChange) > 1, "color"] <- "Log2FC"
results[results$p_value < 1e-5, "color"] <- "p-value"
results[abs(results$log2FoldChange) > 1 & results$p_value < 1e-5, "color"] <- "Log2FC & p-value"
results$color <- factor(results$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey","red","green")

options(ggrepel.max.overlaps = Inf)

#results %>% filter(p_value < 0.01 , log2FoldChange <= -1.5)
results$hgnc_symbol <- rownames(results)

ggplot(results, aes(x=log2FoldChange, y = -log10(p_value))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-1, 1), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(1e-5), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  xlim(-2.5, 2.5) + ## THIS NEEDS TO BE ADAPTED TOO
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = results %>% filter(p_value < 1e-5 , log2FoldChange <= -1), # HERE YOU CAN SPECIFY THE GENE NAMES THAT SHOULD BE LABELLED
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 4) + 
  geom_text_repel(data = results %>% filter(p_value < 1e-5 , log2FoldChange >= 1),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 4)

ggsave("Volcano_R2T7_1a_vs_1b.pdf", width = 6, height = 6, dpi = 600)



############### D3C2R1 - FIGURE ED9E #####################


setwd("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_6/Volcano/DEGs_special_D3C2")
results <- read.csv("markers_R1_0a_vs_0b.txt",sep="\t")
head(results)

colnames(results) <- c("p_value","log2FoldChange","pct.1","pct.2","p_val_adj")

results$color <- "NS"
results[abs(results$log2FoldChange) > 1, "color"] <- "Log2FC"
results[results$p_value < 1e-5, "color"] <- "p-value"
results[abs(results$log2FoldChange) > 1 & results$p_value < 1e-5, "color"] <- "Log2FC & p-value"
results$color <- factor(results$color, levels = c("NS", "Log2FC", "p-value", "Log2FC & p-value"))

col.values <- c("lightgrey", "blue","red","green")

options(ggrepel.max.overlaps = Inf)

#results %>% filter(p_value < 0.01 , log2FoldChange <= -1.5)
results$hgnc_symbol <- rownames(results)

ggplot(results, aes(x=log2FoldChange, y = -log10(p_value))) +
  geom_point(aes(color = color), size = 1) +
  geom_vline(xintercept=c(-1, 1), linetype="dotted", color = "black") +
  geom_hline(yintercept = -log10(1e-5), linetype="dotted", color = "black") +
  theme_classic() +
  # geom_smooth(method=lm, color="darkblue") +
  labs(x="Log2FC",
       y="-Log10(p-value)") + 
  #xlim(-2.5, 2.5) + ## THIS NEEDS TO BE ADAPTED TOO
  theme(legend.position="bottom") +
  scale_color_manual(values = col.values) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold" ),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=8, face="bold")) +
  geom_text_repel(data = results %>% filter(p_value < 1e-5 , log2FoldChange <= -1), # HERE YOU CAN SPECIFY THE GENE NAMES THAT SHOULD BE LABELLED
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = -0.75,
                  direction = "y",
                  fontface = "bold",
                  size = 4) + 
  geom_text_repel(data = results %>% filter(p_value < 1e-5 , log2FoldChange >= 1),
                  aes(label = paste0(hgnc_symbol)) ,
                  hjust = -.35,
                  nudge_x = 0.5,
                  direction = "y",
                  fontface = "bold",
                  size = 4)

ggsave("Volcano_0a_vs_0b.pdf", width = 6, height = 6, dpi = 600)
