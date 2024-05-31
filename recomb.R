setwd("/home/amaniouloux/Documents/data/dataframe")
table_albicollis<-read.table("table_albicollis.txt", header = T)
table_taiga<-read.table("table_taiga.txt", header = T)

library(ggplot2)

table_albicollis_recomb <- table_albicollis[, c("gene", "chr","length","category","new_category","gene_recom_coll","recom_bin_coll","win_cons_pc")]
table_albicollis_recomb<-na.omit(table_albicollis_recomb)
table_albicollis_recomb <- subset(table_albicollis_recomb, chr != "ChrLGE22")
table_albicollis_recomb <- subset(table_albicollis_recomb, chr != "Chr25")

table_taiga_recomb <- table_taiga[, c("gene", "chr","length","category","new_category", "gene_recom_taig","recom_bin_taig","win_cons_pc")]
table_taiga_recomb<-na.omit(table_taiga_recomb)
table_taiga_recomb <- subset(table_taiga_recomb, chr != "ChrLGE22")
table_taiga_recomb <- subset(table_taiga_recomb, chr != "Chr25")

table_recomb<-merge(table_albicollis_recomb, table_taiga_recomb, by = "gene", all = T)

plot(table_taiga_recomb$gene_recom_taig, table_albicollis_recomb$gene_recom_coll)
cor(table_taiga_recomb$gene_recom_taig, table_albicollis_recomb$gene_recom_coll, method = "spearman")

#################################################################################################################################
ggplot(table_taiga_recomb, aes(x = new_category, y = gene_recom_taig, fill = new_category)) +
	geom_boxplot() +
	labs(x = "Category", y = "Recombination rate",
			 fill = "Category",  # Add this line to set the legend title
			 title = "") +
	theme(axis.text.y = element_text(angle = 90, vjust = 0.5))+
	theme(
		legend.title = element_text(size = 15),          
		legend.text = element_text(size = 14),          
		axis.title.x = element_text(size = 16),          
		axis.title.y = element_text(size = 16),        
		axis.text.x = element_text(size = 16),           
		axis.text.y = element_text(size = 12)
	)

anova<-lm(gene_recom_taig ~ new_category, data = table_taiga_recomb)
anova(anova)

ggplot(table_albicollis_recomb, aes(x = new_category, y = gene_recom_coll, fill = new_category )) +
  geom_boxplot() +
  labs(x = "Category", y = "Recombination rate",
       title = "Gene recombination rate for each category (F.albicollis)") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5))

table_taiga_recomb$chr <- reorder(table_taiga_recomb$chr, -table_taiga_recomb$length)
category_colors <- c("macrochromosome" = "green", "intermediate" = "blue", "microchromosome" = "red")
ggplot(table_taiga_recomb, aes(x = chr, y = gene_recom_taig, fill = new_category)) +
  geom_boxplot() +
  labs(x = "Chromosome", y = "Recombination rate",
       title = "Gene recombination rate for each chromosome (F.albicilla)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale

table_albicollis_recomb$chr <- reorder(table_albicollis_recomb$chr, -table_albicollis_recomb$length)
ggplot(table_albicollis_recomb, aes(x = chr, y = gene_recom_coll, fill =new_category)) +
  geom_boxplot() +
  labs(x = "Chromosome", y = "Recombination rate",
       title = "Gene recombination rate for each chromosome (F.albicollis)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale

########################################################################################################################################

table_taiga_recomb$recom_bin_taig <- factor(table_taiga_recomb$recom_bin_taig, levels = c("low", "med", "high"))
ggplot(table_taiga_recomb, aes(x = chr, y = gene_recom_taig, color = recom_bin_taig)) +
  geom_point() +
  labs(title = "Recombination for each chromosome category (F.albicilla)",
       x = "Chromosome category",
       y = "Recombination rate",
       color = "Recombination
    rate level") +
  scale_color_manual(values = c("low" = "green", "med" = "blue", "high" = "red"),
                     labels = c("low", "med", "high")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale


table_albicollis_recomb$recom_bin_coll <- factor(table_albicollis_recomb$recom_bin_coll, levels = c("low", "med", "high"))
ggplot(table_albicollis_recomb, aes(x = category, y = gene_recom_coll, color = recom_bin_coll)) +
  geom_point() +
  labs(title = "Recombination for each chromosome category (F.albicollis)",
       x = "Chromosome category",
       y = "Recombination rate",
       color = "Recombination
    rate level") +
  scale_color_manual(values = c("low" = "green", "med" = "blue", "high" = "red"),
                     labels = c("low", "med", "high")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Orientation verticale

#############################################################################################################################
nombre_high <- sum(table_taiga_recomb$recom_bin_taig == "high")
nombre_med <- sum(table_taiga_recomb$recom_bin_taig == "med")
nombre_low <- sum(table_taiga_recomb$recom_bin_taig == "low")
variance_by_group <- aggregate(gene_recom_taig ~ recom_bin_taig, data = table_taiga_recomb, FUN = var)
variance_by_chr <- aggregate(gene_recom_taig ~ chr, data = table_taiga_recomb, FUN = var)

nombre_high <- sum(table_albicollis_recomb$recom_bin_coll == "high")
nombre_med <- sum(table_albicollis_recomb$recom_bin_coll == "med")
nombre_low <- sum(table_albicollis_recomb$recom_bin_coll == "low")
variance_by_group <- aggregate(gene_recom_coll ~ recom_bin_coll, data = table_albicollis_recomb, FUN = var)
variance_by_chr <- aggregate(gene_recom_coll ~ chr, data = table_albicollis_recomb, FUN = var)
