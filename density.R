
# Read data files
table_albicollis <- read.table("table_albicollis.txt", header = TRUE)
size <- read.table("chr_category.txt", header = TRUE)

# Initialize an empty DataFrame to store counts and densities
counts_df <- data.frame(chr = character(), count = numeric(), stringsAsFactors = FALSE)
densities_df <- data.frame(chr = character(), density = numeric(), stringsAsFactors = FALSE)

# Loop through unique chromosome values to count genes
for (chr_val in unique(table_albicollis$chr)) {
	count <- sum(grepl(paste0('^', chr_val, '$'), table_albicollis$chr))
	new_row <- data.frame(chr = chr_val, count = count)
	counts_df <- rbind(counts_df, new_row)
}

# Merge counts_df with size by 'chr'
colnames(size)[colnames(size) == "chromosome"] <- "chr"
merged_df <- merge(counts_df, size, by = "chr")

# Calculate density (count per length of chromosome)
merged_df$ratio <- merged_df$count / merged_df$length
merged_df$density <- merged_df$ratio * 100000

# Order and filter chromosomes
merged_df <- merged_df[mixedorder(merged_df$chr), ]
filtered_chrs <- c("ChrLGE22", "ChrFal34", "ChrFal36")
merged_df <- merged_df[!merged_df$chr %in% filtered_chrs, ]

# Plot gene density for each chromosome
ggplot(merged_df, aes(x = chr, y = ratio)) +
	geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	labs(x = "Chromosome", y = "Gene density", title = "Gene density for each chromosome (n = 9852 genes associated with chromosomes)") +
	theme_minimal()

library(ggplot2)
merged_df <- merged_df[mixedorder(merged_df$chr), ]
merged_df <- merged_df[merged_df$chr != "ChrLGE22", ]
merged_df <- merged_df[merged_df$chr != "ChrFal34", ]
merged_df <- merged_df[merged_df$chr != "ChrFal36", ]


merged_df$chr <- factor(merged_df$chr, levels = unique(merged_df$chr))



merged_df<-na.omit(merged_df)

color_palette <- c("Macro" = "tomato", "Micro" = "turquoise3")

ggplot(merged_df, aes(x = new_category, y = density, fill = new_category)) +
	geom_boxplot(position = position_dodge(width = 0.75)) +
	scale_fill_manual(values = color_palette) +
	labs(title = "",
			 x = "Category",
			 y = "Gene density",
			 fill = "Category") +
	theme(
		legend.title = element_text(size = 15),          
		legend.text = element_text(size = 14),          
		axis.title.x = element_text(size = 16),          
		axis.title.y = element_text(size = 16),        
		axis.text.x = element_text(size = 16),           
		axis.text.y = element_text(size = 12)
	)

anova<-lm(density ~ new_category, data = merged_df)
anova(anova)


density<-select(merged_df, chr, density)
