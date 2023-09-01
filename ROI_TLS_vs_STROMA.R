####################################################
########### ROI dataset analysis## #################
####################################################
library(visdat) # For visualizing data
library(ggplot2) # for nice plots
library(naniar) # For visualizing and working with missing data
library(GGally) # for combined plots
library(patchwork)
library(rstatix)
library(ggpubr)
library(EnhancedVolcano)

library(caret)
library(simputation) # For handling missing values
library(tidyverse) # Metapackage for dplyr and ggplot2
library("mice")



## Load dataset --------------------------------------------------------
# read.table is a function that allows to read multiple kinds of text files. It admits the following arguments:
df <- read.table(file = "STROMA_TLS_analysis.csv", #Name of text file.
                      sep = ";",                       #Separation character.
                      header = TRUE,                   #If column names are in the first row.
                      na.strings = "",               #Character to be marked as missing value.
                      stringsAsFactors = FALSE)         #żconvert string to factors

#--------------------------------------------------------------
# CREATE MAIN DATASETS ---------------------------------------------
# exclude tumor ROIs
df <- df[df$class != "INFILT",]

# exclude outliers (90 is identified as outlier)
df <- df[df$ROI != 90,]

# converting to factors
df$class = as.factor(df$class)
df$distance_to_tumor = as.factor(df$distance_to_tumor)
df$hist_Hans_B_WHO_2015 = as.factor(df$hist_Hans_B_WHO_2015)
# df$PDL1_status_tumor_1procent = as.factor(df$PDL1_status_tumor_1procent)
# df$PDL1_status_tumor_10procent = as.factor(df$PDL1_status_tumor_10procent)
# df$PDL1_status_tumor_50procent = as.factor(df$PDL1_status_tumor_50procent)
# df$PDL1_status_allcells_50procent = as.factor(df$PDL1_status_allcells_50procent)

#dataset containing only spatial metrics
spatial <- na.omit(df[,3:21])

#dataset used for PCA PC analysis of spatial metrics
spatialPCAset <- na.omit(spatial)

#dataset containing only probe metrics
barcode <- df[,-c(1:20, 22:23)]

#dataset used for PCA plot of all metrics
dfPCAAll <- na.omit(select(df, -Pt_core_scan_ROI, -distance_to_tumor, -hist_Hans_B_WHO_2015))

#dataset used for PCA plot of spatial metrics
dfPCASpatial <- na.omit(df[,2:21])

#dataset used for PCA plot of probe metrics
dfPCABarcode <- na.omit(df[,-c(1, 3:20)])
dfPCABarcode <- select(dfPCABarcode, -distance_to_tumor, -hist_Hans_B_WHO_2015)

#--------------------------------------------------------------------------
# DATA EXPLORATION ----------------------------------------------------
summary(df)
View(df)
head(df)

# graphical data exploration
ix_label <- df[df$ROI %in% c(32, 72, 73),]
ix_label <- ix_label[ix_label$class == "STROMA" & ix_label$CD20 >80,]

# Hide all of the text labels.
dfTemp <- df
dfTemp$ROI <- "" #ifelse((dfPCAResults$PC1 < 1 && dfPCAResults$class == "TLS") || (dfPCAResults$PC1 > 2.5 && dfPCAResults$class == "STROMA"), , "")
dfTemp <- rbind(ix_label, dfTemp)

library(ggrepel)
# Plot PCA results with class and ROI labels
ggplot(dfTemp, aes(x = CD20, y = Ki.67, color = class, label = ROI)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  labs(x = "CD20", y = "Ki-67", title = "CD20 vs Ki-67, investigation of outliers")

hist(df$S6)
ggplot()+geom_histogram(data = df, aes(x = S6, fill = class), alpha = 0.4,bins = 30, position="identity")
ggplot()+geom_histogram(data = df, aes(x = S6, fill = class),bins = 30)+facet_grid(class~.)
ggplot()+geom_histogram(data = df, aes(x = LAG3, fill = class),bins = 30)+facet_grid(class~.)
ggplot()+geom_histogram(data = df, aes(x = CD27, fill = class),bins = 30)+facet_grid(class~.)


#-------------------------------------------------------------------------------
# HANDLING OF NAs--------------------------------------------------------------
# Visualize the amount of missing data as well as the type of data
df %>% 
  vis_dat()

# To focus on the missing data only we can use
df %>% 
  vis_miss() 

# To check for interactions between missing values we do the following:
df %>% 
  gg_miss_upset()

# as a primary step, all RNA-sequence and PDL1 status data is excluded because of NAs and unknown relevance to analysis
# remove the columns with specific names
df <- df[, !names(df) %in% c("PDL1_status_tumor_1procent", "PDL1_status_tumor_10procent", "PDL1_status_tumor_50procent", "PDL1_status_allcells_50procent", "CD4RNAseq", "CD3RNAseq", "CD8RNAseq", "PD1RNAseq", "PDL1RNAseq", "CD44RNAseq")]
str(df)





#-------------------------------------------------------------------------------
# SPATIAL SIGNIFICANCE TESTS BETWEEN TLS AND STROMA & VOLCANO PLOTS -----------------------------
spatial <- na.omit(df[,3:21])
# list all features in the dataframe except "class"
features <- names(spatial)[-which(names(spatial) %in% c("class"))]

# create an empty dataframe to store the results
wilcoxon_spatial_volcano <- data.frame(feature = character(), fold_change = numeric(), p_adj = numeric())
ttest_volcano <- data.frame(feature = features)

# Function to format p-values with significance levels
format_p_value <- function(p_value) {
  if (p_value < 0.0000001) {
    return(" < 0.0000001")
  } else if (p_value < 0.000001) {
    return(" < 0.000001")
  } else if (p_value < 0.00001) {
    return(" < 0.00001")
  } else if (p_value < 0.0001) {
    return(" < 0.0001")
  } else if (p_value < 0.001) {
    return(" < 0.001")
  } else if (p_value < 0.01) {
    return(" < 0.01")
  } else if (p_value < 0.05) {
    return(" < 0.05")
  } else {
    return(paste(" = ", round(p_value, 3)))
  }
}

# loop over each feature and generate the plots
plots <- list()
for (f in features) {
  # Wilcoxon test
  res_wilcox <- wilcox.test(spatial[spatial$class == "TLS", f], spatial[spatial$class == "STROMA", f])
  
  # calculate fold-change
  fc <- log2(mean(spatial[spatial$class == "TLS", f])) - log2(mean(spatial[spatial$class == "STROMA", f]))
  ttest_volcano[ttest_volcano$feature == f, "fold_change"] <- fc
  
  # Get the p-value and adjust for multiple testing
  p_value <- res_wilcox$p.value
  p_adj <- p.adjust(p_value, method = "bonferroni", length(features))
  
  # Add the adjusted p-value to the results dataframe
  wilcoxon_spatial_volcano[nrow(wilcoxon_spatial_volcano) + 1,] <- list(feature = f, fold_change = fc, p_adj = p_adj)
  f_new <- gsub("_", " ", f)
  
  # Format p-values with significance levels
  formatted_p_value <- format_p_value(p_value)
  formatted_p_adj <- format_p_value(p_adj)
  
  # Plot
  plot <- ggboxplot(spatial, x = "class", y = f) +
    #geom_text(aes(label = paste0("p =", formatted_p_value)),
    #          x = 1, y = max(spatial[[f]]), vjust = -0.5, hjust = 1) +
    labs(
      y = gsub("\\.", "+", f_new),
      subtitle = paste("Wilcoxon test (Bonferroni-adjusted): p", formatted_p_adj),
      caption = ""  # remove the pwc label since there are no pairwise comparisons
    )
  
  # Add the plot to the list
  plots[[f]] <- plot
}

# Combine the plots using patchwork
wrap_plots <- wrap_plots(plots, ncol = 6)

# Display the plots
wrap_plots

# Wilcoxon VOLCANO PLOT
ggplot(wilcoxon_spatial_volcano, aes(x = fold_change, y = -log10(p_adj), label = gsub("_", " ", gsub("\\.", "+", feature)))) +
  geom_point(size = 2, alpha = 1) +
  geom_text_repel(size = 2.5, box.padding = 0.2, max.overlaps = Inf) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "Fold Change", y = "-log10(p_adjusted)", title = "Wilcoxon Volcano Plot")

# Wilcoxon VOLCANO PLOT
ggplot(wilcoxon_spatial_volcano, aes(x = fold_change, y = -log10(p_adj), label = gsub("_", " ", gsub("\\.", "+", feature)))) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#CC9999") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#CC9999") +
  geom_point(aes(size = -log10(p_adj), color = abs(fold_change)), alpha = 0.8) +
  geom_text_repel(size = 3, box.padding = 0.5, max.overlaps = Inf, nudge_y = 0) +
  labs(x = "Fold Change", y = "-log10(p-value)", title = "Wilcoxon Volcano Plot") +
  scale_color_gradient(low = "black", high = "#1F77B4") +
  scale_size(range = c(2, 4)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(30, 30, 30, 30, "pt")) +
  annotate("text", label = paste0("p_adjusted = 0.05"),
           x = -1.5, y = -log10(0.05), hjust = 0.2, vjust = -0.2, color = "#CC9999") +
  coord_cartesian(ylim = c(0.8, max(-log10(wilcoxon_spatial_volcano$p_adj), na.rm = TRUE) + 0.5))

# REGULAR VOLCANO PLOT
ggplot(ttest_volcano, aes(x = fold_change, y = -log10(p_value), label = gsub("_", " ", gsub("\\.", "+", feature)))) +
  geom_point(size = 2, alpha = 1) +
  geom_text_repel(size = 2.5, box.padding = 0.5, max.overlaps = Inf) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "Fold Change", y = "-log10(p-value)", title = "t-test Volcano Plot")
  theme(aspect.ratio=1)
  

#-------------------------------------------------------------------------------
# INDIVIDUAL TESTING (this deletes one ROI, imputation is recommended)------------------------------------------------------------------------------
res.kruskal <- spatial %>% kruskal_test(Combined_circularity ~ class)
res.kruskal
# Effect size  0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
spatial %>% kruskal_effsize(Combined_circularity ~ class)
# Pairwise comparisons
pwc <- spatial %>% 
  dunn_test(Combined_circularity ~ class, p.adjust.method = "bonferroni") 
pwc
pwc <- pwc %>% add_xy_position(x = "class")
ggboxplot(df, x = "class", y = "Combined_circularity") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

# combined spatial plot
ggpairs(spatial, aes(color =  class, alpha = 0.3))

#observing high correlations
# ALL VARIABLES HAD SOME SIGNIFICANCE BETWEEN CLASSES




#-------------------------------------------------------------------------------------------------------------------------------
# SPATIAL CORRELATION CHECK without NAs (this deletes one ROI, imputation is recommended)-----------------------------------------------------
# full correlation matrix
corrSpatial <- spatial[,1:18]
# corrSpatial <- df[,-c(1,2,
#                which(names(df) == "class"),
#                which(names(df) == "hist_Hans_B_WHO_2015"),
#                which(names(df) == "PDL1_status_tumor_1procent"),
#                which(names(df) == "PDL1_status_tumor_10procent"),
#                which(names(df) == "PDL1_status_tumor_50procent"),
#                which(names(df) == "PDL1_status_allcells_50procent"))]
# including ordinal distance_to_tumor as numerical to screen if it has any high correlations
#corrSpatial$distance_to_tumor = as.numeric(df$distance_to_tumor)
# remove NAs to do a sloppy but easy analysis
corrSpatial <- na.omit(corrSpatial)
corrSpatial.scaled <- scale(corrSpatial, center = TRUE, scale = TRUE)
cormat <- cor(corrSpatial.scaled)

library(reshape2)
# Reorder the correlation matrix
# reorder_cormat <- function(cormat){
#   # Use correlation between variables as distance
#   dd <- as.dist((1-cormat)/2)
#   hc <- hclust(dd)
#   cormat <-cormat[hc$order, hc$order]
# }
# cormat <- reorder_cormat(cormat)
cormat <- round(cormat, 2)
colnames(cormat) <- gsub("\\.", "+", gsub("_", " ", colnames(cormat)))
rownames(cormat) <- gsub("\\.", "+", gsub("_", " ", rownames(cormat)))

# Get the upper half of the matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#melted_cormat$value[melted_cormat$value < 0.5 & melted_cormat$value > -0.5] <- NA # for highly correlated
#melted_cormat$value[(melted_cormat$value < 0.3 & melted_cormat$value > -0.3) | melted_cormat$value > 0.5 | melted_cormat$value < -0.5] <- NA #highly correlated

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1))+
  coord_fixed()

print(ggheatmap)
# add the new style to the heatmap
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 20, barheight = 2,
                               title.position = "top", title.hjust = 0.5))



#-------------------------------------------------------------------------------------------
# BARCODE ANALYSIS BETWEEN TLS AND STROMA---------------------------------------------------
# find outliers that could potentially be mistakenly classified as stroma by checking CD20
# and mistakenly classified as TLS by checking fibronectin and SMA
barcode <- df[,-c(1:20, 22:23)]
# list all features in the dataframe except "class"
features <- names(barcode)[-which(names(barcode) %in% c("class"))]

# create an empty dataframe to store the results
wilcoxon_barcode_volcano <- data.frame(feature = character(), fold_change = numeric(), p_adj = numeric())
ttest_barcode_volcano <- data.frame(feature = features)

# loop over each feature and generate the plots
plots <- list()
for (f in features) {
  # Wilcoxon test
  res_wilcox <- wilcox.test(barcode[barcode$class == "TLS", f], barcode[barcode$class == "STROMA", f])
  
  # calculate fold-change
  fc <- log2(mean(barcode[barcode$class == "TLS", f])) - log2(mean(barcode[barcode$class == "STROMA", f]))
  ttest_barcode_volcano[ttest_barcode_volcano$feature == f, "fold_change"] <- fc
  
  # Get the p-value and adjust for multiple testing
  p_value <- res_wilcox$p.value
  p_adj <- p.adjust(p_value, method = "bonferroni", length(features))
  
  # Add the adjusted p-value to the results dataframe
  wilcoxon_barcode_volcano[nrow(wilcoxon_barcode_volcano) + 1,] <- list(feature = f, fold_change = fc, p_adj = p_adj)
  f_new <- gsub("_", " ", f)
  
  # Format p-values with significance levels
  formatted_p_value <- format_p_value(p_value)
  formatted_p_adj <- format_p_value(p_adj)
  
  # Plot
  plot <- ggboxplot(barcode, x = "class", y = f) +
    #geom_text(aes(label = paste0("p =", formatted_p_value)),
    #          x = 1, y = max(barcode[[f]]), vjust = -0.5, hjust = 1) +
    labs(
      y = gsub("\\.", "-", f_new),
      subtitle = paste("Wilcoxon test (Bonferroni-adjusted): p", formatted_p_adj),
      caption = ""  # remove the pwc label since there are no pairwise comparisons
    )
  
  # Add the plot to the list
  plots[[f]] <- plot
}
# Combine the plots using patchwork
wrap_plots <- wrap_plots(plots, ncol = 11)

# Display the plots
wrap_plots






#-------------------------------------------------------------------------------------------
# INDIVIDUAL TESTING------------------------------------------------------------------------------
res.kruskal <- barcode %>% kruskal_test(Fibronectin ~ class)
res.kruskal
# Effect size  0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
barcode %>% kruskal_effsize(Fibronectin ~ class)
# Pairwise comparisons
pwc <- barcode %>% 
  dunn_test(Fibronectin ~ class, p.adjust.method = "bonferroni") 
pwc
pwc <- pwc %>% add_xy_position(x = "class")
ggboxplot(df, x = "class", y = "Fibronectin") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )



#--------------------------------------------------------------------------------------------------------------------------------
# PCA ALL MEASUREMENTS---------------------------
#Use hirarchical clustering or PCA with around 10-20 features.
#We can do a substudy and see if patients that have a high ROI count behave 
# similar to those that are sampled only a few times. Not necessary though.
dfPCAAll <- na.omit(select(df, -Pt_core_scan_ROI, -distance_to_tumor, -hist_Hans_B_WHO_2015))

# Perform PCA
pcaAll <- prcomp(select(dfPCAAll, -class, -ROI), center = TRUE, scale. = TRUE)

# Add class and ROI columns back to PCA results
dfPCAAllResults <- data.frame(
  PC1 = pcaAll$x[, 1],
  PC2 = pcaAll$x[, 2],
  class = dfPCAAll$class,
  ROI = dfPCAAll$ROI
)

# USED IF WANTING TO FILTER ON AXIS VALUES
# # Let's just label the outliers.
# ix_label <- dfPCAAllResults[dfPCAAllResults$PC2 < -1,]
# ix_label <- ix_label[ix_label$class == "TLS",]
# ix_label2 <- dfPCAAllResults[dfPCAAllResults$PC1 < -3,]
# ix_label2 <- ix_label2[ix_label2$class == "STROMA",]
# ix_label <- rbind(ix_label, ix_label2)
# # Hide all of the text labels.
# dfPCAAllResults$ROI <- "" #ifelse((dfPCAAllResults$PC1 < 1 && dfPCAAllResults$class == "TLS") || (dfPCAAllResults$PC1 > 2.5 && dfPCAAllResults$class == "STROMA"), , "")
# dfPCAAllResults <- rbind(ix_label, dfPCAAllResults)


#USED WHEN WANTING TO SELECT SPECIFIC ROIs
# Retain specific "ROI" values for class "TLS"
tls_roi <- c(1, 2, 3, 31, 49, 50, 55)
dfTLS <- dfPCAAllResults[dfPCAAllResults$class == "TLS" & dfPCAAllResults$ROI %in% tls_roi, ]
dfTLS <- rbind(dfTLS,dfPCAAllResults[dfPCAAllResults$class == "TLS" & dfPCAAllResults$ROI == 48 & dfPCAAllResults$PC1 < -5,])

# Retain specific "ROI" values for class "STROMA"
stroma_roi <- c(32, 73)
dfSTROMA <- dfPCAAllResults[dfPCAAllResults$class == "STROMA" & dfPCAAllResults$ROI %in% stroma_roi, ]

# Bind the filtered dataframes together
dfPCAAllResults$ROI <- "" #ifelse((dfPCAAllResults$PC1 < 1 && dfPCAAllResults$class == "TLS") || (dfPCAAllResults$PC1 > 2.5 && dfPCAAllResults$class == "STROMA"), , "")
dfPCAAllResults <- rbind(dfTLS, dfSTROMA, dfPCAAllResults)


library(ggrepel)
# Plot PCA results with class and ROI labels
ggplot(dfPCAAllResults, aes(x = PC1, y = PC2, color = class, label = ROI)) +
  geom_point() +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  labs(x = "PC1", y = "PC2", title = "PCA plot with class and ROI labels, all measurements")

# Möjliga felannoteringar: 
# TLS: 1, 2, 3 i scan2
# STROMA: 31, 73 scan 1 & 32 scan2


dfPCAAll <- na.omit(select(df, -Pt_core_scan_ROI, -distance_to_tumor, -hist_Hans_B_WHO_2015))
corrAll <- na.omit(select(dfPCAAll, -class, -ROI))
corrAll.scaled <- scale(corrAll, center = TRUE, scale = TRUE)
cormatAll <- cor(corrAll.scaled)
principalComponentsAll <- princomp(cormatAll)
summary(principalComponentsAll)
principalComponentsAll$loadings[, 1:2]
# Create a data frame with loadings and row names in the original order
loadings_all <- data.frame(PC1 = principalComponentsAll$loadings[, 1],
                       PC2 = principalComponentsAll$loadings[, 2])
loadings_all$row_names <- rownames(principalComponentsAll$loadings)
loadings_all$row_names[1:20] <- gsub("_", " ", gsub("\\.", "+", loadings_all$row_names[1:20]))
loadings_all$row_names[-c(1:20)] <- gsub("_", " ", gsub("\\.", "-", loadings_all$row_names[-c(1:20)]))

# Create a new column with row names as factors
loadings_all$row_names <- factor(loadings_all$row_names, levels = loadings_all$row_names)

# Sort the data frame based on the new factor column
loadings_all <- loadings_all[order(loadings_all$row_names), ]

# Create a histogram using ggplot with the original order of row names
loadings_all_plot <- ggplot(loadings_all, aes(x = row_names)) +
  geom_bar(aes(y = PC1, fill = "PC1"), stat = "identity") +
  geom_bar(aes(y = PC2, fill = "PC2"), stat = "identity") +
  labs(x = "Variable names", y = "Loadings", fill = "Principal Component") +
  ggtitle("Principal Component Loadings") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("blue", "red"),
                    labels = c("PC1", "PC2"),
                    guide = guide_legend(title = NULL))


# Print the histogram
print(loadings_all_plot)
library(factoextra)
fviz_eig(principalComponentsAll, addlabels = TRUE, geom="bar", hjust=0.4, xlab = "Principal Components", ylab = "Proportion of Variance (%)", title = "Proportion of Variance in PCs")


#--------------------------------------------------------------------------------------------------------------------------------
# PCA SPATIALS---------------------------------------------------------------
#Use hirarchical clustering or PCA with around 10-20 features.
#We can do a substudy and see if patients that have a high ROI count behave 
# similar to those that are sampled only a few times. Not necessary though.
dfPCASpatial <- na.omit(df[,2:21])

# Perform PCA
pcaSpatial <- prcomp(select(dfPCASpatial, -class, -ROI), center = TRUE, scale. = TRUE)

# Add class and ROI columns back to PCA results
dfPCASpatialResults <- data.frame(
  PC1 = pcaSpatial$x[, 1],
  PC2 = pcaSpatial$x[, 2],
  class = dfPCASpatial$class,
  ROI = dfPCASpatial$ROI
)

# Let's just label the outliers.
ix_label <- dfPCASpatialResults[dfPCASpatialResults$ROI %in% c(1, 2, 3),]
ix_label <- ix_label[ix_label$class == "TLS",]
ix_label2 <- dfPCASpatialResults[dfPCASpatialResults$ROI %in% c(32, 72, 73) & dfPCASpatialResults$PC1 < 0,]
ix_label2 <- ix_label2[ix_label2$class == "STROMA",]
ix_label <- rbind(ix_label, ix_label2)
# Hide all of the text labels.
dfPCASpatialResults$ROI <- "" #ifelse((dfPCASpatialResults$PC1 < 1 && dfPCASpatialResults$class == "TLS") || (dfPCASpatialResults$PC1 > 2.5 && dfPCASpatialResults$class == "STROMA"), , "")
dfPCASpatialResults <- rbind(ix_label, dfPCASpatialResults)

library(ggrepel)
# Plot PCA results with class and ROI labels
ggplot(dfPCASpatialResults, aes(x = PC1, y = PC2, color = class, label = ROI)) +
  geom_point() +
  geom_text_repel(box.padding = 0.7, max.overlaps = Inf) +
  labs(x = "PC1", y = "PC2", title = "PCA plot with class and ROI labels, spatial features")


#principal component analysis of all spatial PCs
dfPCASpatial <- na.omit(df[,2:21])
corrSpatial <- na.omit(select(dfPCASpatial, -class, -ROI))
corrSpatial.scaled <- scale(corrSpatial, center = TRUE, scale = TRUE)
cormatSpatial <- cor(corrSpatial.scaled)
principalComponentsSpatial <- princomp(cormatSpatial)
summary(principalComponentsSpatial)
principalComponentsSpatial$loadings[, 1:2]
library(factoextra)
POVspatialPlot <- fviz_eig(principalComponentsSpatial, addlabels = TRUE, geom="bar", hjust=0.4, xlab = "Principal Components", ylab = "Proportion of Variance (%)", title = "Proportion of Variance in Spatial PCs")
POVspatialPlot

# Create a data frame with loadings and row names in the original order
loadings_spatial <- data.frame(PC1 = principalComponentsSpatial$loadings[, 1],
                           PC2 = principalComponentsSpatial$loadings[, 2])
loadings_spatial$row_names <- rownames(principalComponentsSpatial$loadings)
loadings_spatial$row_names <- gsub("_", " ", gsub("\\.", "+", loadings_spatial$row_names))

# Create a new column with row names as factors
loadings_spatial$row_names <- factor(loadings_spatial$row_names, levels = loadings_spatial$row_names)

# Sort the data frame based on the new factor column
loadings_spatial <- loadings_spatial[order(loadings_spatial$row_names), ]

# Create a histogram using ggplot with the original order of row names
loadings_spatial_plot <- ggplot(loadings_spatial, aes(x = row_names)) +
  geom_bar(aes(y = PC1, fill = "PC1"), stat = "identity") +
  geom_bar(aes(y = PC2, fill = "PC2"), stat = "identity") +
  labs(x = "Variable names", y = "Loadings", fill = "Principal Component") +
  ggtitle("Principal Component Loadings") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("blue", "red"),
                    labels = c("PC1", "PC2"),
                    guide = guide_legend(title = NULL))


# Print the histogram
print(loadings_spatial_plot)

#only for TLS
dfPCASpatial <- na.omit(df[,2:21])
corrSpatialTLS <- dfPCASpatial[dfPCASpatial$class == "TLS",]
corrSpatialTLS <- select(corrSpatialTLS, -class, -ROI, -Number_of_islands)
corrSpatialTLS.scaled <- scale(corrSpatialTLS, center = TRUE, scale = TRUE)
cormatSpatialTLS <- cor(corrSpatialTLS.scaled)
principalComponentsSpatialTLS <- princomp(cormatSpatialTLS)
summary(principalComponentsSpatialTLS)
principalComponentsSpatialTLS$loadings[, 1:2]
library(factoextra)
POVspatialTLSPlot <- fviz_eig(principalComponentsSpatialTLS, addlabels = TRUE, geom="bar", hjust=0.4, xlab = "Principal Components", ylab = "Proportion of Variance (%)", title = "Proportion of Variance in TLS Spatial PCs")
POVspatialTLSPlot

#only for stroma
dfPCASpatial <- na.omit(df[,2:21])
corrSpatialSTROMA <- dfPCASpatial[dfPCASpatial$class == "STROMA",]
corrSpatialSTROMA <- select(corrSpatialSTROMA, -class, -ROI)
corrSpatialSTROMA.scaled <- scale(corrSpatialSTROMA, center = TRUE, scale = TRUE)
cormatSpatialSTROMA <- cor(corrSpatialSTROMA.scaled)
principalComponentsSpatialSTROMA <- princomp(cormatSpatialSTROMA)
summary(principalComponentsSpatialSTROMA)
principalComponentsSpatialSTROMA$loadings[, 1:2]
library(factoextra)
POVspatialSTROMAPLot <- fviz_eig(principalComponentsSpatialSTROMA, addlabels = TRUE, geom="bar", hjust=0.4, xlab = "Principal Components", ylab = "Proportion of Variance (%)", title = "Proportion of Variance in stromal Spatial PCs")
POVspatialSTROMAPLot

par(mfrow=c(1,3))
plot(POVspatialPlot, main = "Principal components", ylim = c(0, 9))
plot(POVspatialTLSPlot, main = "TLS principal components", ylim = c(0, 9))
plot(POVspatialSTROMAPLot, main = "STROMA principal components", ylim = c(0, 9))


#--------------------------------------------------------------------------------------------------------------------------------
# PCA BARCODE------------------------------------------------------------
#Use hirarchical clustering or PCA with around 10-20 features.
#We can do a substudy and see if patients that have a high ROI count behave 
# similar to those that are sampled only a few times. Not necessary though.
dfPCABarcode <- na.omit(df[,-c(1, 3:20)])
dfPCABarcode <- select(dfPCABarcode, -distance_to_tumor, -hist_Hans_B_WHO_2015)

# Perform PCA
pcaBarcode <- prcomp(select(dfPCABarcode, -class, -ROI), center = TRUE, scale. = TRUE)

# Add class and ROI columns back to PCA results
dfPCABarcodeResults <- data.frame(
  PC1 = pcaBarcode$x[, 1],
  PC2 = pcaBarcode$x[, 2],
  class = dfPCABarcode$class,
  ROI = dfPCABarcode$ROI
)

# Let's just label the outliers.
ix_label <- dfPCABarcodeResults[dfPCABarcodeResults$ROI %in% c(1, 2, 3),]
ix_label <- ix_label[ix_label$class == "TLS",]
ix_label2 <- dfPCABarcodeResults[dfPCABarcodeResults$ROI %in% c(32, 72, 73) & dfPCABarcodeResults$PC1 > 0,]
ix_label2 <- ix_label2[ix_label2$class == "STROMA",]
ix_label <- rbind(ix_label, ix_label2)
# Hide all of the text labels.
dfPCABarcodeResults$ROI <- "" #ifelse((dfPCABarcodeResults$PC1 < 1 && dfPCABarcodeResults$class == "TLS") || (dfPCABarcodeResults$PC1 > 2.5 && dfPCABarcodeResults$class == "STROMA"), , "")
dfPCABarcodeResults <- rbind(ix_label, dfPCABarcodeResults)

library(ggrepel)
# Plot PCA results with class and ROI labels
ggplot(dfPCABarcodeResults, aes(x = PC1, y = PC2, color = class, label = ROI)) +
  geom_point() +
  geom_text_repel(box.padding = 0.65, max.overlaps = Inf) +
  labs(x = "PC1", y = "PC2", title = "PCA plot with class and ROI labels, probe features")



dfPCABarcode <- na.omit(df[,-c(1, 3:20)])
dfPCABarcode <- select(dfPCABarcode, -distance_to_tumor, -hist_Hans_B_WHO_2015)
corrBarcode <- na.omit(select(dfPCABarcode, -class, -ROI))
corrBarcode.scaled <- scale(corrBarcode, center = TRUE, scale = TRUE)
cormatBarcode <- cor(corrBarcode.scaled)
principalComponentsBarcode <- princomp(cormatBarcode)
summary(principalComponentsBarcode)
principalComponentsBarcode$loadings[, 1:2]
library(factoextra)
fviz_eig(principalComponentsBarcode, addlabels = TRUE, geom="bar", hjust=0.4, xlab = "Principal Components", ylab = "Proportion of Variance (%)", title = "Proportion of Variance in Probe data PCs")



#-------------------------------------------------------------------------
# TLS marker Outliers & Sanity checks------------------------------------------------------------------------

#CD20: Höga stroma och låga TLS kan vara outliers (felklassificerade)
ggplot()+geom_histogram(data = df, aes(x = log(CD20), fill = class),bins = 30)+facet_grid(class~.)
#en hög i stroma + två ganska höga, 2 ganska låga i TLS -> kolla uttrycket i Ki67
#Ki67: Höga stroma och låga TLS kan vara outliers
ggplot()+geom_histogram(data = df, aes(x = log(Ki.67), fill = class),bins = 30)+facet_grid(class~.)
#en outlier i stroma
#jämför båda
df[df$class == "TLS",] %>% 
  ggplot(aes(x = log(CD20), y = log(Ki.67), color=class)) +
  geom_miss_point()
# de TLS som är låga i CD20 är också låga i Ki67
df[df$class == "STROMA",] %>% 
  ggplot(aes(x = log(CD20), y = log(Ki.67), color=class)) +
  geom_miss_point()
# vi ser att det är samma outlier, de potentiella stroma-outliers har också högt uttryck av Ki67 men inte mycket mer än de andra
df %>% 
  ggplot(aes(x = log(CD20), y = log(Ki.67), color=class)) +
  geom_miss_point()
# de två misstänkta stroma- och de två misstänkta TLS-outliers passar bra in i varandras fördelningar
# ROI 90 är outlier exkluderas.56 och 73 scan1 är potentiellt TLS. 3 scan2 & 62 scan1 är potentiellt STROMA
df <- df[df$CD20 < 1000,]

#CD45: Höga stroma och låga TLS kan vara outliers (felklassificerade)
ggplot()+geom_histogram(data = df, aes(x = log(CD45), fill = class),bins = 30)+facet_grid(class~.)
#potentiellt en hög i stroma och 3 låga i TLS.

df[df$class == "TLS",] %>% 
  ggplot(aes(x = log(CD20), y = log(CD45), color=class)) +
  geom_miss_point()
# en av dessa är låga även i CD20, de andra två hyfsat lågt men närmare mitten
df[df$class == "TLS",] %>% 
  ggplot(aes(x = log(Ki.67), y = log(CD45), color=class)) +
  geom_miss_point()
# två av dessa är låga även i Ki67

df[df$class == "STROMA",] %>% 
  ggplot(aes(x = log(CD20), y = log(CD45), color=class)) +
  geom_miss_point()
df[df$class == "STROMA",] %>% 
  ggplot(aes(x = log(Ki.67), y = log(CD45), color=class)) +
  geom_miss_point()
# den höga CD45-punkten ligger nära mitten av de andra

df %>% 
  ggplot(aes(x = log(CD20), y = log(CD45), color=class)) +
  geom_miss_point()
# den låga TLS-punkten passar in i STROMA-fördelningen. Denna är ROI 3 scan2
df %>% 
  ggplot(aes(x = log(Ki.67), y = log(CD45), color=class)) +
  geom_miss_point()
# de två låga punkterna passar in i STROMA-fördelningen
# ROI 3 scan2 samt 62 scan1

#Cell density: Höga stroma och låga TLS kan vara outliers
ggplot()+geom_histogram(data = df, aes(x = log(Cell.density.mum..1), fill = class),bins = 30)+facet_grid(class~.)
# potentiellt en låg TLS
df[df$class == "TLS",] %>% 
  ggplot(aes(x = log(CD20), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
# även hyfsat låg i CD20
df[df$class == "TLS",] %>% 
  ggplot(aes(x = log(Ki.67), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
# typ i mitten i Ki67
df[df$class == "TLS",] %>% 
  ggplot(aes(x = log(CD45), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
# hög i CD45

df[df$class == "STROMA",] %>% 
  ggplot(aes(x = log(CD20), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
df[df$class == "STROMA",] %>% 
  ggplot(aes(x = log(Ki.67), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
df[df$class == "STROMA",] %>% 
  ggplot(aes(x = log(CD45), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
# no indication of outlier
df %>% 
  ggplot(aes(x = log(CD20), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
df %>% 
  ggplot(aes(x = log(Ki.67), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
df %>% 
  ggplot(aes(x = log(CD45), y = log(Cell.density.mum..1), color=class)) +
  geom_miss_point()
# passar hyfsat bra in i stromafördelningen men bedöms inte vara outlier.

ggplot()+geom_histogram(data = df, aes(x = centrality_measures_CD45.Ignore, fill = class),bins = 30)+facet_grid(class~.)
#3scan2 is low
ggplot()+geom_histogram(data = df, aes(x = CD45.ratio, fill = class),bins = 30)+facet_grid(class~.)
#31scan2 is low
ggplot()+geom_histogram(data = df, aes(x = Combined_mean_intensity, fill = class),bins = 30)+facet_grid(class~.)
#31scan2, 3scan2 are low TLS


# ROI 90 scan1 is way overexpressed in many metrics, including CD20 -> outlier

# potentally misclassified as STROMA:
# ROI 56 scan1: CD20+, Ki67+, CD45+
# ROI 73 scan1: CD20+, Ki67+, CD45+

# potentially misclassified as TLS:
# ROI 3 scan2: CD20-, Ki67-, CD45-, GDC-, intensity-
# ROI 62 scan1: CD20-, Ki67-, CD45-
# ROI 31 scan2: CD45ratio-, intensity-
# (ROI 48 scan2: SMA+ TLS)


#Sanity checks
df[df$class == "TLS",] %>% 
  ggplot(aes(x = log(CD45._ratio), y = log(Mean_intensity), color=class)) +
  geom_miss_point() +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  labs(y="log(Mean CD45 intensity)", x="log(CD45+ ratio)") +
  annotate("text", label = paste0("R^2 = ", round(summary(lm(log(CD45) ~ log(Mean_intensity), data = df[df$class == "TLS",]))$r.squared, 3)),
           x = -Inf, y = Inf, hjust = 0, vjust = 1.5) +
  annotate("text", label = paste0("p-value = ", round(summary(lm(log(CD45) ~ log(Mean_intensity), data = df[df$class == "TLS",]))$coefficients[8], 5)),
           x = -Inf, y = Inf, hjust = 0, vjust = 3) +
  scale_color_manual(values = "#619CFF")

df[df$class == "STROMA",] %>% 
  ggplot(aes(x = log(CD45._ratio), y = log(Mean_intensity), color=class)) +
  geom_miss_point() +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  labs(y="log(Mean CD45 intensity)", x="log(CD45+ ratio)") +
  annotate("text", label = paste0("R^2 = ", round(summary(lm(log(CD45) ~ log(Mean_intensity), data = df[df$class == "STROMA",]))$r.squared, 3)),
           x = -Inf, y = Inf, hjust = 0, vjust = 1.5) +
  annotate("text", label = paste0("p-value = ", round(summary(lm(log(CD45) ~ log(Mean_intensity), data = df[df$class == "STROMA",]))$coefficients[8], 5)),
           x = -Inf, y = Inf, hjust = 0, vjust = 3)
  scale_color_manual(values = "#F8766D")
df %>% 
  ggplot(aes(x = log(CD45._ratio), y = log(Mean_intensity), color=class)) +
  geom_miss_point() +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth") +
  labs(y="log(Mean CD45 intensity)", x="log(CD45+ ratio)") +
  scale_color_manual(values = c("#F8766D", "#619CFF"))




#-------------------------------------------------------------------------------------------------------------------------------
# BARCODE CORRELATION CHECK without NAs (this deletes one ROI, imputation is recommended------------------------------------------------
corr <- barcode[,-c(1)]
# remove NAs to do a sloppy but easy analysis
corr <- na.omit(corr)
corr.scaled <- scale(corr, center = TRUE, scale = TRUE)
cormat <- cor(corr.scaled)
library(reshape2)
# Reorder the correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
cormat <- round(reorder_cormat(cormat), 2)
# Get the upper half of the matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
print(ggheatmap)
# add the new style to the heatmap
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 17, barheight = 2,
                               title.position = "top", title.hjust = 0.5))

#Outliers from high correlation
df %>% 
  ggplot(aes(x = CD27, y = S6, color=class)) +
  geom_miss_point()
