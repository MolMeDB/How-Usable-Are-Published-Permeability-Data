# Install required packages
install.packages("ggpubr")
install.packages("dplyr")
install.packages("tidyverse")

# Load packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)

# Load MolMeDB dataset
data <- read.csv("inputs/prepared_MolMeDB_dataset.csv", header = TRUE, sep = ";")

config = rev(within(list(), {
  method_groups   = c("PerMM", "COSMOperm", "BLM", "CACO-2", "MDCK", "PAMPA")
  method_labels   = c("PerMM", "COSMOperm", "BLM/Liposomes","CACO-2", "MDCK", "PAMPA")
  # Plot colors
  colors          = list("PerMM" = 'goldenrod', 
                         "COSMOperm" = "forestgreen", 
                         "BLM" = "gray0", 
                         "PAMPA" = "deepskyblue", 
                         "CACO-2" = "deeppink3", 
                         "MDCK" = "darkmagenta")
}))


for(group_i in 1: length(config$method_groups))
{
  group = config$method_groups[[group_i]]
  label = config$method_labels[[group_i]]
  # load data
  g_data <- subset(data, method_label == group)
  
  molecules = g_data[!duplicated(g_data[ , c("SMILES")]), ]
 
  plotData <- ggplot(g_data, aes(x=LogPerm_molmedb_value)) +
           annotate("rect", xmin = -10, xmax = 10, ymin = 150, ymax = 1000, alpha = 0.2, fill = config$colors[[group]]) +
           annotate("text", x=0, y=800, label= paste("datsets =", length(unique(g_data$primary_reference))), size=8) +
           annotate("text", x=0, y=500, label= paste(c("n = "),sep= " ", collaps = length(unique(g_data$SMILES))), size=8) +
           annotate("text", x=0, y=300, label= paste(c("LogP = "),sep= " ", collaps = round(mean(molecules$LogP, na.rm = TRUE), digits = 1)), size=8) +
           annotate("text", x=0, y=175, label= paste(c("MW = "),sep= " ", collaps = round(mean(molecules$MW, na.rm=TRUE), digits = 1), paste(c("Da"))), size=8) +
           theme_bw()+
           geom_histogram(alpha = 0.2, fill= config$colors[[group]], color = config$colors[[group]], binwidth = 0.5) +
           labs(title = label, x = "LogPerm (cm/s)", y = "Occurences")+
           theme(text = element_text(size = 25)) +
           scale_y_log10(lim = c(1,1000)) +
           scale_x_continuous(lim = c(-20,10), breaks = c(-16, -12, -8, -4, 0, 4, 8, 12))+
           geom_vline(xintercept = c(-12), color = "gray31", size = 1)+
           geom_vline(xintercept = c(-8), color = "gray31", size = 0.5)+
           geom_vline(xintercept = c(-4), color = "gray31", size = 0.5, linetype = "longdash")+
           geom_vline(xintercept = c(0), color = "gray31", size = 0.5, linetype = "dashed")+
           geom_vline(xintercept = c(4), color = "gray31", size = 0.5, linetype = "dotdash")
  
  pdf(file=paste("outputs/histograms/", group, ".pdf", sep=""),
      width=16,
      height = 12)
  
  plotData
  
  dev.off()
}


