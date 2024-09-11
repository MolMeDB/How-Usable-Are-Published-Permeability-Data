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
 
  if(group != "PAMPA")
  {
    plotData <- ggplot(g_data, aes(x=LogPerm_molmedb_value)) +
           annotate("rect", xmin = -11.2, xmax = 10, ymin = 950, ymax = 9000, alpha = 0.2, fill = config$colors[[group]]) +
           annotate("text", x=0, y=6400, label= paste("datasets =", length(unique(g_data$primary_reference))), size=8) +
           annotate("text", x=0, y=3800, label= paste(c("n = "),sep= " ", collaps = length(unique(g_data$SMILES))), size=8) +
           annotate("text", x=0, y=2200, label= paste(c("LogP = "),sep= " ", collaps = round(mean(molecules$LogP, na.rm = TRUE), digits = 1)), size=8) +
           annotate("text", x=0, y=1200, label= paste(c("MW = "),sep= " ", collaps = round(mean(molecules$MW, na.rm=TRUE), digits = 1), paste(c("Da"))), size=8) +
           theme_bw()+
           labs(title = label, x = "LogPerm (cm/s)", y = "Occurences")+
           geom_histogram(alpha = 0.2, fill= config$colors[[group]], color = config$colors[[group]], binwidth = 0.5) +
           scale_y_log10(lim = c(1,10000)) +
           scale_x_continuous(lim = c(-20,10), breaks = c(-16, -12, -8, -4, 0, 4, 8, 12))+
           geom_vline(xintercept = c(-12), color = "gray31", size = 1)+
           geom_vline(xintercept = c(-8), color = "gray31", size = 0.5)+
           geom_vline(xintercept = c(-4), color = "gray31", size = 0.5, linetype = "longdash")+
           geom_vline(xintercept = c(0), color = "gray31", size = 0.5, linetype = "dashed")+
           geom_vline(xintercept = c(4), color = "gray31", size = 0.5, linetype = "dotdash")
  }
  else
  {
    plotData<-ggplot(g_data, aes(x=LogPerm_molmedb_value, fill = method, color = method)) +
          theme_bw()+
          geom_histogram(alpha = 0.2, binwidth = 0.5) +
          labs(title = "PAMPA (<span style = 'color:lightskyblue3;'>app.</span>) / (<span style = 'color:navy;'>intr.</span>)")+
          labs(x = "LogPerm (cm/s)", y = "Occurences")+
          theme(text = element_text(size = 25)) +
          scale_y_log10(lim = c(1,10000)) +
          scale_color_manual(values=c("lightskyblue3", "navy"))+
          scale_fill_manual(values=c("lightskyblue", "navy"))+
          scale_x_continuous(lim = c(-20,10), breaks = c(-16, -12, -8, -4, 0, 4, 8, 12))+
          geom_vline(xintercept = c(-12), color = "gray31", size = 1)+
          geom_vline(xintercept = c(-8), color = "gray31", size = 0.5)+
          geom_vline(xintercept = c(-4), color = "gray31", size = 0.5, linetype = "longdash")+
          geom_vline(xintercept = c(0), color = "gray31", size = 0.5, linetype = "dashed")+
          geom_vline(xintercept = c(4), color = "gray31", size = 0.5, linetype = "dotdash") +
          theme(legend.position="none") +
          theme(plot.title = element_markdown()) +
          annotate("rect", xmin = -20, xmax = -7, ymin = 920, ymax = 8500, alpha = 0.2, fill = "deepskyblue")+
          annotate("text", x= -13.5, y=6400, label= paste(c("datasets = "),sep= " ", collaps = length(unique(g_data$primary_reference))), size=7.5)+
          annotate("text", x= -13.5, y=3800, label= paste(c("n = "),sep= " ", collaps = length(unique(g_data$SMILES))), size=7.5)+
          annotate("text", x= -13.5, y=2200, label= paste(c("LogP = "),sep= " ", collaps = round(mean(molecules$LogP, na.rm = TRUE))), size=7.5) +
          annotate("text", x= -13.5, y=1200, label= paste(c("MW = "),sep= " ", collaps = round(mean(molecules$MW, na.rm=TRUE), digits = 1), paste(c("Da"))), size=7.5)
  }
  pdf(file=paste("outputs/histograms/", group, ".pdf", sep=""),
      width=16,
      height = 12)
  
  plotData
  
  dev.off()
}


