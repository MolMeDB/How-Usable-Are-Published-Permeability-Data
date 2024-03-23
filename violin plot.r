# Instal required packages
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

##########################################################################
########################## Graphical abstract ############################
##########################################################################
########################## Prepare violin plot ###########################
##########################################################################
#### -- compares distribution of LogPerm values ##########################
#### -- for MDCK, CACO-2, PAMPA, BLM, COMSMOperm and PerMM methods -- ####
##########################################################################

config = rev(within(list(), {
  method_groups   = c("PerMM", "COSMOperm", "BLM", "PAMPA", "CACO-2", "MDCK")
  labels          = c("PerMM", "COSMOperm", "BLM / Liposomes", "PAMPA", "CACO-2", "MDCK")
  colors          = c("gray0", "deeppink3", "forestgreen", "darkmagenta", "deepskyblue", "goldenrod")
  yBreaks         = c(-16, -12, -8, -4, 0, 4, 8, 12)
}))

### Prepare data for violin plot
data_for_violinplot <- NULL

for (label in config$method_groups)
{
  data_for_violinplot = rbind(data_for_violinplot, subset(data, method_label == label))
}

# Create violinplot
graphical_abstract <- ggplot(data_for_violinplot, aes(x= factor(method_label, level = config$method_groups),y=LogPerm_molmedb_value, fill = method_label, color = method_label, alpha = 0.2))+
  scale_fill_manual(values=config$colors) +
  scale_color_manual(values=config$colors) +
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  labs( y = "LogPerm (cm/s)", x = " ") +
  scale_x_discrete(labels = config$labels) +
  scale_y_continuous(lim = c(-20,10), breaks = config$yBreaks) +
  theme(text = element_text(size = 20, color = "black"), ### Family font must be installed
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  theme(legend.position="none")+
  coord_flip()

# Print final result
pdf(file="outputs/violin plots/01.pdf",
    width=16,
    height = 8)

graphical_abstract

dev.off()

##########################################################################
########################## End of violin plot ############################
##########################################################################


