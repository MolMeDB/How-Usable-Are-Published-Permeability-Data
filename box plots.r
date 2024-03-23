# Install required packages
install.packages("ggpubr")
install.packages("dplyr")
install.packages("tidyverse")

# Load packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)

data_for_boxplot <- read.csv("inputs/box_plot_datset.csv", header = TRUE, sep = ",", blank.lines.skip = TRUE)

# Calculation of median
median_of_logPerm = median(data_for_boxplot$LogPerm)

box_plot_PAMPA <- ggplot(data_for_boxplot, aes(x = Name, y = LogPerm)) +
  geom_boxplot(color = "black", fill = "deepskyblue", alpha = 0.5, outlier.shape = NA)+
  theme_bw()+
  labs( y = "LogPerm (cm/s)", x = " ")+
  theme(text = element_text(size = 25), axis.text.x = element_text(angle = 90), legend.position = "none") +
  geom_hline(yintercept=median_of_logPerm) +
  scale_y_continuous(breaks = seq(-8, 2, by=1))+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  geom_hline(yintercept=median_of_logPerm)
  geom_point(data = subset(data_for_boxplot, LogPerm %in% boxplot.stats(data_for_boxplot$LogPerm)$out),
                                                aes(y = LogPerm, color = method),
                                                position = position_jitter(width = .2), size = 2)+
  scale_color_manual(values=c("dimgrey", "firebrick"))

# Save to file
pdf(file="outputs/box plots/01.pdf",
    width=16,
    height = 12)

box_plot_PAMPA

dev.off()
