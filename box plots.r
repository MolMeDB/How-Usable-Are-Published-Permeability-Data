# Install required packages
install.packages("ggpubr")
install.packages("dplyr")
install.packages("tidyverse")

# Load packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)

data_for_boxplot <- read.csv("inputs/box_plot_datset.csv", header = TRUE, sep = ";", blank.lines.skip = TRUE)

EPAM <- subset(data_for_boxplot, method == "EPAM")

# Calculation of median
EPAM_median_of_logPerm = median(EPAM$LogPerm)

box_plot_PAMPA_EPAM <- ggplot(EPAM, aes(x = Name, y = LogPerm)) +
  geom_boxplot(color = "black", fill = "lightskyblue", alpha = 0.5, outlier.shape = NA)+
  theme_bw()+
  labs( y = "LogPerm (cm/s)", x = " ")+
  theme(text = element_text(size = 25), axis.text.x = element_text(angle = 90), legend.position = "none") +
  stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  labs(title = "A")+
  geom_hline(yintercept=EPAM_median_of_logPerm)

box_plot_PAMPA_EPAM  <- box_plot_PAMPA_EPAM  + geom_point(data = subset(EPAM, LogPerm %in% boxplot.stats(EPAM$LogPerm)$out),
                    aes(y = LogPerm, color = method),
                    position = position_jitter(width = .2), size = 2)+
                    scale_color_manual(values=c("dimgray", "firebrick"))+
                    scale_y_continuous(breaks = seq(-8,3, by=1), limits = c(-8,0))

# Save to file
pdf(file="outputs/box plots/01-A.pdf",
    width=16,
    height = 12)

box_plot_PAMPA_EPAM

dev.off()


## The same with EPAMOL ##
EPAMOL <- subset(data_for_boxplot, method == "EPAMOL")

EMAPOL_median_of_logPerm = median(EPAMOL$LogPerm)

box_plot_PAMPA_EPAMOL <- ggplot(EPAMOL, aes(x = Name, y = LogPerm)) +
  geom_boxplot(color = "black", fill = "navy", alpha = 0.5, outlier.shape = NA)+
  theme_bw()+
  labs( y = "LogPerm (cm/s)", x = " ")+
  theme(text = element_text(size = 25), axis.text.x = element_text(angle = 90), legend.position = "none") +
  stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  labs(title = "B")+
  geom_hline(yintercept=EMAPOL_median_of_logPerm)

box_plot_PAMPA_EPAMOL  <- box_plot_PAMPA_EPAMOL  + geom_point(data = subset(EPAMOL, LogPerm %in% boxplot.stats(EPAMOL$LogPerm)$out),
                    aes(y = LogPerm, color = method),
                    position = position_jitter(width = .2), size = 2)+
                    scale_color_manual(values=c("dimgrey", "firebrick"))+
                    scale_y_continuous(breaks = seq(-8,3, by=1), limits = c(-8,0))


# Save to file
pdf(file="outputs/box plots/01-B.pdf",
    width=16,
    height = 12)

box_plot_PAMPA_EPAMOL

dev.off()
