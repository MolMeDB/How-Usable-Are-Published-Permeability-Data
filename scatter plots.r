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

config = rev(within(list(), {
  method_groups   = c("PerMM", "COSMOperm", "BLM", "CACO-2", "MDCK", "EPAM", "EPAMOL")
  method_labels   = c("PerMM", "COSMOperm", "BLM/Liposomes","CACO-2", "MDCK", "PAMPA (app.)", "PAMPA (intr.)")
  groups_col      = c(rep("1", 5), rep("2", 2))
  datasets_keys   = list("BLM" = c("PERMM_COSMO", "PerMM", "MDCK"), 
                      "PAMPA" = c("MDCK_CACO", "BLM", "PerMM", "MDCK"),
                      "COSMOperm" = c("CACO-2", "PAMPA", "BLM", "PerMM", "MDCK"),
                      "CACO-2"= c("PAMPA", "BLM", "PerMM", "MDCK"),
                      "MDCK"= c("PerMM"))
  # Determines more complex graphs
  lm_model_filters= list("BLM_and_PERMM_COSMO" = c("COSMOperm", "PerMM"),
                         "PAMPA_and_MDCK_CACO" = c("EPAMOL", "CACO-2", "MDCK"))
  # Set dataset colors  
  colors          = list("PerMM" = 'goldenrod', 
                         "COSMOperm" = "forestgreen", 
                         "BLM" = "gray0", 
                         "PAMPA" = "deepskyblue", 
                         "CACO-2" = "deeppink3", 
                         "MDCK" = "darkmagenta")
}))

for (i in 1:length(config$method_groups))
{
  label <- config$method_groups[i]
  var_name = paste("d_",label, sep="")
  
  if(config$groups_col[i] == "1")
    group_values <- subset(data, method_label == label)
  else
    group_values <- subset(data, method == label)
  
  # Compute mean values of permeability coefficients in (cm/s)
  mean <- group_values %>%
    group_by(SMILES) %>%
    summarise(Perm_mean = mean(Perm))
  
  # Calculate logarithmic values of mean values
  mean$LogPerm_mean <- log(mean$Perm_mean, base=10)
  
  # Add column with method label
  mean$method_label <- config$method_labels[i]
  
  # Create dataset with original values and mean (and logarithimc) values
  mean <- merge(x=group_values, y=mean, by="SMILES", all.x=TRUE, all.y=TRUE)
  
  # Calculate logarithmic values of permeability values
  mean$LogPerm <- log(mean$Perm, base=10)
  
  # calculate absolute value of difference between logPerm and logPerm_mean
  mean$diff <- abs(mean$LogPerm - mean$LogPerm_mean)
  
  # If difference between logPerm and logPerm_mean > 2 use logPerm, else logPerm_mean
  mean$value_for_plot <- ifelse(mean$diff > 2,mean$LogPerm ,
                                    ifelse(mean$diff < 2,mean$LogPerm_mean, mean$LogPerm))
  
  # PAMPA methods - combine them...
  if (config$groups_col[i] == "2" && exists("d_PAMPA")) 
  {
    mean$PAMPA_label = label
    assign('d_PAMPA', rbind(d_PAMPA, mean))
  }
  else if(config$groups_col[i] == "2")
  {
    mean$PAMPA_label = label
    assign('d_PAMPA', mean)
  }
  
  # Finally, save as variable
  assign(var_name, mean)
}


# Combine MDCK and CACO-2 to one dataset
d_MDCK$color_label <- 'purple4'
`d_CACO-2`$color_label <- 'deeppink3'
d_MDCK_CACO <- rbind(d_MDCK, `d_CACO-2`);

# Combine PerMM and COSMOperm to one dataset
d_PerMM$color_label = "goldenrod"
d_COSMOperm$color_label = "forestgreen"
d_PERMM_COSMO <- rbind(d_PerMM, d_COSMOperm)


# Create datasets
for(ds1 in names(config$datasets_keys))
{
  var1 = paste('d_', ds1, sep='')
  for(ds2_i in 1:length(config$datasets_keys[[ds1]]))
  {
    ds2 <- config$datasets_keys[[ds1]][ds2_i]
    var2 <- paste('d_', ds2, sep='')
    
    if(!exists(var1) || !exists(var2))
    {
      print('DS1 or DS2 does not exist...');
      return()
    }
    
    # Combine two datasets
    comb <- merge(get(var1), get(var2), by.x="SMILES", by.y = "SMILES")
    
    # Add 2 extra points "x" into dataset because of scatter plot appearance
    comb <- comb %>% add_row(SMILES="x", value_for_plot.x = -16, value_for_plot.y = 7)
    comb <- comb %>% add_row(SMILES="x", value_for_plot.x = 7, value_for_plot.y = 16)
    
    # Color settings
    comb$color_of_point <- with(comb, ifelse(value_for_plot.x == -16 | value_for_plot.x == 7, "white", "black"))
    
    # Parse string as double
    comb$value_for_plot.x <- as.double(comb$value_for_plot.x)
    comb$value_for_plot.y <- as.double(comb$value_for_plot.y)
    
    # Add dashed lines
    comb$x_1_summ <- comb$value_for_plot.x + 1
    comb$x_1_diff <- comb$value_for_plot.x - 1
    
    assign(paste(ds1, ds2, sep="_and_"), comb)
  }
}


##########################################################
############ Plotting scatter plots ######################
##########################################################

get_r2 <- function(model_data)
{
  l_model = lm(value_for_plot.x~value_for_plot.y, data=model_data)
  return(round(summary(l_model)$r.squared,2))
}

### Start with simple plots
for(ds1 in names(config$datasets_keys))
{
  var1 = paste('d_', ds1, sep='')
  for(ds2_i in 1:length(config$datasets_keys[[ds1]]))
  {
    ds2 <- config$datasets_keys[[ds1]][ds2_i]
    var2 <- paste('d_', ds2, sep='')
    
    ds_name = paste(ds1, ds2, sep="_and_")
    
    if(!exists(ds_name))
    {
      print("Dataset does not exist...");
      return()
    }
    
    # Has restrictions? Then it is not simple chart => skip now
    if(ds_name %in% names(config$lm_model_filters))
    {
      next;   
    }
    
    # var_names = NULL
    model_data <- get(ds_name)
    model_data <- model_data [model_data$SMILES != "x",]
    
    # Compute linear model
    r_squared <- get_r2(model_data)
    total_rows = length(unique(model_data$SMILES))
    
    # Prepare plot data
    plotData <- ggplot(get(ds_name), aes(value_for_plot.x, value_for_plot.y,color = color_of_point)) +
      coord_cartesian(xlim = c(-16, 7.5), ylim = c(-16, 7.5)) +
      labs(title = "A")+
      theme_bw()+
      geom_point(size = 2)+
      annotate("text", label= deparse(bquote(R^list(2, n == .(total_rows)) ~"=" ~ .(r_squared))), x=-8, y=5, size = 8, parse = TRUE) +
      scale_color_identity()+
      scale_y_continuous(breaks = seq(-16, 7.5, by=4))+
      scale_x_continuous(breaks = seq(-16, 7.5, by=4))+
      labs(y = paste("LogPerm", ds2, "(cm/s)"), x = paste("LogPerm", ds1, "(cm/s)"))+
      geom_line(aes(x = value_for_plot.x, y = x_1_summ), linetype = "dashed", colour = "black")+
      geom_line(aes(x = value_for_plot.x, y = x_1_diff), linetype = "dashed", colour = "black")+
      geom_line(aes(x = value_for_plot.x, y = value_for_plot.x, colour = "black"))+
      theme(axis.title=element_text(size=29),
            axis.text=element_text(size=20),
            plot.title=element_text(size=30),
            axis.line.x = element_line(linewidth = 2, colour = config$colors[[ds1]], linetype="solid"),
            axis.line.y = element_line(linewidth = 2, colour = config$colors[[ds2]], linetype="solid"))
    
    # save finally to variable
    assign(paste('plot_', ds1,"_x_",ds2, sep=""), plotData)
  }
}
  
  
##### Now, plot two more complex graphs ######

# CACO-2 vs. MDCK vs. PAMPA
model_data = PAMPA_and_MDCK_CACO[PAMPA_and_MDCK_CACO$SMILES != 'x',]

# Distinct between apparent and intr. PAMPA
model_data_app = model_data[model_data$PAMPA_label != "EPAMOL", ]

r_squared <- get_r2(model_data_app)
total_rows = length(unique(model_data_app$SMILES))
# Inclusion of only MDCK
s_data = model_data_app[model_data_app$method_label.x.y != "CACO-2", ]
r_squared_MDCK = get_r2(s_data)
# Inclusion of only CACO-2
s_data = model_data_app[model_data_app$method_label.x.y != "MDCK", ]
r_squared_CACO = get_r2(s_data)

## Add points for proper visualisation
model_data_app = model_data_app %>% add_row(SMILES = "x", value_for_plot.x = -16, value_for_plot.y = 7, x_1_summ = -15, x_1_diff = -17)
model_data_app = model_data_app %>% add_row(SMILES = "x", value_for_plot.x = 7, value_for_plot.y = 16, x_1_summ = 8, x_1_diff = 6)


plot_PAMPA_x_MDCK_x_CACO <-ggplot(model_data_app, aes(value_for_plot.x, value_for_plot.y, color = color_label)) +
  geom_point(aes(color = color_label), size = 3) +
  coord_cartesian(xlim = c(-16, 7.5), ylim = c(-16, 7.5)) +
  scale_shape_manual(values=c(16, 3)) +
  labs(title = "C")+
  theme_bw()+
  geom_point(size = 2)+
  annotate("text", label= deparse(bquote(R[app]^list(2, n==.(total_rows))  ~"=" ~ .(r_squared))), x=-8, y=1.5, size = 8, parse = TRUE) +
  annotate("text", label= deparse(bquote(R[MDCK]^2 ~"=" ~ .(r_squared_MDCK))), x=0.5, y=-12, size = 8, parse = TRUE) +
  annotate("text", label= deparse(bquote(R[CACO-2]^2 ~"=" ~ .(r_squared_CACO))), x=0, y=-14.5, size = 8, parse = TRUE) +
  scale_color_identity()+
  scale_y_continuous(breaks = seq(-16, 7.5, by=4), sec.axis = dup_axis(name = "LogPerm MDCK (cm/s)" ))+
  scale_x_continuous(breaks = seq(-16, 7.5, by=4))+
  labs(y = "LogPerm CACO-2 (cm/s)", x = "LogPerm PAMPA (app.) (cm/s)")+
  geom_line(aes(x = value_for_plot.x, y = x_1_summ), linetype = "dashed", colour = "black", size = 0.75)+
  geom_line(aes(x = value_for_plot.x, y = x_1_diff), linetype = "dashed", colour = "black", size = 0.75)+
  geom_line(aes(x = value_for_plot.x, y = value_for_plot.x, colour = "black"))+
  theme(axis.title=element_text(size=29),
        axis.text=element_text(size=20),
        plot.title=element_text(size=30),
        axis.line.x = element_line(size = 2, colour = "lightskyblue", linetype="solid"),
        axis.line.y.left = element_line(size = 2, colour = "deeppink3", linetype="solid"),
        axis.line.y.right = element_line(size = 2, colour = "darkmagenta", linetype="solid"),
        legend.position = "none")
plot_PAMPA_x_MDCK_x_CACO
  
# Distinct between apparent and intr. PAMPA
model_data_int = model_data[model_data$PAMPA_label != "EPAM", ]

r_squared <- get_r2(model_data_int)
total_rows = length(unique(model_data_int$SMILES))
# Inclusion of only MDCK
s_data = model_data_int[model_data_int$method_label.x.y != "CACO-2", ]
r_squared_MDCK = get_r2(s_data)
# Inclusion of only CACO-2
s_data = model_data_int[model_data_int$method_label.x.y != "MDCK", ]
r_squared_CACO = get_r2(s_data)

## Add points for proper visualisation
model_data_int = model_data_int %>% add_row(SMILES = "x", value_for_plot.x = -16, value_for_plot.y = 7, x_1_summ = -15, x_1_diff = -17)
model_data_int = model_data_int %>% add_row(SMILES = "x", value_for_plot.x = 7, value_for_plot.y = 16, x_1_summ = 8, x_1_diff = 6)

  
plot_PAMPA_x_MDCK_x_CACO_int <-ggplot(model_data_int, aes(value_for_plot.x, value_for_plot.y, color = color_label)) +
  geom_point(aes(color = color_label), size = 3) +
  coord_cartesian(xlim = c(-16, 7.5), ylim = c(-16, 7.5)) +
  scale_shape_manual(values=c(16, 3)) +
  labs(title = "B")+
  theme_bw()+
  geom_point(size = 2)+
  annotate("text", label= deparse(bquote(R[intr]^list(2, n==.(total_rows))  ~"=" ~ .(r_squared))), x=-8, y=1.5, size = 8, parse = TRUE) +
  annotate("text", label= deparse(bquote(R[MDCK]^2 ~"=" ~ .(r_squared_MDCK))), x=0.5, y=-12, size = 8, parse = TRUE) +
  annotate("text", label= deparse(bquote(R[CACO-2]^2 ~"=" ~ .(r_squared_CACO))), x=0, y=-14.5, size = 8, parse = TRUE) +
  scale_color_identity()+
  scale_y_continuous(breaks = seq(-16, 7.5, by=4), sec.axis = dup_axis(name = "LogPerm MDCK (cm/s)" ))+
  scale_x_continuous(breaks = seq(-16, 7.5, by=4))+
  labs(y = "LogPerm CACO-2 (cm/s)", x = "LogPerm PAMPA (intr.) (cm/s)")+
  geom_line(aes(x = value_for_plot.x, y = x_1_summ), linetype = "dashed", colour = "black", size = 0.75)+
  geom_line(aes(x = value_for_plot.x, y = x_1_diff), linetype = "dashed", colour = "black", size = 0.75)+
  geom_line(aes(x = value_for_plot.x, y = value_for_plot.x, colour = "black"))+
  theme(axis.title=element_text(size=29),
        axis.text=element_text(size=20),
        plot.title=element_text(size=30),
        axis.line.x = element_line(size = 2, colour = "navy", linetype="solid"),
        axis.line.y.left = element_line(size = 2, colour = "deeppink3", linetype="solid"),
        axis.line.y.right = element_line(size = 2, colour = "darkmagenta", linetype="solid"),
        legend.position = "none")
plot_PAMPA_x_MDCK_x_CACO_int


# Cosmo vs. PerMM vs BLM
model_data = BLM_and_PERMM_COSMO[BLM_and_PERMM_COSMO$SMILES != 'x',]
r_squared <- get_r2(model_data)
total_rows = length(unique(model_data$SMILES))
# Inclusion of only PerMM
s_data = model_data[model_data$method_label.x.y != "COSMOperm", ]
r_squared_permm = get_r2(s_data)
total_permm = length(unique(s_data$SMILES))
# Inclusion of only COSMOPerm
s_data = model_data[model_data$method_label.x.y != "PerMM", ]
r_squared_cosmo = get_r2(s_data)
total_cosmo = length(unique(s_data$SMILES))

plot_BLM_x_PERMM_x_COSMO <-ggplot(BLM_and_PERMM_COSMO, aes(value_for_plot.x, value_for_plot.y,color = color_label)) +
  geom_point(aes(color = color_label, size = color_label), size = 3) +
  coord_cartesian(xlim = c(-16, 7.5), ylim = c(-16, 7.5)) +
  labs(title = "C")+
  theme_bw()+
  geom_point(size = 2)+
  annotate("text", label= deparse(bquote(R^list(2, n == .(total_rows)) ~"=" ~ .(r_squared))), x=-8, y=5, size = 8, parse = TRUE) +
  annotate("text", label= deparse(bquote(R[PerMM]^2 ~"=" ~ .(r_squared_permm))), x=1.5, y=-12, size = 8, parse = TRUE) +
  annotate("text", label= deparse(bquote(R[COSMOperm]^2 ~"=" ~ .(r_squared_cosmo))), x=0, y=-14.5, size = 8, parse = TRUE) +
  scale_color_identity()+
  scale_y_continuous(breaks = seq(-16, 7.5, by=4), sec.axis = dup_axis(name = "LogPerm COSMOperm (cm/s)" ))+
  scale_x_continuous(breaks = seq(-16, 7.5, by=4))+
  labs(y = "LogPerm PerMM (cm/s)", x = "LogPerm BLM (cm/s)")+
  geom_line(aes(x = value_for_plot.x, y = x_1_summ), linetype = "dashed", colour = "black", size = 0.75)+
  geom_line(aes(x = value_for_plot.x, y = x_1_diff), linetype = "dashed", colour = "black", size = 0.75)+
  geom_line(aes(x = value_for_plot.x, y = value_for_plot.x, colour = "black"))+
  theme(axis.title=element_text(size=29),
        axis.text=element_text(size=20),
        plot.title=element_text(size=30),
        axis.line.x = element_line(size = 2, colour = config$colors$BLM, linetype="solid"),
        axis.line.y.left = element_line(size = 2, colour = config$colors$PerMM, linetype="solid"),
        axis.line.y.right = element_line(size = 2, colour = config$colors$COSMOperm, linetype="solid"))

plot_BLM_x_PERMM_x_COSMO


# pdf(file="outputs/scatter plots/PAMPAapp_x_MDCK_CACO.pdf",
#     width=16,
#     height = 12)

# plot_PAMPA_x_MDCK_x_CACO

# get("plot_COSMOperm_x_CACO-2")

# dev.off()
