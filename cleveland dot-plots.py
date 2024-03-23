import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


# Data columns which should be present in the input data
required_columns = ['MDCK', 'CACO', 'PAMPA', 'COSMO', 'PERMM', 'BLM' ]

# Figures size
figsize = (28,12)

def __main__():
    # Load data
    data = pd.read_csv("inputs/cleveland_plot_dataset.csv", sep= ",", encoding="UTF-8")
    
    # Hold final figures in a list
    figs = []

    # Split data into two datasets
    BLM_subplot = data[~data['BLM'].isna()]
    main_plot = data[~data['BLM'].notna()]
    
    # sort values alphabetically
    BLM_subplot = BLM_subplot.sort_values("name", ascending = False)
    main_plot = main_plot.sort_values("name", ascending = False)
    
    # Identify max nad min values
    max_values_main = main_plot[required_columns].max(axis=1)
    min_values_main = main_plot[required_columns].min(axis=1)

    max_values_BLM = BLM_subplot[required_columns].max(axis=1)
    min_values_BLM = BLM_subplot[required_columns].min(axis=1)
    
    # And save molecule names
    Names_main = main_plot['name'].to_list()
    Names_BLM = BLM_subplot['name'].to_list()
    
    # set x-axis
    list_of_xticks = np.arange(-14, 8, 2)
    
    # main subplot
    fig_main = plt.figure(figsize=figsize)
    plt.grid()
    plt.hlines(y=Names_main, xmin=min_values_main.to_list(), xmax=max_values_main.to_list(), color = "black", label='_nolegend_', linewidth=3)

    plt.scatter(main_plot["MDCK"].to_list(), Names_main, s=300, alpha=1, color = "darkmagenta")
    plt.scatter( main_plot["CACO"].to_list(), Names_main, s=300, alpha=1, color = "deeppink")
    plt.scatter(main_plot["PAMPA"].to_list(), Names_main, s=300, alpha=1, color = "deepskyblue")
    plt.scatter(main_plot["PERMM"].to_list(), Names_main, s=300, alpha=1, color = "goldenrod")
    plt.scatter(main_plot["COSMO"].to_list(), Names_main, s=300, alpha=1, color = "forestgreen" )

    plt.legend(["MDCK", "CACO-2", "PAMPA",'PerMM',"COSMOperm"], loc ="lower left", fontsize = 15)
    plt.xticks(list_of_xticks)
    plt.yticks(Names_main)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.title("B", fontsize = 30, loc = "left")
    plt.xlabel('LogPerm (cm/s)', fontsize=20, labelpad = 10)
    plt.ylabel('Name of molecule', fontsize=20)

    plt.axvline(x=-12, color='dimgray',linewidth=2, linestyle = "solid")
    plt.axvline(x=-8, color='dimgray',linewidth=1.25, linestyle = "solid")
    plt.axvline(x=-4, color='dimgray',linewidth=1.25, linestyle = "dashdot")
    plt.axvline(x=0, color='dimgray',linewidth=1.25, linestyle = "dashed")
    plt.axvline(x=4, color='dimgray',linewidth=1.25, linestyle = "dotted")
    
    figs.append(fig_main)
    plt.close()
    
    # BLM subplot
    fig_blm = plt.figure(figsize=figsize)
    plt.grid()
    plt.hlines(y=Names_BLM, xmin=min_values_BLM.to_list(), xmax=max_values_BLM.to_list(), color = "black", label='_nolegend_', linewidth=3)

    plt.scatter(BLM_subplot["MDCK"].to_list(), Names_BLM, s=300, alpha=1, color = "darkmagenta")
    plt.scatter(BLM_subplot["CACO"].to_list(), Names_BLM, s=300, alpha=1, color = "deeppink")
    plt.scatter(BLM_subplot["PAMPA"].to_list(), Names_BLM, s=300, alpha=1, color = "deepskyblue")
    plt.scatter(BLM_subplot["PERMM"].to_list(), Names_BLM, s=300, alpha=1, color = "goldenrod")
    plt.scatter(BLM_subplot["BLM"].to_list(), Names_BLM, s=300, alpha=1, color = "dimgray")
    plt.scatter(BLM_subplot["COSMO"].to_list(), Names_BLM, s=300, alpha=1, color = "forestgreen" )


    plt.legend(["MDCK", "CACO-2", "PAMPA",'PerMM',"BLM / Liposomes","COSMOperm"], loc ="lower left", ncol = 1, fontsize = 15)
    plt.xticks(list_of_xticks)
    plt.yticks(Names_BLM)
    plt.xlabel('LogPerm (cm/s)', fontsize=20, labelpad= 10)
    plt.ylabel('Name of molecule', fontsize=20)
    fig_blm.text(0.5, 1.1, "Time", fontsize = 20, horizontalalignment="center")
    fig_blm.text(0.17, 0.95, "Weeks",fontsize = 15, horizontalalignment="center")
    fig_blm.text(0.28, 0.95, "Days",fontsize = 15, horizontalalignment="center")
    fig_blm.text(0.42, 0.95, "Hours",fontsize = 15, horizontalalignment="center")
    fig_blm.text(0.56, 0.95, "Seconds",fontsize = 15, horizontalalignment="center")
    fig_blm.text(0.77, 0.95, "Microseconds",fontsize = 15, horizontalalignment="center")
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.title("A", fontsize = 30, loc = "left", pad = 50)
    #plt.rcParams["font.family"] = "Arial"
    plt.margins(y = 0.5)

    plt.axvline(x=-12, color='dimgray',linewidth=2, linestyle = "solid")
    plt.axvline(x=-8, color='dimgray',linewidth=1.25, linestyle = "solid")
    plt.axvline(x=-4, color='dimgray',linewidth=1.25, linestyle = "dashdot")
    plt.axvline(x=0, color='dimgray',linewidth=1.25, linestyle = "dashed")
    plt.axvline(x=4, color='dimgray',linewidth=1.25, linestyle = "dotted")
    
    figs.append(fig_blm)
    plt.close()
    
    
    # Finally, save figures to PDF
    with PdfPages('outputs/cleveland-plots/all.pdf') as pdf:
        names = ['main', 'blm']
        for fig in figs:
            pdf.savefig(fig)
            # Plus, save each figure separately
            fig.savefig('outputs/cleveland-plots/' + names.pop(0) + '.png')
__main__()