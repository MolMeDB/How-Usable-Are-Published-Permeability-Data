import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches


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

    color_list_main_plot = []

    # Iterate through each value in X
    for value in main_plot["method"].to_list():
        if value == 'EPAMOL':
            color_list_main_plot.append('navy')
        elif value == 'EPAM':
            color_list_main_plot.append('skyblue')

    color_list_sub_plot = []

    # Iterate through each value in X
    for value in BLM_subplot["method"].to_list():
        if value == 'EPAMOL':
            color_list_sub_plot.append('navy')
        elif value == 'EPAM':
            color_list_sub_plot.append('skyblue')
    
    # set x-axis
    list_of_xticks = np.arange(-14, 8, 2)
    
    # main subplot
    fig_main = plt.figure(figsize=figsize)
    plt.grid()
    plt.hlines(y=Names_main, xmin=min_values_main.to_list(), xmax=max_values_main.to_list(), color = "black", label='_nolegend_', linewidth=3)

    plt.scatter(main_plot["MDCK"].to_list(), Names_main, s=300, alpha=1, color = "darkmagenta")
    plt.scatter(main_plot["CACO"].to_list(), Names_main, s=300, alpha=1, color = "deeppink")
    plt.scatter(main_plot["PAMPA"].to_list(), Names_main, s=300, alpha=1, color = color_list_main_plot)
    plt.scatter(main_plot["PERMM"].to_list(), Names_main, s=300, alpha=1, color = "goldenrod")
    plt.scatter(main_plot["COSMO"].to_list(), Names_main, s=300, alpha=1, color = "forestgreen" )

    handles, labels = plt.gca().get_legend_handles_labels()

    patch_MDCK = mpatches.Patch(color='darkmagenta', label='MDCK')
    patch_CACO = mpatches.Patch(color='deeppink', label='CACO-2')
    patch_PAMPA_A = mpatches.Patch(color='skyblue', label='PAMPA app.')
    patch_PAMPA_M = mpatches.Patch(color='navy', label='PAMPA intr.')
    patch_PerMM = mpatches.Patch(color='goldenrod', label='PerMM')
    patch_COSMO = mpatches.Patch(color='forestgreen', label='COSMOperm')

    handles.extend([patch_MDCK, patch_CACO, patch_PAMPA_A, patch_PAMPA_M, patch_PerMM, patch_COSMO])

    plt.legend(handles=handles, loc ="lower left", fontsize = 15)

    plt.xticks(list_of_xticks)
    plt.yticks(Names_main)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.title("B", fontsize = 30, loc = "left")
    plt.xlabel('LogPerm (cm/s)', fontsize=20, labelpad = 10)
    plt.ylabel('Name of molecule', fontsize=20)

    figs.append(fig_main)
    plt.close()
    
    # BLM subplot
    fig_blm = plt.figure(figsize=figsize)
    plt.grid()
    plt.hlines(y=Names_BLM, xmin=min_values_BLM.to_list(), xmax=max_values_BLM.to_list(), color = "black", label='_nolegend_', linewidth=3)

    plt.scatter(BLM_subplot["MDCK"].to_list(), Names_BLM, s=300, alpha=1, color = "darkmagenta")
    plt.scatter(BLM_subplot["CACO"].to_list(), Names_BLM, s=300, alpha=1, color = "deeppink")
    plt.scatter(BLM_subplot["PAMPA"].to_list(), Names_BLM, s=300, alpha=1, color = color_list_sub_plot)
    plt.scatter(BLM_subplot["PERMM"].to_list(), Names_BLM, s=300, alpha=1, color = "goldenrod")
    plt.scatter(BLM_subplot["BLM"].to_list(), Names_BLM, s=300, alpha=1, color = "dimgray")
    plt.scatter(BLM_subplot["COSMO"].to_list(), Names_BLM, s=300, alpha=1, color = "forestgreen" )

    handles, labels = plt.gca().get_legend_handles_labels()

    patch_MDCK = mpatches.Patch(color='darkmagenta', label='MDCK')
    patch_CACO = mpatches.Patch(color='deeppink', label='CACO-2')
    patch_PAMPA_A = mpatches.Patch(color='skyblue', label='PAMPA app.')
    patch_PAMPA_M = mpatches.Patch(color='navy', label='PAMPA intr.')
    patch_PerMM = mpatches.Patch(color='goldenrod', label='PerMM')
    patch_BLM = mpatches.Patch(color='dimgray', label='BLM / Liposomes')
    patch_COSMO = mpatches.Patch(color='forestgreen', label='COSMOperm')

    handles.extend([patch_MDCK, patch_CACO, patch_PAMPA_A, patch_PAMPA_M, patch_PerMM, patch_BLM, patch_COSMO])

    plt.legend(handles=handles, loc ="lower left", fontsize = 15)

    plt.xticks(list_of_xticks)
    plt.yticks(Names_BLM)
    plt.xlabel('LogPerm (cm/s)', fontsize=20, labelpad= 10)
    plt.ylabel('Name of molecule', fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.title("A", fontsize = 30, loc = "left")
    plt.margins(y = 0.5)
    
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