############ Loading libraries #####################################
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
import yaml
from matplotlib.patches import Patch

################## Defs to run the script########################################
def compact(df, blocks='block_id', label='c80e3', colormap=None, legend_color=None, title=None, feature_type=['CDS','tRNA','rRNA','tmRNA','PSE']):
    begin = df.groupby('block_id').start.min().to_dict()
    begin = df.block_id.map(begin) - 1
    arrowstyle = mpatches.ArrowStyle.Simple(head_width=14, tail_width=14, head_length=max(0.6*14,5))

    # Setting row colors
    color = '#dddddd'
    if isinstance(colormap, dict):
        color = df[label].map(colormap)
    color = color.fillna('#dddddd')

    # Setting drawing area limits and options
    fig, ax = plt.subplots(figsize=(25, df[blocks].nunique()/1.3))
    ax.set_xlim(0, (df.end - begin).max() + 5)
    ax.set_ylim(0,df[blocks].nunique()+5)
    if title:
        plt.title(title)
    plt.box(False)

    y = 1
    previousBlockID = df.iloc[0][blocks]
    arrows = []
    for idx in range(0,len(df)):
        if previousBlockID != df.iloc[idx][blocks]:
            previousQuery = df[(df.block_id == previousBlockID) & (df['query'] == 1)].drop_duplicates('block_id')
            if not previousQuery.empty:
                previousQuery = previousQuery.iloc[0]
                ax.text(1, y - 0.6, f'{previousQuery.pid} ({previousBlockID}, {previousQuery.organism})')
            y = y + 1
        row = df.iloc[idx]
        if row.type not in feature_type:
            continue
        if row.strand == 1:
            start = row.start - begin.iloc[idx]
            end   = row.end - begin.iloc[idx]
        else:
            end   = row.start - begin.iloc[idx]
            start = row.end - begin.iloc[idx]
        arrow = mpatches.FancyArrowPatch(
            posA=[start, y],
            posB=[end, y],
            shrinkA=0.0, shrinkB=0.0,
            arrowstyle=arrowstyle,
            facecolor=color.iloc[idx],
            zorder=0,
            edgecolor="#000000",
            linewidth=1.0,
            url= f'https://www.ncbi.nlm.nih.gov/protein/{row.pid}',
            label = row[label],
        )
        ax.add_patch(arrow)
        previousBlockID = row[blocks]
    ax.text(1, y - 0.7, f'{df.iloc[-1].pid} ({df.iloc[-1].block_id}, {df.iloc[-1].organism})')

    # Legend
    ax.get_yaxis().set_visible(False)
    ax.axhline(y=0, color="black")
    #frame1.axes.get_xaxis().set_visible(False)
    legend_elements = [ Patch(facecolor=x,edgecolor='black',label=legend_color[x]) for x in legend_color ]
    ax.legend(handles=legend_elements, loc='upper center', ncol=5)

    return fig, ax

#fig, ax = patches(Firmicutes, label='component', colormap=color_d, legend_color=legend_color)
#plt.savefig('./t6ss_neighborhood.04082024_5.svg', format="svg")
