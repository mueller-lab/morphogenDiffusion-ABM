import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from datetime import datetime
import pandas as pd
import os
import time

plt.rcParams["font.family"] = "Times New Roman"
cwdPath=Path(os.path.abspath(os.getcwd()))
fPath = cwdPath/"data"/"forScatterPlot.csv"
parDf = pd.read_csv(fPath)
ltGrey=  (10,10,10)

# make a 3D-scatter plot of parameter-space color coded with output BoundFrac
def plot3D_scatter(x,y,z,m,mCol,lbl,angle, dataPath):
    csfont = {'fontname':'Times New Roman'}
    fig = plt.figure(facecolor=None)

    ax = fig.add_subplot(111, projection='3d')
    # ax.set_facecolor(ltGrey)
    if 'CbyI' in mCol:
        p = ax.scatter(x, y, z, c=m, marker ='o', vmin=0.5, vmax=2.2, s=80, edgecolors= 'none', depthshade=False, cmap= "jet") #*np.array(m)
    elif 'sptD' in mCol:
        p = ax.scatter(x, y, z, c=m, marker ='o', vmin=0, vmax=max(m), s=80, edgecolors= 'none', depthshade=False, cmap= "jet")
    elif 'avD' in mCol:
        p = ax.scatter(x, y, z, c=m, marker ='o', vmin=0, vmax=max(m), s=80, edgecolors= 'none', depthshade=False, cmap= "jet")
    else:
        p = ax.scatter(x, y, z, c=m, marker ='o', vmin=0, vmax=max(m), s=80, edgecolors= 'none', depthshade=False, cmap= "jet")

    cbar= fig.colorbar(p, ax=ax, aspect = 40, shrink= 0.6)
    cbar.set_label(lbl)
    ax.set_xlabel(r'$Narrowness$', **csfont)
    ax.set_ylabel(r'$Receptor\/ Density\/ (log_{10}\/ (\mu m^{-1}))$', **csfont)
    ax.set_zlabel(r'$Residence\/ Time\/ (s)$', **csfont)
    # start, end = ax.get_ylim()

    # ax.set_xlim((0,180))
    ax.set_xlim((-9*20, 180))

    ax.set_ylim((-1.31,0.9))
    ax.set_zlim((0,17))

    # ax.xaxis.set_ticks(np.arange(0, 170, 80))
    ax.xaxis.set_ticks(np.arange(-8*20, 170, 80))

    ax.yaxis.set_ticks(np.arange(-1.31, 0.75, 0.5))
    ax.zaxis.set_ticks(np.arange(0, 17, 4))

    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.view_init(10, angle)
    plt.tight_layout()
    plt.savefig(dataPath/f'scatter_{mCol}.svg')
    plt.close()

# parDf= parDf[parDf['nrw']<9]
# parDf= parDf[parDf['nrw']>(-1)]

parDf= parDf[parDf['nrw']<16]
parDf= parDf[parDf['nrw']>(-9)]
# print(f"head of df is {parDf.head}")

x = parDf['nrw'].tolist()
x = [20*a for a in x]
y = parDf['recepDens'].tolist()
z = parDf['rTm'].tolist()


mCols=['bFrac', 'bFracC', 'bFracI', 'jmpAv', 'jmpAvC', 'jmpAvI', 'sptDC', 'sptDI', 'sptCbyI', \
'avD', 'avD_C', 'avD_I']
colLbl = [r'$Bound\/ fraction$', r'$Bound\/ fraction\/ (C)$', r'$Bound\/ fraction\/ (I)$', r'$Mean\/ jump-distance \/(\mu m)$', \
r'$Mean\/ jump-distance\/ (C) \/(\mu m)$', r'$Mean\/ jump-distance\/ (I) \/(\mu m)$', r'$\rho_c$', r'$\rho_i$', r'$\frac{\rho_c}{\rho_i}$', \
r'$D \/(\frac{\mu m^2}{s})$', r'$D_c \/(\frac{\mu m^2}{s})$', r'$D_i \/(\frac{\mu m^2}{s})$']

for mCol, lbl in zip(mCols,colLbl):
    m= parDf[mCol].tolist()
    plot3D_scatter(x,y,z,m,mCol,lbl,16,cwdPath)
