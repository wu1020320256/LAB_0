import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from torch import arange

# ylabel show %
# 纵坐标显示百分比
def toPercent(temp, position):
    return '%1.0f'%(100*temp) + '%'

# new legend(Nscale & SPR)
# 创建两个legend
def newLegend(legendName):
    tmpList = list(arange(len(legendName)))
    lineStyleList = ['solid','dashed']
    for i in range(len(legendName)):
        tmpList[i] = mlines.Line2D([], [], color='black', lw = 2.8, linestyle = lineStyleList[i], label = legendName[i])
    plt.legend(handles=[tmpList[0],tmpList[1]],loc='lower right',fontsize=13)

#画线
def linePlot(SIR, ax):
    tmpAdd = 6 * (SIR - (-5))/5   #根据SIR取值 选excel中的列 一组SIR对应六列
    pltList = [0 for i in range (8)]
    colSNR = dataRaw.iloc[:,[0]]
    dataRaw.iloc[:,[1]]
    markerList = ['o','^','s']
    colorList = ['tomato','lightskyblue','orange']
    lineStyleList = ['solid','dashed']
    for i in range(1,7):
        iNum = int((i-1)/2)
        labelStr = 'SF=' + '%d'%(8+iNum) if int((i-1)%2)==0 else None
        #plot
        pltList[i+1]=ax.plot(list(colSNR.values),list(dataRaw.iloc[:,[i+tmpAdd]].values), 
        marker=markerList[iNum], 
        color = colorList[iNum], 
        markersize = 9,
        label = labelStr,
        lw = 3.3,
        linestyle = lineStyleList[int((i-1)%2)]
        )

def figPLot(SIR = 5):
    fig, ax = plt.subplots()
    linePlot(SIR, ax)
    l1 = plt.legend(loc = 'best',fontsize=13)
    plt.gca().add_artist(l1)
    newLegend(['Nscale', 'CPR'])
    #yticks = ['0', '20%', '40%', '60%', '80%' '100%']
    plt.grid(ls= '-')
    plt.gca().yaxis.set_major_formatter(FuncFormatter(toPercent))
    plt.xlabel('SNR(dB)',fontsize=15)
    plt.ylabel('Collision Detection Accuracy',fontsize=15)  
    plt.tick_params(labelsize=13)
    plt.tight_layout()
    plt.savefig(r"exp_1_SIR%d.pdf"%SIR)


# dataRaw = pd.read_excel(r'C:/Users/ZKevin/Desktop/exp_1.xlsx',header=2)
# SIRList = [-5, 0, 5]
# for i in SIRList:
#     figPLot(SIR = i)

def linePlotTmp(ax):
    # tmpAdd = 6 * (SIR - (-5))/5   #根据SIR取值 选excel中的列 一组SIR对应六列
    pltList = [0 for i in range (5)]
    x = [1,2,3,4,5,6,7,8,9,10]
    # dataRaw.iloc[:,[1]]
    markerList = ['o','^','s','*']
    colorList = ['tomato','lightskyblue','orange','hotpink']
    # lineStyleList = ['solid','dashed']
    for i in range(1,5):
        labelStr = 'SF=' + '%d'%(6+i)
        colSNR = dataRaw.iloc[:,[i-1]] * 0.5;
        pltList[i]=ax.plot(x, colSNR.values, 
        marker = markerList[i-1], 
        color = colorList[i-1], 
        markersize = 9,
        label = labelStr,
        lw = 3.3,
        # linestyle = lineStyleList[int((i-1)%2)]
        linestyle = 'solid'
        )


def figPLotTmp():
    fig, ax = plt.subplots()
    linePlotTmp(ax)
    l1 = plt.legend(loc = 'best',fontsize=13)
    plt.gca().add_artist(l1)
    # newLegend(['Nscale', 'CPR'])
    x = [1,2,3,4,5,6,7,8,9,10]
    plt.xticks(x, ('-1000', '-800','-600','-400','-200','200','400','600','800','1000'))
    plt.grid(ls= '-')
    plt.ylim((0,100))
    # plt.gca().yaxis.set_major_formatter(FuncFormatter(toPercent))
    plt.xlabel('Hopping Frequency(KHz)',fontsize=15)
    plt.ylabel('Phase Jitter Time(us)',fontsize=15)  
    plt.tick_params(labelsize=13)
    # plt.tight_layout()
    plt.savefig(r"Schemtic_5.pdf")


dataRaw = pd.read_excel(r'C:/Users/ZKevin/Desktop/Schematic_5.xlsx')
figPLotTmp()

