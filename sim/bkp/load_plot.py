"""
script to load sim and plot
"""

from netpyne import sim
from matplotlib import pyplot as plt
import os
import IPython as ipy
import pickle as pkl



poptypeNumber = 61 # max 55 + 6
celltypeNumber = 213 # max 207 + 6

# TO DEBUG - import and simulate only the Cell soma (to study only the Net)
reducedtest = False    

#------------------------------------------------------------------------------  
#------------------------------------------------------------------------------  
# S1 Cells
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop

with open('cells/S1-cells-distributions-Rat.txt') as mtype_file:
    mtype_content = mtype_file.read()       

popNumber = {}
cellNumber = {} 
popLabel = {} 
popParam = []
cellParam = []
meParamLabels = {} 
popLabelEl = {} 
cellLabel = {}

RP_L13 = []
RP_L45 = []
RP_L6 = []

for line in mtype_content.split('\n')[:-1]:
    cellname, mtype, etype, n, m = line.split()
    metype = mtype + '_' + etype[0:3]
    cellNumber[metype] = int(n)
    popLabel[metype] = mtype
    popNumber[mtype] = int(m)
    cellLabel[metype] = cellname

    if mtype not in popParam:
        popParam.append(mtype)
        popLabelEl[mtype] = [] 
               
    popLabelEl[mtype].append(metype)
    
    cellParam.append(mtype + '_' + etype[0:3])

    layernumber = float(metype[1:2])
    if cellNumber[metype]*0.01 > 0.0:
        if int(layernumber) <= 3:
            RP_L13.append(mtype + '_' + etype[0:3])
            # print(layernumber,int(layernumber),mtype + '_' + etype[0:3])
        elif int(layernumber) == 6:
            RP_L6.append(mtype + '_' + etype[0:3])
            # print(layernumber,int(layernumber),mtype + '_' + etype[0:3])
        else:
            RP_L45.append(mtype + '_' + etype[0:3])
            # print(layernumber,int(layernumber),mtype + '_' + etype[0:3])
    
S1pops = popParam[0:55]
S1cells = cellParam[0:207]


###########################
######## MAIN CODE ########
###########################

if __name__ == '__main__':

    dataType = 'spont' #'speech' #'spont'

    if dataType == 'spont':
        filenames = ['../data/v7_batch1/v7_batch1_%d_%d_data.pkl' % (iseed, cseed) for iseed in [0] for cseed in [0]]
        timeRange = [0, 15000]

    layer_bounds= {'L1': 100, 'L2': 160, 'L3': 950, 'L4': 1250, 'L5A': 1334, 'L5B': 1550}#, 'L6': 2000}
    layer_bounds= {'S': 950, 'G': 1250, 'I': 1900}#, 'L6': 2000}


    allData = []

    for filename in filenames:
        sim.load(filename, instantiate = False)
        # sim.load(filename, instantiate=True, instantiateConns = False, instantiateStims = False, instantiateRxD = False, createNEURONObj = False)

        # standardd plots
        # sim.analysis.plotRaster(**{'include': ['allCells'], 'saveFig': True, 'showFig': False, 'labels': None, 'popRates': False,'orderInverse': True, 'timeRange': timeRange, 'figSize': (36,24), 'fontSize':4, 'lw': 5, 'markerSize':10, 'marker': '.', 'dpi': 300})
        # sim.analysis.plotRaster(**{'include': RP_L13, 'saveFig': True, 'showFig': False, 'popRates': 'minimal', 'orderInverse': True, 'timeRange': timeRange, 'orderBy':'y', 'fontSize':8, 'figSize': (24,12), 'lw': 4.0, 'markerSize': 4, 'marker': 'o', 'dpi': 300})
        # timeRange = [10000, 15000]
        sim.analysis.plotRaster(**{'saveFig': True, 'showFig': False, 'labels': None, 'popRates': False, 'orderInverse': True, 'timeRange': timeRange, 'orderBy':'y', 'fontSize':4, 'figSize': (60,30), 'lw': 1.0, 'markerSize': 1, 'marker': 'o', 'dpi': 300})
        # timeRange = [0, 5000]
        # sim.analysis.plotRaster(**{'include': RP_L6, 'saveFig': True, 'showFig': False, 'popRates': 'minimal', 'orderInverse': True, 'timeRange': timeRange, 'orderBy':'y', 'fontSize':8, 'figSize': (48,24), 'lw': 2.0, 'markerSize': 2, 'marker': 'o', 'dpi': 300})
        # sim.analysis.plotRaster(**{'include': RP_L6, 'saveFig': filename[:-4]+'_RP_L6', 'showFig': False, 'popRates': 'minimal', 'orderInverse': True, 'timeRange': timeRange, 'orderBy':'y', 'fontSize':8, 'figSize': (24,12), 'lw': 4.0, 'markerSize': 4, 'marker': 'o', 'dpi': 300})
        # sim.analysis.plotSpikeStats(stats=['rate'],figSize = (6,12), timeRange=[1500, 31500], dpi=300, showFig=0, saveFig=filename[:-4]+'_stats_30sec')
        #sim.analysis.plotSpikeStats(stats=['rate'],figSize = (6,12), timeRange=[1500, 6500], dpi=300, showFig=0, saveFig=filename[:-4]+'_stats_5sec')
        #sim.analysis.plotLFP(**{'plots': ['spectrogram'], 'electrodes': ['avg', [0], [1], [2,3,4,5,6,7,8,9], [10, 11, 12], [13], [14, 15], [16,17,18,19]], 'timeRange': timeRange, 'maxFreq': 50, 'figSize': (8,24), 'saveData': False, 'saveFig': filename[:-4]+'_LFP_spec_7s_all_elecs', 'showFig': False})
        # sim.analysis.plotRaster(**{'include': ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1'], 'saveFig': True, 'showFig': False, 'popRates': 'minimal', 'orderInverse': True, 'timeRange': timeRange, 'orderBy':'y', 'fontSize':8, 'figSize': (48,24), 'lw': 2.0, 'markerSize': 2, 'marker': 'o', 'dpi': 300})
        # sim.analysis.plotRaster(**{'include': S1cells, 'saveFig': True, 'showFig': False, 'labels': None, 'popRates': False,'orderInverse': True, 'timeRange': timeRange, 'figSize': (36,24), 'fontSize':4, 'lw': 5, 'markerSize':10, 'marker': '.', 'dpi': 300})
        
        # sim.analysis.plotLFP(**{'plots': ['locations'], 
        #         'figSize': (24,24), 
        #         'saveData': False, 
        #         'saveFig': True, 'showFig': False, 'dpi': 300})

        # sim.analysis.plotLFP(**{'plots': ['timeSeries'], 
        #         # 'electrodes': 
        #         # [[0,1,2,3]], #'avg', 
        #         'timeRange': timeRange, 
        #         'figSize': (24,12), 'saveFig': True, 'showFig': False})

        # sim.analysis.plotLFP(**{'plots': ['spectrogram'], 
        #         # 'electrodes': 
        #         # [[0,1,2,3]],
        #         'timeRange': timeRange, 
        #         'maxFreq': 400, 
        #         'figSize': (16,12), 
        #         'saveData': False, 
        #         'saveFig': True, 'showFig': False})

        # sim.analysis.plotLFP(**{'plots': ['PSD'], 
        #         # 'electrodes': 
        #         # [[0,1,2,3]],
        #         'timeRange': timeRange, 
        #         'maxFreq': 400, 
        #         'figSize': (8,12), 
        #         'saveData': False, 
        #         'saveFig': True, 'showFig': False})
