"""

OS related procedures for PyDotSim

by Christian Heyn

"""

import os
import datetime
import csv
import numpy as np



baseFolder = ''
tempFolder = ''
dataFolder = ''


def checkFolders():
    global baseFolder, tempFolder, dataFolder
    baseFolder = os.getcwd()+'/'
    tempFolder = baseFolder + '_temp/'
    dataFolder = baseFolder + '_data/'
    if os.path.exists(tempFolder) == False:
        os.makedirs(tempFolder)
    if os.path.exists(dataFolder) == False:
        os.makedirs(dataFolder)


def saveResults(data,name):
    FileName_data = os.getcwd()+'/'+name+', '+datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+'.dat'
    with open(FileName_data, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(results)


def saveWF(WF,eh,CIloop):
    if eh == 'e': fileName = tempFolder+'WFe'+str(CIloop)+'.npy'
    if eh == 'h': fileName = tempFolder+'WFh'+str(CIloop)+'.npy'
    np.save(fileName, WF)


def saveWFflat(WFflat,eh,CIloop):
    if eh == 'e': fileName = tempFolder+'WFeflat'+str(CIloop)+'.txt'
    if eh == 'h': fileName = tempFolder+'WFhflat'+str(CIloop)+'.txt'
    with open(fileName, 'w') as f:
        for i in range(len(WFflat)):
            f.write(str(WFflat[i]))
            f.write('\n')


def loadWF(eh,CIloop):
    if eh == 'e': fileName = tempFolder+'WFe'+str(CIloop)+'.npy'
    if eh == 'h': fileName = tempFolder+'WFh'+str(CIloop)+'.npy'
    WF = np.load(fileName)
    return WF
