"""

Summarized results for PyDotSim

by Christian Heyn

"""


import csv
import matplotlib.lines as lines
import lib.PyDotSim_Constants as cnst
import lib.PyDotSim_Model as model
import lib.PyDotSim_WF as WF


# colMap = 'gnuplot2'
# colMap = 'jet'
colMap = 'gray'

resultsHeader =     ['shape', 'r1', 'r2', 'r3', 'cells', 'cellSize', 'MP', 'MP#',
                     'meQ', 'mhQ', 'meB', 'mhB', 'Fz', 'Fx', 'Bz', 'Le', 'Lh', 'E0e', 'E0h',
                     'rxe', 'rye', 'rze', 'rxh', 'ryh', 'rzh', 'deh',
                     'overlap', 'tau', 'Ceh', ' Ex_single', 'Ep']
resultsSingle =     []
results =           []

Ee = [];    Eh =    []
Ex_single =         0
Ep =                0



# Add QD shape to WF plot
# line = [x1,x2], [y1,y2]
def addQDShape(ax):
    QDshape = model.QDpara[0]
    # Ellipse, so far only sphere
    if QDshape == 1:
        rxQD = 1e9*model.QDpara[1]; ryQD = 1e9*model.QDpara[2]; rzQD = 1e9*model.QDpara[3];
        pass
    # Disk, r1: r, r2: h
    if QDshape == 2:
        rQD = 1e9*model.QDpara[1]; hQD = 1e9*model.QDpara[2]
        pass
    # Cone, r1: r, r2: h
    if QDshape == 3:
        rQD = 1e9*model.QDpara[1]; hQD = 1e9*model.QDpara[2]
        pass
     # Cone shell, r1: r, r2: hHole, r3: hQD
    if QDshape == 4:
        rQD = model.QDpara[1]/model.cellSize; hHole = model.QDpara[2]/model.cellSize; hQD = model.QDpara[3]/model.cellSize
        ix0 = (model.cellNum1D-1)/2
        iz0 = ix0-hHole/2
        line = lines.Line2D([iz0, iz0+hHole], [ix0, ix0+rQD], lw=1, color='red')
        ax.add_line(line)
        line = lines.Line2D([iz0+hQD, iz0+hHole], [ix0, ix0+rQD], lw=1, color='red')
        ax.add_line(line)
        line = lines.Line2D([iz0+hQD, iz0+hHole], [ix0, ix0-rQD], lw=1, color='red')
        ax.add_line(line)
        line = lines.Line2D([iz0, iz0+hHole], [ix0, ix0-rQD], lw=1, color='red')
        ax.add_line(line)


# shows planes through the calculated wave functions
def showPSI(figResults, canvasResults):
    interpol = 'bicubic'
    #interpol = 'none'
    figResults.clf()
    i=model.cellNum1D//2
    psie = WF.psie; psih = WF.psih
    # probability density
    psie = psie**2; psih = psih**2
    figResults.add_subplot(231).imshow(psie[0][i , : , :], interpolation=interpol, cmap=colMap)
    figResults.add_subplot(232).imshow(psie[0][: , i , :], interpolation=interpol, cmap=colMap)
    figResults.add_subplot(233).imshow(psie[0][: , : , i], interpolation=interpol, cmap=colMap)
    figResults.add_subplot(234).imshow(psih[0][i , : , :], interpolation=interpol, cmap=colMap)
    figResults.add_subplot(235).imshow(psih[0][: , i , :], interpolation=interpol, cmap=colMap)
    figResults.add_subplot(236).imshow(psih[0][: , : , i], interpolation=interpol, cmap=colMap)
    # add QD shape
    addQDShape(figResults.axes[0])
    addQDShape(figResults.axes[1])
    addQDShape(figResults.axes[3])
    addQDShape(figResults.axes[4])
    canvasResults.draw()


def E0peak():
    global Ex_single, Ep
    T = 4
    Eg = cnst.GaAs_bandgap(T)
    Ex_single = Eg + abs(Ee[0]) + abs(Eh[0])
    Ex =  Ex_single - WF.Ceh
    particle = model.MultiParticleStr[model.MPpara[0]]
    if particle == 'X':     Ep = Ex
    if particle == 'XX':    Ep = Ex - (2*WF.Ceh - WF.Cee - WF.Chh)
    if particle == 'X+':    Ep = Ex - (WF.Ceh - WF.Chh)
    if particle == 'X-':    Ep = Ex - (WF.Ceh - WF.Cee)


def createResultsSingle():
    global resultsSingle, results
    resultsSingle = []
    resultsSingle.append( model.QDshapeStr[model.QDpara[0]-1] )
    resultsSingle.append( str(model.QDpara[1]*1e9) )
    resultsSingle.append( str(model.QDpara[2]*1e9) )
    resultsSingle.append( str(model.QDpara[3]*1e9) )
    resultsSingle.append( str(model.cellNum1D) )
    resultsSingle.append( str(model.cellSize*1e9) )
    resultsSingle.append( model.MultiParticleStr[model.MPpara[0]] )
    resultsSingle.append( str(model.MultiParticleLoops) )
    resultsSingle.append( str(model.meQ) )
    resultsSingle.append( str(model.mhQ) )
    resultsSingle.append( str(model.meB) )
    resultsSingle.append( str(model.mhB) )
    resultsSingle.append( str(model.EFields[0]) )
    resultsSingle.append( str(model.EFields[1]) )
    resultsSingle.append( str(model.BZfield[0]) )
    resultsSingle.append( str(model.BZfield[1]) )
    resultsSingle.append( str(model.BZfield[2]) )
    resultsSingle.append( str(Ee[0]) )
    resultsSingle.append( str(Eh[0]) )
    resultsSingle.append( str(WF.rxe*1e9) )
    resultsSingle.append( str(WF.rye*1e9) )
    resultsSingle.append( str(WF.rze*1e9) )
    resultsSingle.append( str(WF.rxh*1e9) )
    resultsSingle.append( str(WF.ryh*1e9) )
    resultsSingle.append( str(WF.rzh*1e9) )
    resultsSingle.append( str(WF.deh*1e9) )
    resultsSingle.append( str(WF.overlap) )
    resultsSingle.append( str(WF.tau) )
    resultsSingle.append( str(WF.Ceh) )
    resultsSingle.append( str(str(Ex_single)) )
    resultsSingle.append( str(str(Ep)) )
    results.append(resultsSingle)


def initResults():
    global results
    results = []
    results.append(resultsHeader)


def saveResults(FileName_data):
    #FileName_data = os.getcwd()+'/'+name+', '+datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+'.dat'
    with open(FileName_data, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(results)



