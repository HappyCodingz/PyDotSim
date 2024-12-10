"""

PyDotSim

Computes the quantized energy states of various QD shapes

by Christian Heyn

"""

import timeit
import os
import datetime

import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

import lib.PyDotSim_OS as osprocs
import lib.PyDotSim_FrameMaterials as frameMat
import lib.PyDotSim_FramePointCharges as framePC
import lib.PyDotSim_Constants as cnst
import lib.PyDotSim_Model as model
import lib.PyDotSim_WF as WF
import lib.PyDotSim_Schroedinger as solver
import lib.PyDotSim_Results as results


sTitle = 'PyDotSim (V09.12.2024), by Christian Heyn'

window = tk.Tk()
screenwidth = window.winfo_screenwidth()
screenheight = window.winfo_screenheight()
window.title(sTitle)
window_width = screenwidth-20; window_height = screenheight-100;
window.geometry( str(window_width)+'x'+str(window_height)+'+5+5')             # Window size

# Frames
frameMenue = tk.Frame(window)
frameMenue.pack(side='left', anchor='n', padx=10)
frameModelResults = tk.Frame(window)
frameModelResults.pack(side='left', anchor='n')
frameModel = tk.Frame(frameModelResults)
frameModel.pack(side='top', anchor='w')
frameResults = tk.Frame(frameModelResults)
frameResults.pack(side='top', anchor='w')
frameLog = tk.Frame(window)
frameLog.pack(side='left', anchor='n')

#figModel = Figure(figsize=(12,3.5))
figModel = Figure(figsize=(8,2))
canvasModel = FigureCanvasTkAgg(figModel, frameModel)
canvasModel.get_tk_widget().pack(side=tk.TOP, fill=tk.X)
#figResults = Figure(figsize=(12,7))
figResults = Figure(figsize=(8,4))
canvasResults = FigureCanvasTkAgg(figResults, frameResults)
canvasResults.get_tk_widget().pack(side=tk.TOP, fill=tk.X)


def doSchroedinger():
    Ee, psie = solver.calcSchroedingerFEM('e')
    Eh, psih = solver.calcSchroedingerFEM('h')
    psie = solver.correctInvertedPSI(psie);
    psih = solver.correctInvertedPSI(psih)
    results.Ee = Ee
    results.Eh = Eh
    WF.psie = psie
    WF.psih = psih


def computeSingle():
    start = timeit.default_timer()
    # first run without self-cons.
    model.createModel()
    doSchroedinger()
    # with self-cons.
    if model.MultiParticleLoops > 0:
        for loop in range(model.MultiParticleLoops):
            WFe = WF.normalizeWF(WF.psie[0])
            WFh = WF.normalizeWF(WF.psih[0])
            osprocs.saveWF(WFe,'e',loop)
            osprocs.saveWF(WFh,'h',loop)
            model.createModel()
            WF.addWF2pot(loop)
            model.ShowPotential(figModel, canvasModel)
            doSchroedinger()
        model.ShowPotential(figModel, canvasModel)
        results.showPSI(figResults, canvasResults)
    stop = timeit.default_timer()
    results.E0peak(); WF.E_PL = results.Ex_single
    runTime = stop - start
    return runTime


def doCompute():
    figResults.clf(); canvasResults.draw()
    model.MultiParticleLoops = int(entryMPloops_str.get())
    text_area.insert(tk.INSERT, '\n>> Solver started ('+model.MultiParticleStr[model.MPpara[0]]+': #'+str(model.MultiParticleLoops)+')')
    window.update()
    runTime = computeSingle()
    text_area.insert(tk.INSERT, '\n   Run time: '+str(round(runTime,3))+' s' )
    s = 'e: '
    for i in range(model.statesNum):
        s = s+'{0: >#8.6f}'.format(results.Ee[i])+', '
    text_area.insert(tk.INSERT, '\n   '+s)
    window.update()
    s = 'h: '
    for i in range(model.statesNum):
        s = s+'{0: >#8.6f}'.format(results.Eh[i])+', '
    text_area.insert(tk.INSERT, '\n   '+s)
    results.showPSI(figResults, canvasResults)


def doSeries():
    model.MultiParticleLoops = int(entryMPloops_str.get())
    s = model.seriesStr[combo_SeriesPara.current()]
    text_area.insert(tk.INSERT, '\n>> Series ' + s)
    s = model.QDshapeStr[model.QDpara[0]-1]+', ' + s
    FileName_data = osprocs.dataFolder + s + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+'.dat'
    paraIndex = combo_SeriesPara.current()
    seriesFrom = float(seriesFrom_str.get())
    seriesTo = float(seriesTo_str.get())
    seriesStep = float(seriesStep_str.get())
    OK = False
    paraVal = seriesFrom
    model.readInputVals()
    results.initResults()
    while OK == False:
        text_area.insert(tk.INSERT, '\n   ' + str(paraVal))
        window.update()
        if paraIndex in [0,1,2]:         # r1,r2,r3
            model.QDpara[paraIndex+1] = paraVal/1e9
        if paraIndex in [3,4]:          # Fz,Fx
            model.EFields[paraIndex-3] = paraVal
        if paraIndex in [5]:            # Bz
            model.BZfield[0] = paraVal

        model.setInputVals()
        if autoCell_bool.get(): model.autoCellSize()
        model.createModel()
        model.ShowPotential(figModel, canvasModel)
        computeSingle()
        WF.analyzeWF()
        results.E0peak()
        results.createResultsSingle()
        results.saveResults(FileName_data)
        paraVal += seriesStep
        OK = paraVal > seriesTo
    results.saveResults(FileName_data)
    text_area.insert(tk.INSERT, '\n   Series finished')


#############################################################################################
#### GUI

# Button materials
def btnMaterials_callback():
    frameMat.frameMaterials(window)
button_Materials = tk.Button(master=frameMenue, text = "Materials", width=15, height=1,
                           command = btnMaterials_callback)
button_Materials.pack(side='top',  anchor='w', padx = 5, pady = 10)

# Frame QD shape
frame_Shape = tk.Frame(frameMenue)
frame_Shape.pack(side='top', anchor='w', padx = (5,0), pady = 5)
label_Shape = tk.Label(master=frame_Shape, text = 'QD shape ')
label_Shape.pack(side='left', anchor='w')

def callback_combo_Shape(event):
    model.readInputVals()
    if combo_Shape.current() == 0:      # Ellipse
        model.QDpara[0] = 1
        model.QDpara[1] = 5e-9; model.QDpara[2] = 5e-9; model.QDpara[3] = 5e-9;
        label_ShapePara['text'] = model.QDshapeParaStr[0]
    if combo_Shape.current() == 1:      # Disk
        model.QDpara[0] = 2
        model.QDpara[1] = 10e-9; model.QDpara[2] = 5e-9; model.QDpara[3] = 0;
        label_ShapePara['text'] = model.QDshapeParaStr[1]
    if combo_Shape.current() == 2:      # Cone
        model.QDpara[0] = 3
        model.QDpara[1] = 10e-9; model.QDpara[2] = 10e-9; model.QDpara[3] = 0;
        label_ShapePara['text'] = model.QDshapeParaStr[2]
    if combo_Shape.current() == 3:      # Cone-Shell
        model.QDpara[0] = 4
        model.QDpara[1] = 35e-9; model.QDpara[2] = 19e-9; model.QDpara[3] = 13.6e-9;
        label_ShapePara['text'] = model.QDshapeParaStr[3]
    model.setInputVals()
    model.autoCellSize()

combo_Shape = ttk.Combobox(master=frame_Shape, width = 12)
combo_Shape['values'] = model.QDshapeStr
combo_Shape.pack(side='left', anchor='w', padx = (5,0))
combo_Shape.current(0)
combo_Shape.bind("<<ComboboxSelected>>", callback_combo_Shape)

# Frame parameters
label_ShapePara = tk.Label(master=frame_Shape, text = model.QDshapeParaStr[0])
label_ShapePara.pack(side='left', anchor='w')

framePara = tk.Frame(frameMenue)
framePara.pack(side='top', anchor='w')

def createParaInputFrames():
    strList = model.inputStr
    for widget in framePara.winfo_children():
        widget.destroy()
    for i in range(len(strList)):
        frame = tk.Frame(framePara)
        frame.pack(side='top', anchor='w', pady = (5,0))
        label_name = tk.Label(master=frame, width=10, text = strList[i][0])
        label_name.pack(side='left', anchor='w')
        for j in range(len(strList[i])-1):
            name = strList[i][j+1]
            label_name = tk.Label(master=frame, width=4, text = name)
            label_name.pack(side='left', anchor='w', padx = (2,0))
            entry_str = tk.StringVar()
            entry = tk.Entry(master=frame, width=6, textvariable=entry_str)
            entry.pack(side='left', anchor='w', padx = (2,0))
            entry_str.set(str(strList[i][j+1]))
            model.inputVals.append(entry_str)
    model.setInputVals()

createParaInputFrames()

# Frame point charges
framePC.createFramePC(frameMenue)

# Frame MultiParticles
frame_MultiParticles = tk.Frame(frameMenue)
frame_MultiParticles.pack(side='top', anchor='w', padx = (5,0), pady = 10)
label_MultiParticles = tk.Label(master=frame_MultiParticles, text = 'Particle')
label_MultiParticles.pack(side='left', anchor='w')
def callback_MultiParticles(event):
    model.MPpara[0] = combo_MultiParticles.current()
combo_MultiParticles = ttk.Combobox(master=frame_MultiParticles, width = 3)
combo_MultiParticles['values'] = model.MultiParticleStr
combo_MultiParticles.pack(side='left', anchor='w', padx = (5,0))
combo_MultiParticles.current(0)
combo_MultiParticles.bind("<<ComboboxSelected>>", callback_MultiParticles)
label_MultiParticleLoops = tk.Label(master=frame_MultiParticles, text = 'Self-consist. loops ')
label_MultiParticleLoops.pack(side='left', padx = (5,0))
entryMPloops_str = tk.StringVar()
entryMPloops_str.set('0')
entryMPloops = tk.Entry(master=frame_MultiParticles, width=2, textvariable=entryMPloops_str)
entryMPloops.pack(side='left')

# Frame WF threshold
frame_WFthreshold = tk.Frame(frameMenue)
frame_WFthreshold.pack(side='top', anchor='w', padx = (5,0), pady = 5)
label_WFthreshold = tk.Label(master=frame_WFthreshold, text = 'WF threshold')
label_WFthreshold.pack(side='left', anchor='w')
entryWFthreshold_str = tk.StringVar()
entryWFthreshold_str.set(str(model.WFflat_threshold))
entryWFthreshold = tk.Entry(master=frame_WFthreshold, width=5, textvariable=entryWFthreshold_str)
entryWFthreshold.pack(side='left', padx = (5,0))

# Frame buttons
frameExec = tk.Frame(frameMenue)
frameExec.pack(side='top', anchor='w', pady = (10,0))

# Button create model
def btnCreateModel_callback():
    model.readInputVals()
    framePC.setListPC()
    model.createModel()
    setSeriesPara()
    shapeStr = model.QDshapeStr[model.QDpara[0]-1]
    QDpara = model.QDpara
    text_area.insert(tk.INSERT, '\n>> Model created: '+shapeStr)
    text_area.insert(tk.INSERT, '\n   QD size   : ' +str(QDpara[1]*1e9)+', '+str(QDpara[2]*1e9)+', '+str(QDpara[3]*1e9) )
    text_area.insert(tk.INSERT, '\n   Cells (1D): ' +str(model.cellNum1D) )
    text_area.insert(tk.INSERT, '\n   Cell size : ' +'{0: >#6.4f}'.format(model.cellSize*1e9))
    text_area.insert(tk.INSERT, '\n   T: ' +str(model.T)+', xAl: '+str(model.xAl) )
    text_area.insert(tk.INSERT, '\n   Masses : ' +str(model.meQ)+', '+str(model.mhQ)+', '+str(model.meB)+', '+str(model.mhB)  )
    text_area.insert(tk.INSERT, '\n   E-field Fz: ' +str(model.EFields[0])+' Fx: '+str(model.EFields[1]) )
    text_area.insert(tk.INSERT, '\n   Bz-field: ' +str(model.BZfield[0])+' ('+str(int(model.BZfield[1]))+','+str(int(model.BZfield[2]))+')' )
    model.ShowPotential(figModel, canvasModel)


button_CreateModel = tk.Button(master=frameExec, text = "Create \nmodel", width=8, height=2,
                           command = btnCreateModel_callback)
button_CreateModel.pack(side='left')


# Button auto cell size
def btnAutoCell_callback():
    if model.modelCreated:
        model.autoCellSize()
        btnCreateModel_callback()
button_AutoCell = tk.Button(master=frameExec, text = "Auto \ncell", width=4, height=2,
                           command = btnAutoCell_callback)
button_AutoCell.pack(side='left', padx = 5)


# Button compute
def btnCompute_callback():
    if model.modelCreated: doCompute()
button_Compute = tk.Button(master=frameExec, text = "Compute", width=10, height=2,
                           command = btnCompute_callback)
button_Compute.pack(side='left', padx = 5)


# Button analyze WF
def btnAnalyzeWF_callback():
    if model.modelCreated:
        model.WFflat_threshold = float(entryWFthreshold_str.get())
        text_area.insert(tk.INSERT, '\n>> Analyze WF (threshold: '+str(model.WFflat_threshold)+')')
        window.update()
        runTime = WF.analyzeWF()
        text_area.insert(tk.INSERT, '\n   Run time: '+'{0: >#10.6f}'.format(runTime)+' s' )
        text_area.insert(tk.INSERT, '\n   deh [nm]       : '+'{0: >#10.8f}'.format(WF.deh*1e9))
        text_area.insert(tk.INSERT, '\n   Overlap int.   : '+'{0: >#10.8f}'.format(WF.overlap))
        text_area.insert(tk.INSERT, '\n   Lifetime [ns]  : '+'{0: >#10.8f}'.format(WF.tau))
        text_area.insert(tk.INSERT, '\n   Ceh [eV]       : '+'{0: >#10.8f}'.format(WF.Ceh))
        particle = model.MultiParticleStr[model.MPpara[0]]
        if particle in ['XX', 'X+', 'X-']:
            text_area.insert(tk.INSERT, '\n   Cee [eV]       : '+'{0: >#10.8f}'.format(WF.Cee))
            text_area.insert(tk.INSERT, '\n   Chh [eV]       : '+'{0: >#10.8f}'.format(WF.Chh))
        results.E0peak()
        text_area.insert(tk.INSERT, '\n   E [eV]         : '+'{0: >#10.8f}'.format(results.Ep))


button_AnalyzeWF = tk.Button(master=frameExec, text = "Analyze \nWF", width=10, height=2,
                           command = btnAnalyzeWF_callback)
button_AnalyzeWF.pack(side='left', padx = 5)


####################################
# Series

frameSeries = tk.Frame(frameMenue)
frameSeries.pack(side='top', anchor='w', pady = (10,0))

label_Series = tk.Label(master=frameSeries, text = 'Series: ')
label_Series.pack(side='top', anchor='w')

def callback_combo_SeriesPara(event):
    setSeriesPara()
combo_SeriesPara = ttk.Combobox(master=frameSeries, width = 12)
combo_SeriesPara['values'] = model.seriesStr
combo_SeriesPara.pack(side='top', anchor='w', padx = (5,0))
combo_SeriesPara.current(0)
combo_SeriesPara.bind("<<ComboboxSelected>>", callback_combo_SeriesPara)

frameSeriesRange = tk.Frame(frameSeries)
frameSeriesRange.pack(side='top', anchor='w', pady = (10,0))
label_SeriesFrom = tk.Label(master=frameSeriesRange, text = 'From: ')
label_SeriesFrom.pack(side='left', anchor='w')
seriesFrom_str = tk.StringVar(); seriesFrom_str.set('0')
entry_seriesFrom = tk.Entry(master=frameSeriesRange, width=6, textvariable=seriesFrom_str)
entry_seriesFrom.pack(side='left', anchor='w', padx = (2,0))
label_SeriesTo = tk.Label(master=frameSeriesRange, text = ' to: ')
label_SeriesTo.pack(side='left', anchor='w')
seriesTo_str = tk.StringVar(); seriesTo_str.set('1')
entry_seriesTo = tk.Entry(master=frameSeriesRange, width=6, textvariable=seriesTo_str)
entry_seriesTo.pack(side='left', anchor='w', padx = (2,0))
label_SeriesStep = tk.Label(master=frameSeriesRange, text = ' step: ')
label_SeriesStep.pack(side='left', anchor='w')
seriesStep_str = tk.StringVar(); seriesStep_str.set('1')
entry_seriesStep = tk.Entry(master=frameSeriesRange, width=6, textvariable=seriesStep_str)
entry_seriesStep.pack(side='left', anchor='w', padx = (2,0))
autoCell_bool = tk.BooleanVar()
Checkbutton_autoCell = tk.Checkbutton(master=frameSeriesRange, text="auto cell", variable=autoCell_bool)
Checkbutton_autoCell.pack(side='left', anchor='w', padx = (5,0))
autoCell_bool.set(True)

def setSeriesPara():
    paraIndex = combo_SeriesPara.current()
    if paraIndex in [0,1,2]:
        paraVal = model.QDpara[paraIndex+1]*1e9
    if paraIndex in [3,4]:
        paraVal = model.EFields[paraIndex-3]
    if paraIndex in [5]:
        paraVal = model.BZfield[paraIndex-5]
    seriesFrom_str.set(str(paraVal))
    seriesTo_str.set(str(paraVal+1))
    model.WFflat_threshold = float(entryWFthreshold_str.get())
setSeriesPara()


def btnGo_callback():
    doSeries()
button_Go = tk.Button(master=frameSeries, text = "Series go", width=12, height=2,
                           command = btnGo_callback)
button_Go.pack(side='top', anchor='w', padx = (5,0), pady = 10)


####################################
# Quit

def btnQuit_callback():
    window.quit()
    window.destroy()
button_Quit = tk.Button(master=frameSeries, text = "Quit", width=12, height=1,
                           command = btnQuit_callback)
button_Quit.pack(side='top', anchor='w', padx = (5,0), pady = 10)


####################################
# Results

# Textarea Memo and results
text_area = scrolledtext.ScrolledText(frameLog,
                                          wrap = tk.WORD,
                                          width = 70,
                                          height =45
                                          )
text_area.insert(tk.INSERT, '-- Log --')
text_area.pack(side='top', anchor='n', padx = (5,0))

# Button clear
def btnClear_callback():
    text_area.delete('1.0', tk.END)

button_Clear = tk.Button(master=frameLog, text = "Clear", width=12, height=1,
                           command = btnClear_callback)
button_Clear.pack(side='top', pady = 5)


####################################
# Main programm

osprocs.checkFolders()
window.mainloop()
