"""

Frame for materials for PyDotSim

by Christian Heyn

"""


import tkinter as tk
from tkinter import ttk
import lib.PyDotSim_Model as model
import lib.PyDotSim_Constants as cnst


windowIsOpen =          False           # Flag if window is already open
windowLeft =            350             # Left position setup window
windowTop =             35              # Top position setup window



def frameMaterialsIsOpen():
    return windowIsOpen


# New window with GUI for frame materials
def frameMaterials(window):
    global windowIsOpen

    if windowIsOpen: return
    newWindow = tk.Toplevel(window)
    newWindow.geometry('260x320+'+str(windowLeft)+'+'+str(windowTop))
    windowIsOpen = True

    # Label
    label_01 = tk.Label(master=newWindow, text = 'Materials')
    label_01.pack(side='top', pady = (5,0))

    # Frame Materials
    frame_Materials = tk.Frame(newWindow)
    frame_Materials.pack(side='top', anchor='w', padx = (5,0), pady = 15)
    label_Materials = tk.Label(master=frame_Materials, text = 'Materials')
    label_Materials.pack(side='left', anchor='w')
    def callback_Materials(event):
        pass
    combo_Materials = ttk.Combobox(master=frame_Materials, width = 20)
    combo_Materials['values'] = ['GaAs-AlGaAs']
    combo_Materials.pack(side='left', anchor='w', padx = (5,0))
    combo_Materials.current(0)
    combo_Materials.bind("<<ComboboxSelected>>", callback_Materials)

    x_str = tk.StringVar()
    T_str = tk.StringVar()
    dEc_str = tk.StringVar()
    dEv_str = tk.StringVar()
    meQ_str = tk.StringVar()
    mhQ_str = tk.StringVar()
    meB_str = tk.StringVar()
    mhB_str = tk.StringVar()

    # Frame x, T
    frame_xT = tk.Frame(newWindow)
    frame_xT.pack(side='top', anchor='w', padx = (5,0), pady = 5)
    label_x = tk.Label(master=frame_xT, text = 'x')
    label_x.pack(side='left')
    x_str.set(str(model.xAl))
    entry_x = tk.Entry(master=frame_xT, width=5, textvariable=x_str)
    entry_x.pack(side='left', padx = 5)
    label_T = tk.Label(master=frame_xT, text = 'T [K]')
    label_T.pack(side='left', padx = 0)
    T_str.set(str(model.T))
    entry_T = tk.Entry(master=frame_xT, width=5, textvariable=T_str)
    entry_T.pack(side='left', padx = 15)
    # Button Apply
    def btnApply_callback():
        model.xAl = float(x_str.get())
        model.T =   float(T_str.get())
        dEc, dEv = cnst.bandgapDiscontinuity(model.T, model.xAl)
        dEc_str.set(str(round(dEc,4)))
        dEv_str.set(str(round(dEv,4)))
        meQ, mhQ =  cnst.AlGaAs_effMass(0)
        meQ_str.set(str(round(meQ,4)))
        mhQ_str.set(str(round(mhQ,4)))
        meB, mhB =  cnst.AlGaAs_effMass(model.xAl)
        meB_str.set(str(round(meB,4)))
        mhB_str.set(str(round(mhB,4)))


    button_Apply = tk.Button(master=frame_xT, text = "Apply", width=6,
                            command =  btnApply_callback)
    button_Apply.pack(side='left', anchor='w', padx = 10)

    # Frame dEc, dEv
    dEc, dEv = cnst.bandgapDiscontinuity(model.T, model.xAl)
    frame_dE = tk.Frame(newWindow)
    frame_dE.pack(side='top', anchor='w', padx = (5,0), pady = 5)
    label_dEc = tk.Label(master=frame_dE, text = 'dEc [eV]')
    label_dEc.pack(side='left')
    dEc_str.set(str(round(dEc,4)))
    entry_dEc = tk.Entry(master=frame_dE, width=7, textvariable=dEc_str)
    entry_dEc.pack(side='left', padx = 5)
    label_dEv = tk.Label(master=frame_dE, text = 'dEv [eV]')
    label_dEv.pack(side='left', padx = 10)
    dEv_str.set(str(round(dEv,4)))
    entry_dEv = tk.Entry(master=frame_dE, width=7, textvariable=dEv_str)
    entry_dEv.pack(side='left', padx = 0)


    # Frame eff. mass QD
    meQ, mhQ =  cnst.AlGaAs_effMass(0)
    frame_mQD = tk.Frame(newWindow)
    frame_mQD.pack(side='top', anchor='w', padx = (5,0), pady = 5)
    label_meQ = tk.Label(master=frame_mQD, text = 'meQD   ')
    label_meQ.pack(side='left')
    meQ_str.set(str(round(meQ,4)))
    entry_meQ = tk.Entry(master=frame_mQD, width=7, textvariable=meQ_str)
    entry_meQ.pack(side='left', padx = 5)
    label_mhQ = tk.Label(master=frame_mQD, text = 'mhQD   ')
    label_mhQ.pack(side='left', padx = 10)
    mhQ_str.set(str(round(mhQ,4)))
    entry_mhQ = tk.Entry(master=frame_mQD, width=7, textvariable=mhQ_str)
    entry_mhQ.pack(side='left', padx = 0)


    # Frame eff. mass Barrier
    meB, mhB =  cnst.AlGaAs_effMass(model.xAl)
    frame_mBar = tk.Frame(newWindow)
    frame_mBar.pack(side='top', anchor='w', padx = (5,0), pady = 5)
    label_meB = tk.Label(master=frame_mBar, text = 'meBar   ')
    label_meB.pack(side='left')
    meB_str.set(str(round(meB,4)))
    entry_meB = tk.Entry(master=frame_mBar, width=7, textvariable=meB_str)
    entry_meB.pack(side='left', padx = 5)
    label_mhB = tk.Label(master=frame_mBar, text = 'mhBar   ')
    label_mhB.pack(side='left', padx = 10)
    mhB_str.set(str(round(mhB,4)))
    entry_mhB = tk.Entry(master=frame_mBar, width=7, textvariable=mhB_str)
    entry_mhB.pack(side='left', padx = 0)


    frame_smoothBorders = tk.Frame(newWindow)
    frame_smoothBorders.pack(side='top', anchor='w', padx=(5,5), pady=(5,0))
    withSmoothBorders = tk.BooleanVar()
    withSmoothBorders.set(model.smoothMaterialBorders)
    def _callback():
        model.smoothMaterialBorders = withSmoothBorders.get()
    Checkbutton_withSmoothBorders = tk.Checkbutton(master=frame_smoothBorders,text="With smooth borders", variable=withSmoothBorders, command=_callback)
    Checkbutton_withSmoothBorders.pack(side='top', anchor='w', padx=5)

    # Button Quit
    def btnQuit_callback():
        global windowIsOpen
        model.xAl = float(x_str.get())
        model.T =   float(T_str.get())
        model.meQ = float(meQ_str.get())
        model.mhQ = float(mhQ_str.get())
        model.meB = float(meB_str.get())
        model.mhB = float(mhB_str.get())
        newWindow.destroy()
        windowIsOpen = False
    button_Quit = tk.Button(master=newWindow, text = "Quit", width=12,
                            command =  btnQuit_callback)
    button_Quit.pack(side='top', anchor='w', padx = 5, pady = (15,0))





