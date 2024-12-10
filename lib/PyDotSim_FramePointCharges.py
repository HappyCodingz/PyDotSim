"""

Frame for point charges for PyDotSim

by Christian Heyn

"""


import tkinter as tk
import lib.PyDotSim_Model as model


PClist = []     # List of point charges


# sets the defined point charges to the variable model.PointCharges
def setListPC():
    model.PointCharges = []
    for i in range(len(PClist)):
        PCx = float(PClist[i][0].get())*1e-9
        PCy = float(PClist[i][1].get())*1e-9
        PCz = float(PClist[i][2].get())*1e-9
        model.PointCharges.append([PCx,PCy,PCz])


# Creates a frame for the definition of point chardes
def createFramePC(frameMenue):
    frame_PointCharges = tk.Frame(frameMenue)
    frame_PointCharges.pack(side='top', anchor='w', padx = 5, pady = 5)
    frame_PCMenue = tk.Frame(frame_PointCharges)
    frame_PCMenue.pack(side='top', anchor='w', padx = 0, pady = 0)
    label_PC = tk.Label(master=frame_PCMenue, text = 'Point charges: ')
    label_PC.pack(side='left', anchor='w')

    def btnPCadd_callback():
        global PClist, frame_PClist
        if len(PClist) == 0:
            frame_PClist = tk.Frame(frame_PointCharges)
            frame_PClist.pack(side='top', anchor='w', padx = 0, pady = 0)
        frame_singlePC = tk.Frame(frame_PClist)
        frame_singlePC.pack(side='top', anchor='w', padx = 0, pady = 0)
        label_singlePC = tk.Label(master=frame_singlePC)
        label_singlePC['text'] = str(len(PClist)+1)+') x,y,z [nm] :'
        label_singlePC.pack(side='left')
        PCx_str = tk.StringVar(); PCx_str.set('0')
        PCy_str = tk.StringVar(); PCy_str.set('0')
        PCz_str = tk.StringVar(); PCz_str.set('0')
        entryPCx = tk.Entry(master=frame_singlePC, width=5, textvariable=PCx_str)
        entryPCx.pack(side='left', padx = 5)
        entryPCy = tk.Entry(master=frame_singlePC, width=5, textvariable=PCy_str)
        entryPCy.pack(side='left', padx = 5)
        entryPCz = tk.Entry(master=frame_singlePC, width=5, textvariable=PCz_str)
        entryPCz.pack(side='left', padx = 5)
        PClist.append([entryPCx,entryPCy,entryPCz])
    button_PCadd = tk.Button(master=frame_PCMenue, text = "Add", width=7, height=1,
                           command = btnPCadd_callback)
    button_PCadd.pack(side='left', padx = 5)

    def btnPCclear_callback():
        global PClist
        if len(PClist) > 0:
            PClist = []
            model.PointCharges = []
            frame_PClist.destroy()
    button_PCclear = tk.Button(master=frame_PCMenue, text = "Clear all", width=7, height=1,
                           command = btnPCclear_callback)
    button_PCclear.pack(side='left',  padx = 15)


