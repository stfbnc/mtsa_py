# -*- coding: iso-8859-15 -*-

import MTsA_fun_app as MF
import Outputs_win_app as OW
import Tkinter as tk
import ttk
import tkMessageBox
from tkFileDialog import askopenfilename
#import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import os
### the window shows up, otherwise you need to click on icon in the dock (only for app, not if launched form command line)
#os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
### more general
#script = 'tell application "System Events" to set frontmost of the first process whose unix id is {pid} to true'.format(pid=os.getpid())
#os.system("/usr/bin/osascript -e '{script}'".format(script=script))

### TS_PLOT
def ts_plot():

    Frame_fig = tk.Frame(root)
    Frame_fig.grid(row = 0,column = 1,padx = int(root.winfo_screenwidth() * 0.01))

    fig_ts = Figure(figsize = (int(root.winfo_screenwidth() * 0.01),int(root.winfo_screenheight() * 0.008)))
    a = fig_ts.add_subplot(111)
    a.plot(-1,-1,'r')
    a.set_ylabel('$X_t$',fontsize = 15)
    a.set_xlabel('t',fontsize = 15)
    a.set_title('Time series',fontsize = 15)
    a.set_xlim((0,1000))
    a.set_ylim((0,1000))
    fig_ts.tight_layout()
    canvas = FigureCanvasTkAgg(fig_ts,master = Frame_fig)
    canvas.get_tk_widget().grid(row = 0,column = 0)
    canvas.draw()

### FILE_PLOT
def file_plot(file_name,data_col,time_col,Delta_t,typeoffit):

    pn,time,useless_var,t = MF.load_file(file_name,data_col,time_col,Delta_t)
    
    if typeoffit != '':
        useless_pn,fitted_curve = MF.trend_detrend(pn,t,typeoffit)

    ticks = np.arange(0,len(t),len(t) / 10,dtype = int)
    t = np.array(t,dtype = int)
    time = np.array(time,dtype = str)
    ticks_vec = t[ticks]
    time_label = time[ticks]

    Frame_fig = tk.Frame(root)
    Frame_fig.grid(row = 0,column = 1,padx = int(root.winfo_screenwidth() * 0.01))
    fig_ts = Figure(figsize = (int(root.winfo_screenwidth() * 0.01),int(root.winfo_screenheight() * 0.008)))
    a = fig_ts.add_subplot(111)
    a.plot(t,pn,'r')
    if typeoffit != '' and not np.isscalar(fitted_curve):
        a.plot(t,fitted_curve,'k')
        a.plot(t,pn - fitted_curve,'b')
    a.set_ylabel('$X_t$',fontsize = 15)
    a.set_xlabel('t',fontsize = 15)
    a.set_title('Time series',fontsize = 15)
    a.set_xlim((t[0],t[-1]))
    a.set_xticks(ticks_vec)
    a.set_xticklabels(time_label)
    fig_ts.tight_layout()
    canvas = FigureCanvasTkAgg(fig_ts,master = Frame_fig)
    canvas.get_tk_widget().grid(row = 0,column = 0)
    canvas.draw()

### BUTTONS_STATE
def buttons_state(state_string):
    
    file_button.config(state = state_string)
    plot_button.config(state = state_string)
    for i in range(len(entry_1)):
        entry_1[i].config(state = state_string)
    for i in range(len(entry_2)):
        entry_2[i].config(state = state_string)
    unit.config(state = state_string)
    tf.config(state = state_string)
    ok_button.config(state = state_string)
    for i in range(len(entry_3)):
        entry_3[i].config(state = state_string)

### CHOOSE_FILE
def choose_file():
    
    global file_name

    root.update()
    file_name = askopenfilename()
    if file_name != '':
        buttons_state('normal')

### PLOT_BUTTON
def PLOT_button():

    if 'file_name' not in globals():
        tkMessageBox.showerror("","Please select a file!",icon = 'error')
    elif file_name == '':
        tkMessageBox.showerror("","Please select a file!",icon = 'error')
    else:
        entries_check = 0
        try:
            time_col = int(entry_1[0].get())
            if time_col <= 0:
                entries_check = 1
                tkMessageBox.showerror("","Time column must be a positive number greater than zero! Please reinsert value",icon = 'error')
            data_col = int(entry_1[1].get())
            if data_col <= 0:
                entries_check = 1
                tkMessageBox.showerror("","Data column must be a positive number greater than zero! Please reinsert value",icon = 'error')
            Delta_t = float(entry_1[2].get())
            if Delta_t <= 0:
                entries_check = 1
                tkMessageBox.showerror("","Sampling time must be a positive number greater than zero! Please reinsert value",icon = 'error')
            typeoffit = var_fit.get()
            if typeoffit not in fit_type:
                typeoffit = ''
            if entries_check == 0:
                file_plot(file_name,data_col,time_col,Delta_t,typeoffit)
        except:
            tkMessageBox.showerror("","Not all variables inserted!",icon = 'error')

### OK_BUTTON
def OK_button():
    
    global time_col
    global data_col
    global Delta_t
    global typeoffit

    def out_win():
        top = tk.Toplevel()
        top.geometry("680x350")
        top.wm_title("Outputs")
        top.resizable(width = False,height = False)
        
        top_frame_1 = tk.Frame(top)
        top_frame_1.grid(row = 0,column = 0)
        tk.Label(top_frame_1,text = "Periodicities",font = "Verdana 20 bold").grid(row = 0,column = 0)
        
        top_frame_2 = tk.Frame(top)
        top_frame_2.grid(row = 1,column = 0)
        tk.Label(top_frame_2,text = "Filter round     ",font = "Verdana 13 bold").grid(row = 0,column = 0,sticky = tk.W)
        tk.Label(top_frame_2,text = "Frequency          ",font = "Verdana 13 bold").grid(row = 0,column = 1,sticky = tk.W)
        tk.Label(top_frame_2,text = "Period          ",font = "Verdana 13 bold").grid(row = 0,column = 2,sticky = tk.W)
        tk.Label(top_frame_2,text = "Percentage",font = "Verdana 13 bold").grid(row = 0,column = 3,sticky = tk.W)
        tk.Label(top_frame_2,text = "").grid(row = 1,column = 0)
        tk.Label(top_frame_2,text = "[" + time_units + "^-1]",font = "Verdana 13 bold").grid(row = 1,column = 1,sticky = tk.W)
        tk.Label(top_frame_2,text = "[" + time_units + "]",font = "Verdana 13 bold").grid(row = 1,column = 2,sticky = tk.W)
        tk.Label(top_frame_2,text = "").grid(row = 1,column = 3,sticky = tk.W)
        freq_list = tk.Text(top_frame_2,font = "Verdana 13",width = 48,height = 6,borderwidth = 3,relief = tk.RIDGE)
        for i in range(len(percent)):
            freq_list.insert(tk.END,"%7d%25.8f%16.2f%21.2f\n" % (num_filter[i],1.0 / freq_fig[i],freq_fig[i],percent[i]))
        freq_list.grid(row = 2,column = 0,columnspan = 4)
        freq_list.config(state = 'disabled')
        scrollbar = tk.Scrollbar(top_frame_2,command = freq_list.yview)
        scrollbar.grid(row = 2,column = 4,sticky = 'ns')
        freq_list.config(yscrollcommand = scrollbar.set)

        path_only,file_only = os.path.split(file_name)
        top_frame_5 = tk.Frame(top)
        top_frame_5.grid(row = 3,column = 0)
        tk.Button(top_frame_5,text = "Save in a file",font = "Verdana 13 bold",command = lambda:OW.save_freq(file_only,num_filter,freq_fig,percent)).grid(row = 0,column = 0)
        
        top_frame_3 = tk.Frame(top)
        top_frame_3.grid(row = 0,column = 1,padx = 10)
        tk.Label(top_frame_3,text = "Graphic outputs",font = "Verdana 20 bold").grid(row = 0,column = 0)
        
        top_frame_4 = tk.Frame(top)
        top_frame_4.grid(row = 1,column = 1,padx = 10)
        tk.Label(top_frame_4,text = "").grid(row = 0,column = 0)
        tk.Button(top_frame_4,text = "Time Series",font = "Verdana 13 bold",width = 14,command = lambda:OW.detrend_plot(root,pn_plot,fitted_curve,time)).grid(row = 1,column = 0)
        tk.Button(top_frame_4,text = "Spectrum",font = "Verdana 13 bold",width = 14,command = lambda:OW.spectrum_plot(root,PNT_plot,freq,pth)).grid(row = 2,column = 0)
        tk.Button(top_frame_4,text = "Final spectrum",font = "Verdana 13 bold",width = 14,command = lambda:OW.spectrum_plot(root,PNT_single,freq_single,0)).grid(row = 3,column = 0)
        tk.Button(top_frame_4,text = "Residuals",font = "Verdana 13 bold",width = 14,command = lambda:OW.res_plot(root,pn,s_res,pn_norm,time)).grid(row = 4,column = 0)
        tk.Button(top_frame_4,text = "DFA",font = "Verdana 13 bold",width = 14,command = lambda:OW.dfa_plot(root,s,F,log_fit,H_mono)).grid(row = 5,column = 0)
        tk.Button(top_frame_4,text = "MFDFA",font = "Verdana 13 bold",width = 14,command = lambda:OW.mdfa_plot(root,s2,F2,MDFA_fit,q,H,H_mono,alpha,sing_spec)).grid(row = 6,column = 0)
        tk.Button(top_frame_4,text = "H(t)",font = "Verdana 13 bold",width = 14,command = lambda:OW.MFDFA2_plot(root,Ht_plot,Htbin,Ph,fit_gauss,mu,sigma)).grid(row = 7,column = 0)
    
    entries_check = 0
    try:
        time_col = int(entry_1[0].get())
        if time_col <= 0:
            entries_check = 1
            tkMessageBox.showerror("","Time column must be a positive number greater than zero! Please reinsert value",icon = 'error')
        data_col = int(entry_1[1].get())
        if data_col <= 0:
            entries_check = 1
            tkMessageBox.showerror("","Data column must be a positive number greater than zero! Please reinsert value",icon = 'error')
        Delta_t = float(entry_1[2].get())
        if Delta_t <= 0:
            entries_check = 1
            tkMessageBox.showerror("","Sampling time must be a positive number greater than zero! Please reinsert value",icon = 'error')
        year_in = int(entry_2[0].get())
        year_fin = int(entry_2[1].get())
        if year_in <= 0 or year_fin <= 0 or year_fin <= year_in:
            entries_check = 1
            tkMessageBox.showerror("","Invalid years! Please reinsert value",icon = 'error')
        ofac = int(entry_2[2].get())
        if ofac <= 0:
            entries_check = 1
            tkMessageBox.showerror("","Oversampling factor must be a positive number greater than zero! Please reinsert value",icon = 'error')
        time_units = var_unit.get()
        units_choices = {'seconds':31536000,'minutes':525600,'hours':8760,'days':365,'weeks':52,'months':12}
        units_converter = units_choices.get(time_units)
        if units_converter not in units_choices.values():
            entries_check = 1
            tkMessageBox.showerror("","Invalid time units! Please select one",icon = 'error')
        typeoffit = var_fit.get()
        if typeoffit not in fit_type:
            entries_check = 1
            tkMessageBox.showerror("","Invalid fit function! Please select one",icon = 'error')
        scale_min = int(entry_3[0].get())
        if scale_min < 3:
            entries_check = 1
            tkMessageBox.showerror("","Minimum scale for DFA should be at least 3! Please reinsert value",icon = 'error')
        scale_MFDFA = [int(i) for i in entry_3[1].get().split(',')]
        for i in scale_MFDFA:
            if i < 3:
                entries_check = 1
                tkMessageBox.showerror("","Scales for local Hurst exponent can not be less than 3! Please reinsert value",icon = 'error')
    except:
        entries_check = 1
        tkMessageBox.showerror("","Not all variables inserted!",icon = 'error')

    if entries_check == 0:
        buttons_state('disabled')
        ### replot time series in case plot button wasn't pressed
        file_plot(file_name,data_col,time_col,Delta_t,typeoffit)
        
        top_bar = tk.Toplevel()
        top_bar.geometry("400x120")
        top_bar.wm_title("")
        top_bar.resizable(width = False,height = False)
        tk.Label(top_bar,text = "Analysing...",font = "Verdana 20 bold").grid(row = 0,column = 0,pady = 10)
        style_bar = ttk.Style()
        style_bar.theme_use("classic")
        style_bar.configure("TProgressbar",thickness = 25)
        progress_var = tk.DoubleVar()
        progressbar = ttk.Progressbar(top_bar,variable = progress_var,maximum = 8,length = 350,style = "TProgressbar")
        progressbar.grid(row = 1,column = 0,pady = 10,padx = 25)
        
        root.update()
        ### read file
        pn,time,nan_pos_pn,t = MF.load_file(file_name,data_col,time_col,Delta_t)
        progress_var.set(0)
        top_bar.update_idletasks()
        ### detrending
        pn,fitted_curve = MF.trend_detrend(pn,t,typeoffit)
        pn_plot = pn.copy()
        progress_var.set(1)
        top_bar.update_idletasks()
        ### lomb spectrum
        PNT,freq,pth = MF.lomb_scargle(pn,t,ofac,Delta_t)
        PNT_plot = PNT.copy()
        progress_var.set(2)
        top_bar.update_idletasks()
        ### finding all frequencies
        freq_fig,percent,num_filter = MF.peaks_filter_all(len(freq),PNT,freq,pth,units_converter,year_in,year_fin,t,nan_pos_pn,Delta_t,pn,ofac)
        progress_var.set(3)
        top_bar.update_idletasks()
        ### filter peaks
        pn = MF.peaks_filter(1,PNT,freq,pth,units_converter,year_in,year_fin,t,nan_pos_pn,Delta_t,pn,ofac)
        progress_var.set(4)
        top_bar.update_idletasks()
        ### spectrum of residuals and outliers
        s_res,pn_norm,freq_single,PNT_single = MF.res_lomb(pn,Delta_t,t)
        progress_var.set(5)
        top_bar.update_idletasks()
        ### detrended fluctuation analysis
        s,F,log_fit,H_mono = MF.dfa(pn,scale_min,1)
        progress_var.set(6)
        top_bar.update_idletasks()
        ### multifractal detrended fluctuation analysis
        s2,F2,MDFA_fit,q,H,alpha,sing_spec = MF.mdfa(H_mono,pn,scale_min,3,1)
        progress_var.set(7)
        top_bar.update_idletasks()
        ### local hurst exponent
        Ht_plot,Htbin,Ph,fit_gauss,mu,sigma = MF.MFDFA2(pn,scale_MFDFA,1)
        progress_var.set(8)
        top_bar.update_idletasks()
        top_bar.destroy()
        buttons_state('normal')
        out_win()

################################ MAIN ###################################################
root = tk.Tk()
root.wm_title("MTsA")
root.geometry("%dx%d" % (int(root.winfo_screenwidth() * 0.93),int(root.winfo_screenheight() * 0.67)))
root.resizable(width = False,height = False)

Frame_entries = tk.Frame(root)
Frame_entries.grid(row = 0,column = 0)

file_button = tk.Button(Frame_entries,text = "Choose file",font = "Verdana 10 bold",command = choose_file)
file_button.grid(row = 0,column = 0,sticky = tk.W)
tk.Label(Frame_entries,text = "").grid(row = 1,column = 0)

labels_1 = [" Times column"," Data column"," Sampling time"]
entry_1 = []
for i in range(2,2 * len(labels_1) + 2,2):
    tk.Label(Frame_entries,text = labels_1[(i - 2) / 2],font = "Verdana 10 bold").grid(row = i,column = 0,sticky = tk.W)
    entry_1.append(tk.Entry(Frame_entries,width = 4))
    entry_1[-1].grid(row = i,column = 1)
    entry_1[-1].config(state = 'disabled')

labels_2 = [" Initial year"," Final year"," Oversampling factor"]
entry_2 = []
for i in range(2 * len(labels_1) + 2,2 * len(labels_1) + 2 + 2 * len(labels_2),2):
    tk.Label(Frame_entries,text = labels_2[(i - (2 * len(labels_1) + 2)) / 2],font = "Verdana 10 bold").grid(row = i,column = 0,sticky = tk.W)
    entry_2.append(tk.Entry(Frame_entries,width = 4))
    entry_2[-1].grid(row = i,column = 1)
    entry_2[-1].config(state = 'disabled')

for i in range(1,2 * len(labels_1) + 3 + 2 * len(labels_2),2):
    tk.Label(Frame_entries,text = "").grid(row = i,column = 0)

unit_label = tk.Label(Frame_entries,text = " Time units",font = "Verdana 10 bold")
unit_label.grid(row = 18,column = 0,sticky = tk.W)
var_unit = tk.StringVar(root)
var_unit.set("")
unit = tk.OptionMenu(Frame_entries,var_unit,"seconds","minutes","hours","days","weeks","months")
unit.config(width = 5)
unit.grid(row = 18,column = 1)
unit.config(state = 'disabled')

tk.Label(Frame_entries,text = "").grid(row = 19,column = 0)

tf_label = tk.Label(Frame_entries,text = " Type of fit",font = "Verdana 10 bold")
tf_label.grid(row = 20,column = 0,sticky = tk.W)
var_fit = tk.StringVar(Frame_entries)
var_fit.set("")
fit_type = ['none','poly1','poly2','poly3','poly4','poly5','exp']
tf = tk.OptionMenu(Frame_entries,var_fit,*fit_type)
tf.config(width = 5)
tf.grid(row = 20,column = 1)
tf.config(state = 'disabled')

tk.Label(Frame_entries,text = "").grid(row = 21,column = 0)

labels_3 = [" Minimum scale"," H(t) scales"]
entry_3 = []
for i in range(22,25,2):
    tk.Label(Frame_entries,text = labels_3[(i - 22) / 2],font = "Verdana 10 bold").grid(row = i,column = 0,sticky = tk.W)
    entry_3.append(tk.Entry(Frame_entries,width = 4))
    entry_3[-1].grid(row = i,column = 1)
    entry_3[-1].config(state = 'disabled')

for i in range(23,26,2):
    tk.Label(Frame_entries,text = "").grid(row = i,column = 0)

plot_button = tk.Button(Frame_entries,text = "Plot!",font = "Verdana 10 bold",command = PLOT_button)
plot_button.grid(row = 26,column = 0,sticky = tk.W)
plot_button.config(state = 'disabled')

ok_button = tk.Button(Frame_entries,text = "Analyse!",font = "Verdana 10 bold",command = OK_button)
ok_button.grid(row = 26,column = 1)
ok_button.config(state = 'disabled')

ts_plot()

root.mainloop()
