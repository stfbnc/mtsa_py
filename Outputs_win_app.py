import sys
if sys.version_info[0] == 2:
    import Tkinter as tk
    from tkFileDialog import askdirectory
else:
    import tkinter as tk
    from tkinter.filedialog import askdirectory
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

### set size for figures
x_size = 7
y_size = 5

#############################################################################
############################## SAVE_FREQ ####################################
#############################################################################
def save_freq(arg0,arg2,arg3,arg4):
    
    path_file = askdirectory()
    path_file = os.path.join(path_file,'freq_perc.txt')
    f = open(path_file,'w')
    f.write('Frequencies for file %s\n\n' % (arg0))
    f.write('Filter_round Frequency Period Percentage\n\n')
    for i in range(len(arg2)):
        f.write('%d %.8f %.2f %.2f\n' % (arg2[i],1.0 / arg3[i],arg3[i],arg4[i]))
    f.close()
#############################################################################
#############################################################################
############################## DETREND_PLOT #################################
#############################################################################
def detrend_plot(main_win,arg0,arg1,arg2):
    
    ticks = np.arange(0,len(arg0),len(arg0) / 7,dtype = int)
    t = np.arange(1,len(arg0) + 1,dtype = int)
    time = np.array(arg2,dtype = str)
    ticks_vec = t[ticks]
    time_label = time[ticks]
    
    def save_fig():
        path_tot = askdirectory()
        plt.rc('text',usetex = True)
        plt.rc('font',family = 'serif')
        plt.plot(t,arg0 + arg1,'r',label = 'original')
        if not np.isscalar(arg1):
            plt.plot(t,arg1,'k',label = 'trend')
            plt.plot(t,arg0,'b',label = 'detrended')
        plt.legend(loc = 0)
        plt.xlim(float(entries[0].get()),float(entries[1].get()))
        plt.ylim(float(entries[2].get()),float(entries[3].get()))
        plt.xlabel(entries[4].get())
        plt.ylabel(entries[5].get())
        plt.title(entries[6].get())
        plt.xticks(ticks,time_label)
        plt.margins(0.2)
        plt.subplots_adjust(bottom = 0.2)
        plt.savefig(os.path.join(path_tot,'ts.pdf'))
        plt.close()

    def screen_fig():
        fig_ts = Figure(figsize = (x_size,y_size))
        a = fig_ts.add_subplot(111)
        a.plot(t,arg0 + arg1,'r',label = 'original')
        if not np.isscalar(arg1):
            a.plot(t,arg1,'k',label = 'trend')
            a.plot(t,arg0,'b',label = 'detrended')
        a.legend(loc = 0)
        a.set_xlim(float(entries[0].get()),float(entries[1].get()))
        a.set_ylim(float(entries[2].get()),float(entries[3].get()))
        a.set_xlabel(entries[4].get(),fontsize = 15)
        a.set_ylabel(entries[5].get(),fontsize = 15)
        a.set_title(entries[6].get(),fontsize = 15)
        a.set_xticks(ticks_vec)
        a.set_xticklabels(time_label)
        fig_ts.tight_layout()
        canvas = FigureCanvasTkAgg(fig_ts,master = frame_1)
        canvas.get_tk_widget().grid(row = 0,column = 0)
        canvas.draw()

    def reset_fig():
        for i in range(len(entries)):
            entries[i].delete(0,tk.END)
            entries[i].insert(0,values[i])
        screen_fig()
    
    top = tk.Toplevel(main_win)
    top.geometry("%dx%d" % (int(main_win.winfo_screenwidth() * 0.93 * 0.85),
                            int(main_win.winfo_screenheight() * 0.65)))
    top.wm_title("Time Series")
    top.resizable(width = False,height = False)
    
    frame_1 = tk.Frame(top)
    frame_1.grid(row = 0,column = 0)
    
    frame_2 = tk.Frame(top)
    frame_2.grid(row = 0,column = 1)
    names = ["X Limit (left)","X Limit (right)","Y Limit (bottom)","Y Limit (top)","X Label","Y Label","Title"]
    if not np.isscalar(arg1):
        values = [t[0],t[-1],np.min([np.min(arg0[~np.isnan(arg0)]),np.min(arg1[~np.isnan(arg1)]),
                                     np.min(arg0[~np.isnan(arg0)] + arg1[~np.isnan(arg1)])]) - 1.0,
                  np.max([np.max(arg0[~np.isnan(arg0)]),np.max(arg1[~np.isnan(arg1)]),
                          np.max(arg0[~np.isnan(arg0)] + arg1[~np.isnan(arg1)])]) + 1.0,'t','$X_t$',
                  'Time Series']
    else:
        values = [t[0],t[-1],np.min(arg0[~np.isnan(arg0)]) - 1.0,np.max(arg0[~np.isnan(arg0)]) + 1.0,
                  't','$X_t$','Time Series']
    entries = []
    for i in range(len(names)):
        tk.Label(frame_2,text = names[i],font = "Verdana 13 bold").grid(row = 2 * i,column = 0,
                                                                        padx = int(main_win.winfo_screenwidth() * 0.01))
        entries.append(tk.Entry(frame_2,width = 18))
        entries[-1].insert(0,values[i])
        entries[-1].grid(row = 2 * i,column = 1)
    for i in range(len(names)):
        tk.Label(frame_2,text = "").grid(row = 2 * i + 1,column = 0)
    screen_fig()
    tk.Button(frame_2,text = "Replot",font = "Verdana 13 bold",command = screen_fig).grid(row = 2 * len(names),column = 0)
    tk.Button(frame_2,text = "Save",font = "Verdana 13 bold",command = save_fig).grid(row = 2 * len(names),column = 1)
    tk.Label(frame_2,text = "").grid(row = 2 * len(names) + 1,column = 0)
    tk.Button(frame_2,text = "Reset",font = "Verdana 13 bold",command = reset_fig).grid(row = 2 * len(names) + 2,column = 0)
#############################################################################
#############################################################################
############################## SPECTRUM_PLOT ################################
#############################################################################
def spectrum_plot(main_win,arg0,arg1,arg2):
    
    def save_fig():
        path_tot = askdirectory()
        plt.rc('text',usetex = True)
        plt.rc('font',family = 'serif')
        plt.plot(arg1,arg0,'b')
        if arg2 != 0:
            plt.plot((arg1[0],arg1[-1]),(arg2,arg2),'r')
        plt.xlabel(entries[0].get())
        plt.ylabel(entries[1].get())
        plt.xlim(float(entries[2].get()),float(entries[3].get()))
        plt.ylim(float(entries[4].get()),float(entries[5].get()))
        plt.title(entries[6].get())
        plt.savefig(os.path.join(path_tot,'spectrum_in.pdf'))
        plt.close()
    
    def screen_fig():
        fig_ts = Figure(figsize = (x_size,y_size))
        a = fig_ts.add_subplot(111)
        a.plot(arg1,arg0,'b')
        if arg2 != 0:
            a.plot((arg1[0],arg1[-1]),(arg2,arg2),'r')
        a.set_xlabel(entries[0].get(),fontsize = 15)
        a.set_ylabel(entries[1].get(),fontsize = 15)
        a.set_xlim(float(entries[2].get()),float(entries[3].get()))
        a.set_ylim(float(entries[4].get()),float(entries[5].get()))
        a.set_title(entries[6].get(),fontsize = 15)
        fig_ts.tight_layout()
        canvas = FigureCanvasTkAgg(fig_ts,master = frame_1)
        canvas.get_tk_widget().grid(row = 0,column = 0)
        canvas.draw()
    
    def reset_fig():
        for i in range(len(entries)):
            entries[i].delete(0,tk.END)
            entries[i].insert(0,values[i])
        screen_fig()
    
    top = tk.Toplevel(main_win)
    top.geometry("%dx%d" % (int(main_win.winfo_screenwidth() * 0.93 * 0.85),
                            int(main_win.winfo_screenheight() * 0.65)))
    if arg2 != 0:
        top.wm_title("Spectrum")
    else:
        top.wm_title("Spectrum of residuals")
    top.resizable(width = False,height = False)

    frame_1 = tk.Frame(top)
    frame_1.grid(row = 0,column = 0)
    
    frame_2 = tk.Frame(top)
    frame_2.grid(row = 0,column = 1)
    names = ["X Label","Y Label","X Limit (left)","X Limit (right)","Y Limit (bottom)","Y Limit (top)","Title"]
    if arg2 != 0:
        values = ['$\\nu$','$P(\\nu)$',0,arg1[-1],0,np.max(arg0) + 10.0,'LS spectrum (initial)']
    else:
        values = ['$\\nu$','$P(\\nu)$',0,arg1[-1],0,np.max(arg0) + 10.0,'LS spectrum of residuals']
    entries = []
    for i in range(len(names)):
        tk.Label(frame_2,text = names[i],font = "Verdana 13 bold").grid(row = 2 * i,column = 0,
                                                                        padx = int(main_win.winfo_screenwidth() * 0.01))
        entries.append(tk.Entry(frame_2,width = 18))
        entries[-1].insert(0,values[i])
        entries[-1].grid(row = 2 * i,column = 1)
    for i in range(len(names)):
        tk.Label(frame_2,text = "").grid(row = 2 * i + 1,column = 0)
    screen_fig()
    tk.Button(frame_2,text = "Replot",font = "Verdana 13 bold",command = screen_fig).grid(row = 2 * len(names),column = 0)
    tk.Button(frame_2,text = "Save",font = "Verdana 13 bold",command = save_fig).grid(row = 2 * len(names),column = 1)
    tk.Label(frame_2,text = "").grid(row = 2 * len(names) + 1,column = 0)
    tk.Button(frame_2,text = "Reset",font = "Verdana 13 bold",command = reset_fig).grid(row = 2 * len(names) + 2,column = 0)
#############################################################################
#############################################################################
############################## RES_PLOT #####################################
#############################################################################
def res_plot(main_win,arg0,arg1,arg2,arg3):
    
    ticks = np.arange(0,len(arg0),len(arg0) / 7,dtype = int)
    t = np.arange(1,len(arg0) + 1,dtype = int)
    time = np.array(arg3,dtype = str)
    ticks_vec = t[ticks]
    time_label = time[ticks]
    
    pn_norm_notnan = arg2[~np.isnan(arg2)]
    outlier_lim = 3.0
    num_outliers_max = len(pn_norm_notnan[pn_norm_notnan > outlier_lim])
    num_outliers_min = len(pn_norm_notnan[pn_norm_notnan < -outlier_lim])
    num_outliers = num_outliers_max + num_outliers_min
    
    def save_fig():
        path_tot = askdirectory()
        plt.figure(figsize = (12,9))
        plt.rc('text',usetex = True)
        plt.rc('font',family = 'serif')
    
        plt.subplot(2,1,1)
        plt.plot(t,arg0)
        plt.xlim(int(entries[0].get()),int(entries[1].get()))
        plt.ylim(float(entries[5].get()),float(entries[6].get()))
        plt.xticks(ticks,'')
        plt.ylabel(entries[2].get())
        plt.title(entries[4].get())
        plt.margins(0.2)
        plt.subplots_adjust(hspace = 0.0)
    
        plt.subplot(2,1,2)
        sigma = '%.2f' % arg1
        if int(matplotlib.__version__.split('.')[0]) == 2:
            plt.bar(t,arg2,width = 10,label = 'num outl = ' + str(num_outliers))
        else:
            plt.bar(t,arg2,width = 0.1,label = 'num outl = ' + str(num_outliers))
        plt.plot((t[0],t[-1]),(outlier_lim,outlier_lim),'r',label = '$\sigma$ = ' + sigma)
        plt.plot((t[0],t[-1]),(-outlier_lim,-outlier_lim),'r')
        plt.legend(loc = 0)
        plt.xlim(int(entries[0].get()),int(entries[1].get()))
        plt.ylim(float(entries[7].get()),float(entries[8].get()))
        plt.xticks(ticks,time_label)
        plt.ylabel(entries[3].get())
        plt.margins(0.2)
        plt.subplots_adjust(hspace = 0.0)

        plt.savefig(os.path.join(path_tot,'res.pdf'))
        plt.close()

    def screen_fig():
        fig_ts = Figure(figsize = (x_size,y_size))
        a = fig_ts.add_subplot(211)
        a.plot(t,arg0)
        a.set_xlim(int(entries[0].get()),int(entries[1].get()))
        a.set_ylim(float(entries[5].get()),float(entries[6].get()))
        a.set_xticks(ticks_vec)
        a.set_xticklabels('')
        a.set_ylabel(entries[2].get(),fontsize = 15)
        a.set_title(entries[4].get(),fontsize = 15)

        b = fig_ts.add_subplot(212)
        sigma = '%.2f' % arg1
        if int(matplotlib.__version__.split('.')[0]) == 2:
            b.bar(t,arg2,width = 10,label = 'num outl = ' + str(num_outliers))
        else:
            b.bar(t,arg2,width = 0.1,label = 'num outl = ' + str(num_outliers))
        b.plot((t[0],t[-1]),(outlier_lim,outlier_lim),'r',label = '$\sigma$ = ' + sigma)
        b.plot((t[0],t[-1]),(-outlier_lim,-outlier_lim),'r')
        b.legend(loc = 0)
        b.set_xlim(int(entries[0].get()),int(entries[1].get()))
        b.set_ylim(float(entries[7].get()),float(entries[8].get()))
        b.set_xticks(ticks)
        b.set_xticklabels(time_label)
        b.set_ylabel(entries[3].get(),fontsize = 15)
        
        fig_ts.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig_ts,master = frame_1)
        canvas.get_tk_widget().grid(row = 0,column = 0)
        canvas.draw()

    def reset_fig():
        for i in range(len(entries)):
            entries[i].delete(0,tk.END)
            entries[i].insert(0,values[i])
        screen_fig()
    
    top = tk.Toplevel(main_win)
    top.geometry("%dx%d" % (int(main_win.winfo_screenwidth() * 0.93 * 0.85),
                            int(main_win.winfo_screenheight() * 0.65)))
    top.wm_title("Residuals")
    top.resizable(width = False,height = False)
    
    frame_1 = tk.Frame(top)
    frame_1.grid(row = 0,column = 0)
    
    frame_2 = tk.Frame(top)
    frame_2.grid(row = 0,column = 1)
    names = ["X Limit (left)","X Limit (right)","Y Label (top)","Y Label (bottom)","Title",
             "Y1 Limit (bottom)","Y1 Limit (top)","Y2 Limit (bottom)","Y2 Limit (top)"]
    values = [t[0],t[-1],'$N_t$','$N_t^{norm}$','Residuals / Normalised residuals',np.min(arg0[~np.isnan(arg0)]) - 10.0,
              np.max(arg0[~np.isnan(arg0)]) + 10.0,np.min(arg2[~np.isnan(arg0)]) - 1.0,np.max(arg2[~np.isnan(arg0)]) + 1.0]
    entries = []
    for i in range(len(names)):
        tk.Label(frame_2,text = names[i],font = "Verdana 13 bold").grid(row = 2 * i,column = 0,
                                                                        padx = int(main_win.winfo_screenwidth() * 0.01))
        entries.append(tk.Entry(frame_2,width = 18))
        entries[-1].insert(0,values[i])
        entries[-1].grid(row = 2 * i,column = 1)
    for i in range(len(names)):
        tk.Label(frame_2,text = "").grid(row = 2 * i + 1,column = 0)
    screen_fig()
    tk.Button(frame_2,text = "Replot",font = "Verdana 13 bold",command = screen_fig).grid(row = 2 * len(names),column = 0)
    tk.Button(frame_2,text = "Save",font = "Verdana 13 bold",command = save_fig).grid(row = 2 * len(names),column = 1)
    tk.Label(frame_2,text = "").grid(row = 2 * len(names) + 1,column = 0)
    tk.Button(frame_2,text = "Reset",font = "Verdana 13 bold",command = reset_fig).grid(row = 2 * len(names) + 2,column = 0)
#############################################################################
#############################################################################
############################## DFA_PLOT #####################################
#############################################################################
def dfa_plot(main_win,arg0,arg1,arg2,arg3):
    
    def save_fig():
        path_tot = askdirectory()
        plt.rc('text',usetex = True)
        plt.rc('font',family = 'serif')
        plt.plot(np.log(arg0),np.log(arg1),'o',label = '$H$ = ' + arg3)
        plt.plot(np.log(arg0),arg2,'r')
        plt.legend(loc = 0)
        plt.xlim(float(entries[0].get()),float(entries[1].get()))
        plt.ylim(float(entries[2].get()),float(entries[3].get()))
        plt.xlabel(entries[4].get())
        plt.ylabel(entries[5].get())
        plt.title(entries[6].get())
        plt.savefig(os.path.join(path_tot,'dfa.pdf'))
        plt.close()
    
    def screen_fig():
        fig_ts = Figure(figsize = (x_size,y_size))
        a = fig_ts.add_subplot(111)
        a.plot(np.log(arg0),np.log(arg1),'o',label = '$H$ = ' + arg3)
        a.plot(np.log(arg0),arg2,'r')
        a.legend(loc = 0)
        a.set_xlim(float(entries[0].get()),float(entries[1].get()))
        a.set_ylim(float(entries[2].get()),float(entries[3].get()))
        a.set_xlabel(entries[4].get())
        a.set_ylabel(entries[5].get())
        a.set_title(entries[6].get())
        fig_ts.tight_layout()
        canvas = FigureCanvasTkAgg(fig_ts,master = frame_1)
        canvas.get_tk_widget().grid(row = 0,column = 0)
        canvas.draw()

    def reset_fig():
        for i in range(len(entries)):
            entries[i].delete(0,tk.END)
            entries[i].insert(0,values[i])
        screen_fig()
    
    top = tk.Toplevel(main_win)
    top.geometry("%dx%d" % (int(main_win.winfo_screenwidth() * 0.93 * 0.85),
                            int(main_win.winfo_screenheight() * 0.65)))
    top.wm_title("DFA")
    top.resizable(width = False,height = False)
    
    frame_1 = tk.Frame(top)
    frame_1.grid(row = 0,column = 0)
    
    frame_2 = tk.Frame(top)
    frame_2.grid(row = 0,column = 1)
    names = ["X Limit (left)","X Limit (right)","Y Limit (bottom)","Y Limit (top)","X Label","Y Label","Title"]
    values = [np.log(arg0[0]) - 0.3,np.log(arg0[-1]) + 0.3,np.min(np.log(arg1)) - 1.0,np.max(np.log(arg1)) + 1.0,
              'log$(F(n))$','log$(n)$','DFA fit']
    entries = []
    for i in range(len(names)):
        tk.Label(frame_2,text = names[i],font = "Verdana 13 bold").grid(row = 2 * i,column = 0,
                                                                        padx = int(main_win.winfo_screenwidth() * 0.01))
        entries.append(tk.Entry(frame_2,width = 18))
        entries[-1].insert(0,values[i])
        entries[-1].grid(row = 2 * i,column = 1)
    for i in range(len(names)):
        tk.Label(frame_2,text = "").grid(row = 2 * i + 1,column = 0)
    screen_fig()
    tk.Button(frame_2,text = "Replot",font = "Verdana 13 bold",command = screen_fig).grid(row = 2 * len(names),column = 0)
    tk.Button(frame_2,text = "Save",font = "Verdana 13 bold",command = save_fig).grid(row = 2 * len(names),column = 1)
    tk.Label(frame_2,text = "").grid(row = 2 * len(names) + 1,column = 0)
    tk.Button(frame_2,text = "Reset",font = "Verdana 13 bold",command = reset_fig).grid(row = 2 * len(names) + 2,column = 0)
#############################################################################
#############################################################################
############################## MDFA_PLOT ####################################
#############################################################################
def mdfa_plot(main_win,arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7):
    
    def save_fig():
        path_tot = askdirectory()
        plt.figure(figsize = (11,11))
        plt.rc('text',usetex = True)
        plt.rc('font',family = 'serif')

        plt.subplot(2,2,1)
        plt.plot(np.log(arg0),np.log(arg1[0,:]),'b.')
        plt.plot(np.log(arg0),arg2[:,0],'b',label = 'q = -3')
        plt.plot(np.log(arg0),np.log(arg1[50,:]),'r.')
        plt.plot(np.log(arg0),arg2[:,50],'r',label = 'q = 0')
        plt.plot(np.log(arg0),np.log(arg1[-1,:]),'g.')
        plt.plot(np.log(arg0),arg2[:,-1],'g',label = 'q = 3')
        plt.legend(loc = 0)
        plt.xlim(float(entries[0].get()),float(entries[1].get()))
        plt.ylim(float(entries[2].get()),float(entries[3].get()))
        plt.xlabel(entries[4].get())
        plt.ylabel(entries[5].get())
        plt.title(entries[6].get())
        plt.margins(0.2)
        plt.subplots_adjust(bottom = 0.2)

        plt.subplot(2,2,2)
        plt.plot(arg3,arg4,'b',label = 'h(q)')
        plt.plot((arg3[0],arg3[-1]),(arg5,arg5),'k',label = 'H')
        plt.legend(loc = 0)
        plt.xlim(float(entries[7].get()),float(entries[8].get()))
        plt.ylim(float(entries[9].get()),float(entries[10].get()))
        plt.xlabel(entries[11].get())
        plt.ylabel(entries[12].get())
        plt.title(entries[13].get())
        plt.margins(0.2)
        plt.subplots_adjust(bottom = 0.2)

        plt.subplot(2,2,3)
        plt.plot(arg6,arg7,'b')
        plt.xlim(float(entries[14].get()),float(entries[15].get()))
        plt.ylim(float(entries[16].get()),float(entries[17].get()))
        plt.xlabel(entries[18].get())
        plt.ylabel(entries[19].get())
        plt.title(entries[20].get())
        plt.margins(0.2)
        plt.subplots_adjust(bottom = 0.2)

        plt.savefig(os.path.join(path_tot,'mdfa.pdf'))
        plt.close()

    def screen_fig():
        fig_ts = Figure(figsize = (x_size,y_size))
        a = fig_ts.add_subplot(221)
        a.plot(np.log(arg0),np.log(arg1[0,:]),'b.')
        a.plot(np.log(arg0),arg2[:,0],'b',label = 'q = -3')
        a.plot(np.log(arg0),np.log(arg1[50,:]),'r.')
        a.plot(np.log(arg0),arg2[:,50],'r',label = 'q = 0')
        a.plot(np.log(arg0),np.log(arg1[-1,:]),'g.')
        a.plot(np.log(arg0),arg2[:,-1],'g',label = 'q = 3')
        a.legend(loc = 0)
        a.set_xlim(float(entries[0].get()),float(entries[1].get()))
        a.set_ylim(float(entries[2].get()),float(entries[3].get()))
        a.set_xlabel(entries[4].get())
        a.set_ylabel(entries[5].get())
        a.set_title(entries[6].get())

        b = fig_ts.add_subplot(222)
        b.plot(arg3,arg4,'b',label = 'H(q)')
        b.plot((arg3[0],arg3[-1]),(arg5,arg5),'k',label = 'H')
        b.legend(loc = 0)
        b.set_xlim(float(entries[7].get()),float(entries[8].get()))
        b.set_ylim(float(entries[9].get()),float(entries[10].get()))
        b.set_xlabel(entries[11].get())
        b.set_ylabel(entries[12].get())
        b.set_title(entries[13].get())

        c = fig_ts.add_subplot(223)
        c.plot(arg6,arg7,'b')
        c.set_xlim(float(entries[14].get()),float(entries[15].get()))
        c.set_ylim(float(entries[16].get()),float(entries[17].get()))
        c.set_xlabel(entries[18].get())
        c.set_ylabel(entries[19].get())
        c.set_title(entries[20].get())

        fig_ts.tight_layout()

        canvas = FigureCanvasTkAgg(fig_ts,master = frame_1)
        canvas.get_tk_widget().grid(row = 0,column = 0)
        canvas.draw()

    def reset_fig():
        for i in range(len(entries)):
            entries[i].delete(0,tk.END)
            entries[i].insert(0,values[i])
        screen_fig()
    
    top = tk.Toplevel(main_win)
    top.geometry("%dx%d" % (int(main_win.winfo_screenwidth() * 0.93 * 0.85),
                            int(main_win.winfo_screenheight() * 0.75)))
    top.wm_title("MFDFA")
    top.resizable(width = False,height = False)
    
    frame_1 = tk.Frame(top)
    frame_1.grid(row = 0,column = 0)
    
    frame_2 = tk.Frame(top)
    frame_2.grid(row = 0,column = 1)
    names = ["X1 Limit (left)","X1 Limit (right)","Y1 Limit (bottom)","Y1 Limit (top)","X1 Label","Y1 Label","Title1",
             "X2 Limit (left)","X2 Limit (right)","Y2 Limit (bottom)","Y2 Limit (top)","X2 Label","Y2 Label","Title2",
             "X3 Limit (left)","X3 Limit (right)","Y3 Limit (bottom)","Y3 Limit (top)","X3 Label","Y3 Label","Title3"]
    values = [np.log(arg0[0]),np.log(arg0[-1]),np.min([np.min(np.log(arg1[0,:])),np.min(arg2[:,0]),
              np.min(np.log(arg1[50,:])),np.min(arg2[:,50]),np.min(np.log(arg1[-1,:])),np.min(arg2[:,-1])]) - 1.0,
              np.max([np.max(np.log(arg1[0,:])),np.max(arg2[:,0]),np.max(np.log(arg1[50,:])),np.max(arg2[:,50]),
              np.max(np.log(arg1[-1,:])),np.max(arg2[:,-1])]) + 1.0,'log(n)','log(F(n))','MDFA fit',
              arg3[0],arg3[-1],np.min(arg4) - 0.1,np.max(arg4) + 0.1,'q','H(q)','Generalised Hurst exponent',
              np.min(arg6) - 0.2,np.max(arg6) + 0.2,np.min(arg7) - 0.2,1.2,'$\\alpha$','$f(\\alpha)$',
              'Singularity spectrum']
    entries = []
    for i in range(len(names)):
        tk.Label(frame_2,text = names[i],font = "Verdana 13 bold").grid(row = i,column = 0,
                                                                        padx = int(main_win.winfo_screenwidth() * 0.01))
        entries.append(tk.Entry(frame_2,width = 18))
        entries[-1].insert(0,values[i])
        entries[-1].grid(row = i,column = 1)
    screen_fig()
    tk.Label(frame_2,text = "").grid(row = len(names),column = 0)
    tk.Button(frame_2,text = "Replot",font = "Verdana 13 bold",command = screen_fig).grid(row = len(names) + 1,column = 0)
    tk.Button(frame_2,text = "Save",font = "Verdana 13 bold",command = save_fig).grid(row = len(names) + 1,column = 1)
    tk.Button(frame_2,text = "Reset",font = "Verdana 13 bold",command = reset_fig).grid(row = len(names) + 2,column = 0)
#############################################################################
#############################################################################
############################## MFDFA2_PLOT ##################################
#############################################################################
def MFDFA2_plot(main_win,arg0,arg1,arg2,arg3,arg4,arg5):
    
    def save_fig():
        path_tot = askdirectory()
        plt.figure(figsize = (12,9))
        plt.rc('text',usetex = True)
        plt.rc('font',family = 'serif')

        plt.subplot(2,1,1)
        ax = plt.gca()
        if int(matplotlib.__version__.split('.')[0]) == 2:
            ax.set_facecolor('black')
        else:
            ax.set_axis_bgcolor('black')
        plt.plot(arg0,'y')
        plt.plot(0.5 * np.ones((len(arg0),)),'w')
        plt.plot(np.ones((len(arg0),)),'m')
        plt.plot(1.5 * np.ones((len(arg0),)),'r')
        plt.xlim(float(entries[0].get()),float(entries[1].get()))
        plt.ylim(float(entries[2].get()),float(entries[3].get()))
        plt.xlabel(entries[4].get())
        plt.ylabel(entries[5].get())
        plt.title(entries[6].get())
        plt.margins(0.2)
        plt.subplots_adjust(hspace = 0.3)

        plt.subplot(2,1,2)
        plt.plot(arg1,arg2,'b',label = '$\mu$ = ' + arg4)
        plt.plot(arg1,arg3,'r',linewidth = 2.0,label = '$\sigma$ = ' + arg5)
        plt.legend(loc = 0)
        plt.xlim(float(entries[7].get()),float(entries[8].get()))
        plt.ylim(float(entries[9].get()),float(entries[10].get()))
        plt.ylabel(entries[11].get())
        plt.xlabel(entries[12].get())
        plt.title(entries[13].get())
        plt.margins(0.2)
        plt.subplots_adjust(hspace = 0.3)

        plt.savefig(os.path.join(path_tot,'MFDFA2.pdf'))
        plt.close()

    def screen_fig():
        fig_ts = Figure(figsize = (x_size,y_size))
        a = fig_ts.add_subplot(211)
        ax = fig_ts.gca()
        if int(matplotlib.__version__.split('.')[0]) == 2:
            ax.set_facecolor('black')
        else:
            ax.set_axis_bgcolor('black')
        a.plot(arg0,'y')
        a.plot(0.5 * np.ones((len(arg0),)),'w')
        a.plot(np.ones((len(arg0),)),'m')
        a.plot(1.5 * np.ones((len(arg0),)),'r')
        a.set_xlim(float(entries[0].get()),float(entries[1].get()))
        a.set_ylim(float(entries[2].get()),float(entries[3].get()))
        a.set_xlabel(entries[4].get())
        a.set_ylabel(entries[5].get())
        a.set_title(entries[6].get())

        b = fig_ts.add_subplot(212)
        b.plot(arg1,arg2,'b',label = '$\mu$ = ' + arg4)
        b.plot(arg1,arg3,'r',linewidth = 2.0,label = '$\sigma$ = ' + arg5)
        b.legend(loc = 0)
        b.set_xlim(float(entries[7].get()),float(entries[8].get()))
        b.set_ylim(float(entries[9].get()),float(entries[10].get()))
        b.set_ylabel(entries[11].get())
        b.set_xlabel(entries[12].get())
        b.set_title(entries[13].get())
        
        fig_ts.tight_layout()

        canvas = FigureCanvasTkAgg(fig_ts,master = frame_1)
        canvas.get_tk_widget().grid(row = 0,column = 0)
        canvas.draw()

    def reset_fig():
        for i in range(len(entries)):
            entries[i].delete(0,tk.END)
            entries[i].insert(0,values[i])
        screen_fig()
    
    top = tk.Toplevel(main_win)
    top.geometry("%dx%d" % (int(main_win.winfo_screenwidth() * 0.93 * 0.85),
                            int(main_win.winfo_screenheight() * 0.65)))
    top.wm_title("DFA")
    top.resizable(width = False,height = False)
    
    frame_1 = tk.Frame(top)
    frame_1.grid(row = 0,column = 0)
    
    frame_2 = tk.Frame(top)
    frame_2.grid(row = 0,column = 1)
    names = ["X1 Limit (left)","X1 Limit (right)","Y1 Limit (bottom)","Y1 Limit (top)","X1 Label (top)",
             "Y1 Label (top)","Title1 (top)","X2 Limit (left)","X2 Limit (right)","Y2 Limit (bottom)",
             "Y2 Limit (top)","X2 Label (top)","Y2 Label (bottom)","Title2 (bottom)"]
    values = [0,len(arg0),0,3,'time','$H_t$','local Hurst exponent',np.min(arg1) - 0.2,np.max(arg1) + 0.2,
              0,np.max(arg2) * 11 / 10,'P($H_t$)','$H_t$','Prob distr of $H_t$']
    entries = []
    for i in range(len(names)):
        tk.Label(frame_2,text = names[i],font = "Verdana 13 bold").grid(row = i,column = 0,
                                                                        padx = int(main_win.winfo_screenwidth() * 0.01))
        entries.append(tk.Entry(frame_2,width = 18))
        entries[-1].insert(0,values[i])
        entries[-1].grid(row = i,column = 1)
    screen_fig()
    tk.Label(frame_2,text = "").grid(row = len(names),column = 0)
    tk.Button(frame_2,text = "Replot",font = "Verdana 13 bold",command = screen_fig).grid(row = len(names) + 1,column = 0)
    tk.Button(frame_2,text = "Save",font = "Verdana 13 bold",command = save_fig).grid(row = len(names) + 1,column = 1)
    tk.Label(frame_2,text = "").grid(row = len(names) + 2,column = 0)
    tk.Button(frame_2,text = "Reset",font = "Verdana 13 bold",command = reset_fig).grid(row = len(names) + 3,column = 0)
#############################################################################
