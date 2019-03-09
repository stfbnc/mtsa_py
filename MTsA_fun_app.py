import sys
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelmax
if sys.version_info[0] == 2:
    import tkMessageBox as tkM
else:
    import tkinter.messagebox as tkM
from scipy.stats import norm

#############################################################################
############################## LOAD_FILE ####################################
###################### last modified 12/03/2017 #############################
### Reads file and prepares vectors for the analysis
### INPUT ARGUMENTS
### path_tot -> path to the main folder
### file_name -> file name
### data_col -> numer of data column
### time_col -> number of times column
### Delta_t -> sampling time
### in_flag -> 0 for main window plot, 1 for analysis
### OUTPUT ARGUMENTS
### pn -> data vector
### time -> times vector
### nan_pos_pn -> position of missing data
### t -> vector of subsequent integers spaced by sampling time and with the
###      same length of time
### time_label -> sub-vector of time for the x axis of plots
### ticks -> sub-vector of t with positions of time_label elements
#############################################################################
def load_file(file_name,data_col,time_col,Delta_t):

    file = open(file_name,'r')
    pn = []
    time = []
    for line in file:
        cols = line.split()
        pn.append(cols[data_col-1])
        time.append(cols[time_col-1])
    file.close()
    pn = np.array(pn,dtype = float)

    nan_pos_pn = []
    not_nan_pos = []
    for i in range(len(pn)):
        if np.isnan(pn[i]):
            nan_pos_pn.append(i)
        else:
            not_nan_pos.append(i)

    if not_nan_pos[-1] != len(pn):
        pn = pn[0:not_nan_pos[-1] + 1]
        time = time[0:not_nan_pos[-1] + 1]

    if not_nan_pos[0] != 0:
        pn = pn[not_nan_pos[0]:]
        time = time[not_nan_pos[0]:]

    nan_percentage = (float(len(nan_pos_pn) / len(pn))) * 100.0
    if nan_percentage > 20.0:
        tkM.showerror("",
        "Too much missing data for an accurate analysis! Please consider a resample",
        icon='error')

    t_fin = int((len(pn) - 1) * Delta_t + 1)
    t = np.arange(1,t_fin,Delta_t,dtype = int)
    t = np.append(t,t_fin)
    
    return pn,time,nan_pos_pn,t
#############################################################################
#############################################################################
############################ TREND_DETREND ##################################
###################### last modified 28/03/2017 #############################
### Detrends time series and plots data and detrended data all together
### INPUT ARGUMENTS
### time_col -> number of times column
### time -> times vector as in the file
### pn -> data vector
### t -> vector of subsequent integers spaced by sampling time and with the
###      same length of time
### typeoffit -> type of fit to be performed on data (none,poly1-5,exp)
### time_label -> sub-vector of time for the x axis of plots
### ticks -> sub-vector of t with positions of time_label elements
### OUTPUT ARGUMENTS
### pn -> detrended data
### fitted_curve -> trend
#############################################################################
def trend_detrend(pn,t,typeoffit):

    fit_type = {'none':-1,'poly1':1,'poly2':2,'poly3':3,'poly4':4,
                'poly5':5,'exp':0}
    deg = fit_type.get(typeoffit)

    if deg != -1:
        x_fit = []
        y_fit = []
        for i in range(len(pn)):
            if not np.isnan(pn[i]):
                x_fit.append(t[i])
                y_fit.append(pn[i])
        x_fit = np.array(x_fit,dtype = float)
        y_fit = np.array(y_fit,dtype = float)
        if deg != 0:
            fitted_curve = np.polyfit(x_fit,y_fit,deg)
            fitted_curve = np.poly1d(fitted_curve)
            fitted_curve = fitted_curve(t)
        else:
            def exp_func(x,a,b):
                return a * np.exp(b * x)
            popt,pcov = curve_fit(exp_func,x_fit,y_fit,p0 = (1,1e-6))
            fitted_curve = exp_func(t,*popt)
        for i in range(len(pn)):
            if np.isnan(pn[i]):
                fitted_curve[i] = np.nan
        pn = pn - fitted_curve
    else:
        fitted_curve = 0.0

    return pn,fitted_curve
#############################################################################
#############################################################################
################################# GLS #######################################
###################### last modified 03/02/2017 #############################
### Coumputes the generalised lomb spectrum
### INPUT ARGUMENTS
### ts_vector -> data vector
### t_vector -> times vector (subsequent integers spaced by sampling time)
### frequencies -> vector of frequencies for spectrum computation
### OUTPUT ARGUMENTS
### P -> generalised lomb spectrum
#############################################################################
def gls(ts_vector,t_vector,frequencies):

    ts_vector_not_nan = []
    t_vector_not_nan = []
    for i in range(len(ts_vector)):
        if not np.isnan(ts_vector[i]):
            ts_vector_not_nan.append(ts_vector[i])
            t_vector_not_nan.append(t_vector[i])
    ts_vector_not_nan = np.array(ts_vector_not_nan,dtype = float)
    t_vector_not_nan = np.array(t_vector_not_nan,dtype = float)
    N = len(ts_vector_not_nan)
    err_vector_not_nan = np.ones((N,),dtype = float)
    W = np.sum(1.0 / (err_vector_not_nan ** 2.0))
    w_err = 1.0 / (W * err_vector_not_nan ** 2.0)
    ts_vector_not_nan -= np.mean(ts_vector_not_nan)
    sum_dev = np.sum(w_err * (ts_vector_not_nan ** 2.0))
    P = np.zeros(len(frequencies),dtype = float)
    for i in range(len(frequencies)):
        wt = 2.0 * np.pi * frequencies[i] * t_vector_not_nan
        swt = np.sin(wt)
        cwt = np.cos(wt)
        Ss2wt = 2.0 * np.sum(w_err * cwt * swt) - 2.0 * np.sum(w_err * cwt) * np.sum(w_err * swt)
        Sc2wt = np.sum(w_err * (cwt - swt) * (cwt + swt)) - np.sum(w_err * cwt) ** 2.0 + np.sum(w_err * swt) ** 2.0
        wtau = 0.5 * np.arctan2(Ss2wt,Sc2wt)
        swtau = np.sin(wtau)
        cwtau = np.cos(wtau)
        swttau = swt * cwtau - cwt * swtau
        cwttau = cwt * cwtau + swt * swtau
        P[i] = (np.sum(w_err * ts_vector_not_nan * cwttau) ** 2.0) / (np.sum(w_err * cwttau * cwttau) - np.sum(w_err * cwttau) ** 2.0) + (np.sum(w_err * ts_vector_not_nan * swttau) ** 2.0) / (np.sum(w_err * swttau * swttau) - np.sum(w_err * swttau) ** 2.0)
        P[i] = N * P[i] / (2.0 * sum_dev)

    return P
#############################################################################
#############################################################################
############################# LOMB_SCARGLE ##################################
###################### last modified 03/02/2017 #############################
### Computes threshold for the lomb spectrum and plots it
### INPUT ARGUMENTS
### pn -> data vector
### t -> times vector (subsequent integers spaced by sampling time)
### ofac -> oversampling factor
### Delta_t -> sampling time
### OUTPUT ARGUMENTS
### PNT -> generalised lomb spectrum
### freq -> frequencies vector
### pth -> threshold value
#############################################################################
def lomb_scargle(pn,t,ofac,Delta_t):

    if (len(pn) * ofac) % 2 == 1:
        freq = np.linspace(1,np.floor(0.5 * len(pn) * ofac),
                           np.floor(0.5 * len(pn) * ofac)) / (len(pn) * Delta_t * ofac)
    else:
        freq = np.linspace(1,np.floor(0.5 * len(pn) * ofac) - 1,
                           np.floor(0.5 * len(pn) * ofac) - 1) / (len(pn) * Delta_t * ofac)
    PNT = gls(pn,t,freq)
    M = 2.0 * len(freq) / ofac
    peak_prob = 0.95
    pth = (-np.log(1.0 - peak_prob ** (1.0 / M)))

    return PNT,freq,pth
#############################################################################
#############################################################################
########################### PEAKS_FILTER_ALL ################################
###################### last modified 28/03/2017 #############################
### Finds and filters all peaks in the spectrum above a threshold, since
### no more peaks are found
### INPUT ARGUMENTS
### max_iter -> maximum number of iteration for peaks
### PNT -> lomb spectrum
### freq -> frequencies vector
### pth -> threshold value for the lomb spectrum
### units_converter -> factor to express periods in years
### year_in -> initial year of data
### year_fin -> final year of data
### t -> times vector (subsequent integers spaced by sampling time)
### nan_pos_pn -> position of missing data in data vector
### Delta_t -> sampling time
### pn -> data vector
### ofac -> oversampling factor
### OUTPUT ARGUMENTS
### freq_fig -> vector with frequency associated to periodicities
### percent -> percentage of periodic peaks
### num_filter -> number of filter round
#############################################################################
def peaks_filter_all(max_iter,PNT,freq,pth,units_converter,year_in,year_fin,
                     t,nan_pos_pn,Delta_t,pn,ofac):

    tot_spectrum = np.sum(PNT)
    PNT_all = PNT.copy()
    PNT_single = PNT.copy()
    freq_single = freq.copy()
    part_over_tot = 0.0
    iter_peaks = 0
    freq_fig = []
    percent = []
    num_filter = []
    
    while iter_peaks < max_iter:
        iter_peaks += 1
        pks_ind = argrelmax(PNT_single)
        pks_ind = np.array(pks_ind[0],dtype = int)
        pks = []
        for i in pks_ind:
            pks.append(PNT_single[i])
        pks = np.array(pks,dtype = float)
        if PNT_single[1] < PNT_single[0]:
            pks = np.append(pks,PNT_single[0])
        if PNT_single[-2] < PNT_single[-1]:
            pks = np.append(pks,PNT_single[-1])
        num_peaks = len(pks[pks > pth])
        if num_peaks != 0:
            ord_pk = pks[pks > pth]
            interval = []
            for it_peak in range(len(ord_pk)):
                for i in range(len(PNT_single)):
                    if ord_pk[it_peak] == PNT_single[i]:
                        locs_new = freq_single[i]
                        j = i
                if j == 0:
                    interval.append(freq_single[0])
                    x1 = 0
                    for k in np.arange(j + 1,len(freq_single)):
                        if PNT_single[k] > PNT_single[k - 1]:
                            interval.append(freq_single[k - 1])
                            x2 = k - 1
                            break
                        if k == len(PNT_single) - 1:
                            interval.append(freq_single[-1])
                            x2 = len(freq_single) - 1
                elif j == len(freq_single) - 1:
                    for k in np.arange(j - 1,-1,-1):
                        if PNT_single[k] > PNT_single[k + 1]:
                            interval.append(freq_single[k + 1])
                            x1 = k + 1
                            break
                        if k == 0:
                            interval.append(freq_single[0])
                            x1 = 0
                    interval.append(freq_single[-1] + (freq_single[-1] - interval[0]) + freq_single[0])
                    x2 = len(freq_single) - 1
                else:
                    for k in np.arange(j - 1,-1,-1):
                        if PNT_single[k] > PNT_single[k + 1]:
                            interval.append(freq_single[k + 1])
                            x1 = k + 1
                            break
                        if k == 0:
                            interval.append(freq_single[0])
                            x1 = 0
                    for k in np.arange(j + 1,len(freq_single)):
                        if PNT_single[k] > PNT_single[k - 1]:
                            interval.append(freq_single[k - 1])
                            x2 = k - 1
                            break
                        if k == len(PNT_single) - 1:
                            interval.append(freq_single[-1])
                            x2 = len(freq_single) - 1
            
                sum_PNT = np.sum(PNT_all[x1:x2 + 1])
                ratio = (sum_PNT / tot_spectrum) * 100.0
                part_over_tot += ratio
                for i in range(len(PNT_all)):
                    if i >= x1 and i <= x2:
                        PNT_all[i] = 0.0
    
                if (1.0 / locs_new) > (units_converter * (year_fin - year_in + 1)):
                    freq_fig.append(units_converter * (year_fin - year_in + 1))
                    percent.append(ratio)
                    num_filter.append(iter_peaks)
                    if iter_peaks > 1:
                        if freq_fig[iter_peaks - 1] == freq_fig[iter_peaks - 2]:
                            break
                else:
                    freq_fig.append(1.0 / locs_new)
                    percent.append(ratio)
                    num_filter.append(iter_peaks)
                    if iter_peaks > 1:
                        if freq_fig[iter_peaks - 1] == freq_fig[iter_peaks - 2]:
                            break
        
            interval = np.array(interval,dtype = float)
            
            it_filt = 0
            while it_filt <= 10:
                it_filt += 1
                m = np.mean(pn[~np.isnan(pn)])
                dataContent = pn.copy()
                nandata=list(nan_pos_pn)
                if len(nandata) != 0:
                    dataContent[np.isnan(dataContent)] = np.interp(t[np.isnan(dataContent)],
                                t[~np.isnan(dataContent)],dataContent[~np.isnan(dataContent)])
                Ts = Delta_t
                data = dataContent.copy()
                data -= m
                sz = len(data)
                idata = np.fft.fft(data,ofac * sz)
                cont_freq_fft = [0]
                if (len(pn) * ofac) % 2 == 1:
                    fdata = np.concatenate([cont_freq_fft,freq_single,freq_single[::-1]])
                else:
                    freq_plus = [np.floor(0.5 * len(pn) * ofac) / (len(pn) * Delta_t * ofac)]
                    fdata = np.concatenate([cont_freq_fft,freq_single,freq_plus,freq_single[::-1]])
                I = np.ones(len(fdata))
                for i in range(0,it_peak + 1):
                    i_int = np.nonzero(np.logical_and(fdata >= interval[2 * i],
                                                      fdata <= interval[2 * i + 1]))[0]
                    I[i_int] = 0.0
                I = np.concatenate([I,np.zeros(len(idata) - len(I))])
                idata = idata * I
                pn = np.fft.ifft(idata)
                pn = pn.real
                pn += m
                if sys.version_info[0] == 2:
                    pn = pn[0:(len(pn) / ofac)]
                else:
                    pn = pn[0:(len(pn) // ofac)]
                pn[nan_pos_pn] = np.nan
                PNT_single = gls(pn,t,freq_single)
        else:
            break

    return freq_fig,percent,num_filter
#############################################################################
#############################################################################
############################# PEAKS_FILTER ##################################
###################### last modified 28/03/2017 #############################
### Finds and filters peaks in the spectrum above a threshold
### INPUT ARGUMENTS
### max_iter -> maximum number of iteration for peaks
### PNT -> lomb spectrum
### freq -> frequencies vector
### pth -> threshold value for the lomb spectrum
### units_converter -> factor to express periods in years
### year_in -> initial year of data
### year_fin -> final year of data
### t -> times vector (subsequent integers spaced by sampling time)
### nan_pos_pn -> position of missing data in data vector
### Delta_t -> sampling time
### pn -> data vector
### ofac -> oversampling factor
### OUTPUT ARGUMENTS
### pn -> residual data after filtering
#############################################################################
def peaks_filter(max_iter,PNT,freq,pth,units_converter,year_in,year_fin,t,
                 nan_pos_pn,Delta_t,pn,ofac):

    PNT_single = PNT.copy()
    freq_single = freq.copy()
    iter_peaks = 0

    while iter_peaks < max_iter:
        iter_peaks += 1
        pks_ind = argrelmax(PNT_single)
        pks_ind = np.array(pks_ind[0],dtype = int)
        pks = []
        for i in pks_ind:
            pks.append(PNT_single[i])
        pks = np.array(pks,dtype = float)
        if PNT_single[1] < PNT_single[0]:
            pks = np.append(pks,PNT_single[0])
        if PNT_single[-2] < PNT_single[-1]:
            pks = np.append(pks,PNT_single[-1])
        num_peaks = len(pks[pks > pth])
        if num_peaks != 0:
            ord_pk = pks[pks > pth]
            interval = []
            for it_peak in range(len(ord_pk)):
                for i in range(len(PNT_single)):
                    if ord_pk[it_peak] == PNT_single[i]:
                        locs_new = freq_single[i]
                        j = i
                if j == 0:
                    interval.append(freq_single[0])
                    x1 = 0
                    for k in np.arange(j + 1,len(freq_single)):
                        if PNT_single[k] > PNT_single[k - 1]:
                            interval.append(freq_single[k - 1])
                            x2 = k - 1
                            break
                        if k == len(PNT_single) - 1:
                            interval.append(freq_single[-1])
                            x2 = len(freq_single) - 1
                elif j == len(freq_single) - 1:
                    for k in np.arange(j - 1,-1,-1):
                        if PNT_single[k] > PNT_single[k + 1]:
                            interval.append(freq_single[k + 1])
                            x1 = k + 1
                            break
                        if k == 0:
                            interval.append(freq_single[0])
                            x1 = 0
                    interval.append(freq_single[-1] + (freq_single[-1] - interval[0]) + freq_single[0])
                    x2 = len(freq_single) - 1
                else:
                    for k in np.arange(j - 1,-1,-1):
                        if PNT_single[k] > PNT_single[k + 1]:
                            interval.append(freq_single[k + 1])
                            x1 = k + 1
                            break
                        if k == 0:
                            interval.append(freq_single[0])
                            x1 = 0
                    for k in np.arange(j + 1,len(freq_single)):
                        if PNT_single[k] > PNT_single[k - 1]:
                            interval.append(freq_single[k - 1])
                            x2 = k - 1
                            break
                        if k == len(PNT_single) - 1:
                            interval.append(freq_single[-1])
                            x2 = len(freq_single) - 1

            interval = np.array(interval,dtype = float)

            it_filt = 0
            while it_filt <= 10:
                it_filt += 1
                m = np.mean(pn[~np.isnan(pn)])
                dataContent = pn.copy()
                nandata=list(nan_pos_pn)
                if len(nandata) != 0:
                    dataContent[np.isnan(dataContent)] = np.interp(t[np.isnan(dataContent)],
                                t[~np.isnan(dataContent)],dataContent[~np.isnan(dataContent)])
                Ts = Delta_t
                data = dataContent.copy()
                data -= m
                sz = len(data)
                idata = np.fft.fft(data,ofac * sz)
                cont_freq_fft = [0]
                if (len(pn) * ofac) % 2 == 1:
                    fdata = np.concatenate([cont_freq_fft,freq_single,freq_single[::-1]])
                else:
                    freq_plus = [np.floor(0.5 * len(pn) * ofac) / (len(pn) * Delta_t * ofac)]
                    fdata = np.concatenate([cont_freq_fft,freq_single,freq_plus,freq_single[::-1]])
                I = np.ones(len(fdata))
                for i in range(0,it_peak + 1):
                    i_int = np.nonzero(np.logical_and(fdata >= interval[2 * i],
                                                      fdata <= interval[2 * i + 1]))[0]
                    I[i_int] = 0.0
                I = np.concatenate([I,np.zeros(len(idata) - len(I))])
                idata = idata * I
                pn = np.fft.ifft(idata)
                pn = pn.real
                pn += m
                if sys.version_info[0] == 2:
                    pn = pn[0:(len(pn) / ofac)]
                else:
                    pn = pn[0:(len(pn) // ofac)]
                pn[nan_pos_pn] = np.nan
        else:
            break

    return pn
#############################################################################
#############################################################################
############################## RES_LOMB #####################################
###################### last modified 03/02/2017 #############################
### Computes lomb spectrum of residuals and finds outliers
### INPUT ARGUMENTS
### pn -> residuals vector
### Delta_t -> sampling time
### t -> times vector (subsequent integers spaced by sampling time)
### time -> times vector as in the file
### ticks -> sub-vector of t with positions of time_label elements
### time_label -> sub-vector of time for the x axis of plots
### pth -> threshold value
### OUTPUT ARGUMENTS
### s_res -> standard deviation of residuals
### pn_norm -> normalised residuals
### freq_single -> frequencies
### PNT_single -> spectrum of residuals
#############################################################################
def res_lomb(pn,Delta_t,t):

    if len(pn) % 2 == 1:
        freq_single = np.linspace(1,np.floor(0.5 * len(pn)),
                                  np.floor(0.5 * len(pn))) / (len(pn) * Delta_t)
    else:
        freq_single = np.linspace(1,np.floor(0.5 * len(pn)) - 1,
                                  np.floor(0.5 * len(pn)) - 1) / (len(pn) * Delta_t)
    PNT_single = gls(pn,t,freq_single)

    m_res = np.nanmean(pn)
    s_res = np.nanstd(pn)
    pn_norm = (pn - m_res) / s_res

    return s_res,pn_norm,freq_single,PNT_single
#############################################################################
#############################################################################
################################# DFA #######################################
###################### last modified 12/03/2017 #############################
### Computes detrended fluctuations analysis coefficient and plots it
### INPUT ARGUMENTS
### pn -> data vector
### min_win -> minimum scale
### rev_seg -> windows forward and backward
### OUTPUT ARGUMENTS
### s -> scales fo fluctuations
### F -> fluctuations
### log_fit -> fit to fluctuations
### H_mono -> Hurst exponent
#############################################################################
def dfa(pn,min_win,rev_seg):

    nan_pos = []
    for i in range(len(pn)):
        if np.isnan(pn[i]):
            nan_pos.append(i)
    N = len(pn)
    t = np.arange(1,N + 1)
    a_ave = np.nanmean(pn)
    pn -= a_ave
    y = np.zeros((N,),dtype = float)
    for i in range(N):
        y[i] = np.nansum(pn[0:i + 1])
    y[nan_pos] = np.nan
    max_win = 10
    end_dfa = np.floor(N / max_win)
    s = np.arange(min_win,end_dfa + 1,dtype = int)

    F = np.zeros((len(s),),dtype = float)
    for i in range(len(s)):
        N_s = int(N / s[i])
        F_nu1 = np.zeros((N_s,),dtype = float)
        if rev_seg == 1:
            F_nu2 = np.zeros((N_s,),dtype = float)
        for v in range(N_s):
            start_lim = v * s[i]
            end_lim = (v + 1) * s[i]
            t_fit = t[start_lim:end_lim]
            y_fit = y[start_lim:end_lim]
            if len(y_fit[np.isnan(y_fit)]) / len(y_fit) < 0.2:
                n_fit = np.polyfit(t_fit[~np.isnan(y_fit)],
                                   y_fit[~np.isnan(y_fit)],1)
                n_fit = np.poly1d(n_fit)
                n_fit = n_fit(t_fit)
                F_nu1[v] = np.nansum((y_fit - n_fit) ** 2.0) / float(len(y_fit[~np.isnan(y_fit)]))
            else:
                F_nu1[v] = np.nan
        if rev_seg == 1:
            for v in range(N_s):
                start_lim = v * s[i] + (N - N_s * s[i])
                end_lim = (v + 1) * s[i] + (N - N_s * s[i])
                t_fit = t[start_lim:end_lim]
                y_fit = y[start_lim:end_lim]
                if float(len(y_fit[np.isnan(y_fit)])) / float(len(y_fit)) < 0.2:
                    n_fit = np.polyfit(t_fit[~np.isnan(y_fit)],
                                       y_fit[~np.isnan(y_fit)],1)
                    n_fit = np.poly1d(n_fit)
                    n_fit = n_fit(t_fit)
                    F_nu2[v] = np.nansum((y_fit - n_fit) ** 2.0) / float(len(y_fit[~np.isnan(y_fit)]))
                else:
                    F_nu2[v] = np.nan
            F_nu = np.concatenate([F_nu1,F_nu2])
        else:
            F_nu = F_nu1
        F[i] = np.sqrt(np.nansum(F_nu) / float(len(F_nu[~np.isnan(F_nu)])))

    log_fit = np.polyfit(np.log(s),np.log(F),1)
    H_mono = '%.2f' % log_fit[0]
    log_fit = np.poly1d(log_fit)
    log_fit = log_fit(np.log(s))

    return s,F,log_fit,H_mono
#############################################################################
#############################################################################
################################ MDFA #######################################
###################### last modified 16/03/2017 #############################
### Multifractal detrended fluctuations analysis
### INPUT ARGUMENTS
### H_mono -> Hurst exponent for dfa
### pn -> data vector
### min_win -> minimum scale
### q_max -> absolute value of maximum order q
### rev_seg -> windows forward and backward
### OUTPUT ARGUMENTS
### s -> scales fo fluctuations
### F -> fluctuations
### MDFA_fit -> fit to fluctuations
### q -> vector of q orders
### H -> Hurst exponent for every q
### alpha -> singularity indexes
### sing_spec -> spectrum of singularities
#############################################################################
def mdfa(H_mono,pn,min_win,q_max,rev_seg):
    
    nan_pos = []
    for i in range(len(pn)):
        if np.isnan(pn[i]):
            nan_pos.append(i)
    N = len(pn)
    t = np.arange(1,N + 1)
    a_ave = np.nanmean(pn)
    pn -= a_ave
    y = np.zeros((N,),dtype = float)
    for i in range(N):
        y[i] = np.nansum(pn[0:i + 1])
    y[nan_pos] = np.nan
    max_win = 10
    end_dfa = np.floor(N / max_win)
    s = np.arange(min_win,end_dfa + 1,dtype = int)
    q = np.linspace(-q_max,q_max,101)
    
    F = np.zeros((len(q),len(s)),dtype = float)
    for i in range(len(s)):
        N_s = int(N / s[i])
        F_nu1 = np.zeros((N_s,),dtype = float)
        if rev_seg == 1:
            F_nu2 = np.zeros((N_s,),dtype = float)
        for v in range(N_s):
            start_lim = v * s[i]
            end_lim = (v + 1) * s[i]
            t_fit = t[start_lim:end_lim]
            y_fit = y[start_lim:end_lim]
            if float(len(y_fit[np.isnan(y_fit)])) / float(len(y_fit)) < 0.2:
                n_fit = np.polyfit(t_fit[~np.isnan(y_fit)],y_fit[~np.isnan(y_fit)],1)
                n_fit = np.poly1d(n_fit)
                n_fit = n_fit(t_fit)
                F_nu1[v] = np.nansum((y_fit - n_fit) ** 2.0) / float(len(y_fit[~np.isnan(y_fit)]))
            else:
                F_nu1[v] = np.nan
        if rev_seg == 1:
            for v in range(N_s):
                start_lim = v * s[i] + (N - N_s * s[i])
                end_lim = (v + 1) * s[i] + (N - N_s * s[i])
                t_fit = t[start_lim:end_lim]
                y_fit = y[start_lim:end_lim]
                if float(len(y_fit[np.isnan(y_fit)])) / float(len(y_fit)) < 0.2:
                    n_fit = np.polyfit(t_fit[~np.isnan(y_fit)],y_fit[~np.isnan(y_fit)],1)
                    n_fit = np.poly1d(n_fit)
                    n_fit = n_fit(t_fit)
                    F_nu2[v] = np.nansum((y_fit - n_fit) ** 2.0) / float(len(y_fit[~np.isnan(y_fit)]))
                else:
                    F_nu2[v] = np.nan
            F_nu = np.concatenate([F_nu1,F_nu2])
        else:
            F_nu = F_nu1
        for k in range(len(q)):
            if q[k] == 0.0:
                F[k,i]=np.exp(np.nansum(np.log(F_nu)) / (2.0 * float(len(F_nu[~np.isnan(F_nu)]))))
            else:
                F[k,i]=(np.nansum(F_nu ** (q[k] / 2.0)) / float(len(F_nu[~np.isnan(F_nu)]))) ** (1.0/q[k])

    H = np.zeros((len(q),),dtype = float)
    MDFA_fit = np.zeros((len(s),len(q)),dtype = float)
    for i in range(len(q)):
        log_fit = np.polyfit(np.log(s),np.log(F[i,:]),1)
        H[i] = log_fit[0]
        log_fit = np.poly1d(log_fit)
        MDFA_fit[:,i] = log_fit(np.log(s))

    tau = H * q - 1
    alpha = np.diff(tau) / (q[1] - q[0])
    sing_spec = q[0 : -1] * alpha - tau[0 : -1]

    return s,F,MDFA_fit,q,H,alpha,sing_spec
#############################################################################
#############################################################################
################################ MFDFA2 #####################################
###################### last modified 23/03/2017 #############################
### Computes Local Hurst exponent (python version of Ihlen's Matlab MFDFA2.m
### INPUT ARGUMENTS
### signal -> data vector
### scale -> scales at which compute H(t)
### m -> order for fluctuations polynomial fit
### OUTPUT ARGUMENTS
### Ht_plot -> H(t) for minimum scale
### Htbin -> bins for histogram
### Ph -> values of histogram
### fit_gauss -> gaussian fit to the histogram
### mu -> mu of fit
### sigma -> sigma of fit
#############################################################################
def MFDFA2(signal,scale,m):

    nan_pos = []
    for i in range(len(signal)):
        if np.isnan(signal[i]):
            nan_pos.append(i)
    X = np.zeros((len(signal),),dtype = float)
    for i in range(len(signal)):
        X[i] = np.nansum(signal[0:i + 1])
    X[nan_pos] = np.nan

    scmin = 10
    scmax = len(signal) / 10
    scale0 = np.arange(scmin,scmax + 1,dtype = int)

    Fq0 = np.zeros((len(scale0),),dtype = float)
    for ns in range(len(scale0)):
        if sys.version_info[0] == 2:
            segments = len(X) / scale0[ns]
        else:
            segments = len(X) // scale0[ns]
        RMS0 = np.zeros((segments,),dtype = float)
        for v in range(segments):
            Index0 = np.arange((v * scale0[ns]),((v + 1) * scale0[ns]),dtype = int)
            X_fit = X[Index0]
            if float(len(X_fit[np.isnan(X_fit)])) / float(len(X_fit)) < 0.2:
                C0 = np.polyfit(Index0[~np.isnan(X_fit)],X_fit[~np.isnan(X_fit)],m)
                C0 = np.poly1d(C0)
                fit0 = C0(Index0)
                RMS0[v] = np.sqrt(np.nanmean((X_fit - fit0) ** 2.0))
            else:
                RMS0[v] = np.nan
        Fq0[ns] = np.exp(0.5 * np.nanmean(np.log(RMS0 ** 2)))
    C = np.polyfit(np.log(scale0),np.log(Fq0),1)
    Hq0 = C[0]
    C = np.poly1d(C)
    Regfit = C(np.log(scale))

    halfmax = int(np.max(scale) / 2.0)
    Time_index = np.arange(halfmax,len(X) - halfmax,dtype = int)
    maxL = len(Time_index)
    RMS = np.zeros((len(scale),len(Time_index)),dtype = float)
    for ns in range(len(scale)):
        halfseg = int(scale[ns] / 2.0)
        for v in Time_index:
            Index = np.arange(v - halfseg,v + halfseg + 1,dtype = int)
            X_fit = X[Index]
            if float(len(X_fit[np.isnan(X_fit)])) / float(len(X_fit)) < 0.2:
                C = np.polyfit(Index[~np.isnan(X_fit)],X_fit[~np.isnan(X_fit)],m)
                C = np.poly1d(C)
                fitt = C(Index)
                RMS[ns,v - Time_index[0]] = np.sqrt(np.nanmean((X_fit - fitt) ** 2.0))
            else:
                RMS[ns,v - Time_index[0]] = np.nan

    Ht = np.zeros((len(scale),len(RMS[0,:])),dtype = float)
    Ht_row = np.zeros((len(scale) * len(RMS[0,:]),),dtype = float)
    for ns in range(len(scale)):
        RMSt = RMS[ns,:]
        resRMS = Regfit[ns] - np.log(RMSt)
        logscale = np.log(maxL) - np.log(scale[ns])
        Ht[ns,:] = resRMS / float(logscale) + Hq0
        Ht_row[(ns * len(resRMS)):((ns + 1) * len(resRMS))] = Ht[ns,:]

    BinNumb = int(np.sqrt(len(Ht_row[~np.isnan(Ht_row)])))
    freq,Htbin = np.histogram(Ht_row[~np.isnan(Ht_row)],bins = BinNumb)
    Htbin = (Htbin[:-1] + Htbin[1:]) / 2.0
    Ph = freq / (float(freq.sum()) * float(Htbin[1] - Htbin[0]))

    param = norm.fit(Ht_row[~np.isnan(Ht_row)])
    fit_gauss = norm.pdf(Htbin,loc = param[0],scale = param[1])
    mu = '%.2f' % param[0]
    sigma = '%.2f' % param[1]
    Ht_plot = Ht[0,:]

    return Ht_plot,Htbin,Ph,fit_gauss,mu,sigma
#############################################################################
