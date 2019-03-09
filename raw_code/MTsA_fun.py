import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelmax
import os
import sys
import shutil
from scipy.stats import norm,genpareto,invgauss,lognorm,t,rayleigh,weibull_min,weibull_max,gamma,burr,expon,kstest

#############################################################################
############################# ARG_PARSER ####################################
###################### last modified 01/04/2017 #############################
### Parses and converts inputs from command line
### OUTPUT ARGUMENTS
### file_name -> path to file and name of the file to analyse
### time_col -> number of times column
### data_col -> numer of data column
### Delta_t -> sampling time
### time_units -> time units
### ofac -> oversampling factor
### typeoffit -> function for fit
### year_in -> initial year
### year_fin -> final year
### rev_seg -> 0 (1) for forward (backward and forward) dfa
### scale_min -> minimum scale for dfa
### scale_MFDFA -> scales for H(t)
#############################################################################
def args_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument("file_name",help = "path to file + file name",type = str)
    parser.add_argument("time_col",help = "number of time column",type = int)
    parser.add_argument("data_col",help = "number of data column",type = int)
    parser.add_argument("Delta_t",help = "sampling time",type = float)
    parser.add_argument("time_units",help = "units of sampling time (seconds, minutes, hours, days, weeks, months)",type = str)
    parser.add_argument("ofac",help = "oversampling factor",type = int)
    parser.add_argument("typeoffit",help = "choice of fit for trend removal (none, poly1, poly2, poly3, poly4, poly5, exp) ",type = str)
    parser.add_argument("year_in",help = "initial year of time series",type = int)
    parser.add_argument("year_fin",help = "final year of time series",type = int)
    parser.add_argument("rev_seg",help = "0 or 1 for DFA and MDFA windows backward and forward",type = int)
    parser.add_argument("scale_min",help = "smaller scale for DFA and MDFA",type = int)
    parser.add_argument("scale_MFDFA",help = "scales for local hurst exponent",type = lambda x: [int(i) for i in x.split(',')])
    args = parser.parse_args()
    file_name = args.file_name
    time_col = args.time_col
    data_col = args.data_col
    Delta_t = args.Delta_t
    time_units = args.time_units
    ofac = args.ofac
    typeoffit = args.typeoffit
    year_in = args.year_in
    year_fin = args.year_fin
    rev_seg = args.rev_seg
    scale_min = args.scale_min
    scale_MFDFA = args.scale_MFDFA

    return file_name,time_col,data_col,Delta_t,time_units,ofac,typeoffit,year_in,year_fin,rev_seg,scale_min,scale_MFDFA
#############################################################################
#############################################################################
############################ INPUT_CHECKS ###################################
###################### last modified 01/04/2017 #############################
### Checks if inputs are valid
### INPUT ARGUMENTS
### time_col -> number of times column
### data_col -> numer of data column
### Delta_t -> sampling time
### year_in -> initial year
### year_fin -> final year
### ofac -> oversampling factor
### time_units -> time units
### typeoffit -> function for fit
### rev_seg -> 0 (1) for forward (backward and forward) dfa
### scale_min -> minimum scale for dfa
### scale_MFDFA -> scales for H(t)
### OUTPUT ARGUMENTS
### units_converter -> converter for time units
#############################################################################
def input_checks(time_col,data_col,Delta_t,year_in,year_fin,ofac,time_units,typeoffit,rev_seg,scale_min,scale_MFDFA):

    if time_col <= 0:
        print "Time column must be a positive number greater than zero!"
        sys.exit()
    if data_col <= 0:
        print "Data column must be a positive number greater than zero!"
        sys.exit()
    if Delta_t <= 0:
        print "Sampling time must be a positive number greater than zero!"
        sys.exit()
    if year_in <= 0 or year_fin <= 0 or year_fin <= year_in:
        print "Invalid years!"
        sys.exit()
    if ofac <= 0:
        print "Oversampling factor must be a positive number greater than zero!"
        sys.exit()
    units_choices = {'seconds':31536000,'minutes':525600,'hours':8760,'days':365,'weeks':52,'months':12}
    units_converter = units_choices.get(time_units)
    if units_converter not in units_choices.values():
        print "Invalid time units!"
        sys.exit()
    fit_type = ['none','poly1','poly2','poly3','poly4','poly5','exp']
    if typeoffit not in fit_type:
        print "Invalid fit function!"
        sys.exit()
    if rev_seg not in [0,1]:
        print "!"
        sys.exit()
    if scale_min < 3:
        print "Minimum scale for dfa should be at least 3!"
        sys.exit()
    for i in scale_MFDFA:
        if i < 3:
            print "Scales for local Hurst exponent can not be less than 3!"
            sys.exit()

    return units_converter
#############################################################################
#############################################################################
############################ MAIN_FOLDER ####################################
###################### last modified 01/04/2017 #############################
### Creates the main folder where to save analisys' outputs
### INPUT ARGUMENTS
### file_name -> path to file and name of the file to analyse
### data_col -> numer of data column
### OUTPUT ARGUMENTS
### path_tot -> path to the main folder
#############################################################################
def main_folder(file_name,data_col):

    try:
        slash = file_name[::-1].index('/')
        slash = len(file_name) - 1 - slash
        path_fm = file_name[0:slash]
        file_name_slash = file_name[slash + 1:]
    except:
        path_fm = '.'

    try:
        last_dot = file_name_slash[::-1].index('.')
    except:
        print "No file extension found!"
        sys.exit()
    last_dot = len(file_name_slash) - 1 - last_dot
    file_name_noext = file_name_slash[0:last_dot]
    path_tot = os.path.join(path_fm,file_name_noext + '_col' + str(data_col))
    if not os.path.exists(path_tot):
        os.makedirs(path_tot)
    else:
        overwrite = raw_input("You already analyzed this time series or another time series in the same file. Continue? [Y/n]  ")
        if overwrite == 'Y':
            shutil.rmtree(path_tot)
            os.makedirs(path_tot)
        else:
            sys.exit()

    return path_tot
#############################################################################
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
### OUTPUT ARGUMENTS
### pn -> data vector
### time -> times vector
### nan_pos_pn -> position of missing data
### t -> vector of subsequent integers spaced by sampling time and with the
###      same length of time
### time_label -> sub-vector of time for the x axis of plots
### ticks -> sub-vector of t with positions of time_label elements
#############################################################################
def load_file(path_tot,file_name,data_col,time_col,Delta_t):

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
        print('Too much missing data for an accurate analysis!\n')
        sys.exit()

    t_fin = (len(pn) - 1) * Delta_t + 1
    t = np.arange(1,t_fin,Delta_t,dtype = float)
    t = np.append(t,t_fin)
    time_label = []
    ticks = []
    for i in range(len(time)):
        if i % np.ceil(len(time) * 7.0 / 100.0) == 0:
            time_label = np.append(time_label,time[i])
            ticks = np.append(ticks,int(t[i]))
    time_label = np.array(time_label,dtype = str)
    ticks = np.array(ticks,dtype = int)
    t[nan_pos_pn] = np.nan

    return pn,time,nan_pos_pn,t,time_label,ticks
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
### typeoffit -> type of fit to be performed on data (linear,2nd order
###              polynomial,5th order polynomial,exponential)
### path_tot -> path to the main folder
### time_label -> sub-vector of time for the x axis of plots
### ticks -> sub-vector of t with positions of time_label elements
### OUTPUT ARGUMENTS
### pn -> detrended data
#############################################################################
def trend_detrend(time_col,time,pn,t,typeoffit,path_tot,time_label,ticks):

    fit_type = {'none':-1,'poly1':1,'poly2':2,'poly3':3,'poly4':4,'poly5':5,'exp':0}
    deg = fit_type.get(typeoffit)

    if deg != -1:
        x_fit = []
        y_fit = []
        for i in range(len(t)):
            if not np.isnan(t[i]):
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

    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')
    plt.plot(t,pn,'r',label = 'original')
    if deg != -1:
        plt.plot(t,fitted_curve,'k',label = 'trend')
        plt.plot(t,pn-fitted_curve,'b',label = 'detrended')
    plt.legend(loc = 0)
    plt.ylabel('$X_t$')
    plt.xlabel('t')
    plt.title('time series')
    plt.xlim((t[0],t[-1]))
    plt.xticks(ticks,time_label,rotation = 'vertical')
    plt.margins(0.2)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig(os.path.join(path_tot,'ts.pdf'))
    plt.close()

    path_file = os.path.join(path_tot,'detrend.txt')
    f = open(path_file,'w')
    for i in range(len(pn)):
        if deg != -1:
            f.write('%s %f %f %f\n' % (time[i],t[i],pn[i],fitted_curve[i]))
        else:
            f.write('%s %f %f\n' % (time[i],t[i],pn[i]))
    f.close()

    pn = pn-fitted_curve

    return pn
#############################################################################
#############################################################################
############################## NORM_LOMB ####################################
###################### last modified 13/09/2016 #############################
### Coumputes the normalised lomb spectrum
### INPUT ARGUMENTS
### ts_vector -> data vector
### t_vector -> times vector (subsequent integers spaced by sampling time)
### frequencies -> vector of frequencies for spectrum computation
### OUTPUT ARGUMENTS
### P -> normalised lomb spectrum
#############################################################################
def norm_lomb(ts_vector,t_vector,frequencies):

	ts_vector_not_nan = []
	t_vector_not_nan = []
	for i in range(len(ts_vector)):
		if not np.isnan(ts_vector[i]):
			ts_vector_not_nan.append(ts_vector[i])
			t_vector_not_nan.append(t_vector[i])
	ts_vector_not_nan = np.array(ts_vector_not_nan,dtype = float)
	t_vector_not_nan = np.array(t_vector_not_nan,dtype = float)
	sigma2 = np.var(ts_vector_not_nan,dtype = np.float)
	ts_vector_not_nan -= np.mean(ts_vector_not_nan)
	P = np.zeros(len(frequencies),dtype = float)
	for i in range(len(frequencies)):
		wt = 2.0 * np.pi * frequencies[i] * t_vector_not_nan
		swt = np.sin(wt)
		cwt = np.cos(wt)
		Ss2wt = 2.0 * np.sum(cwt * swt)
		Sc2wt = np.sum((cwt - swt) * (cwt + swt))
		wtau = 0.5 * np.arctan2(Ss2wt,Sc2wt)
		swtau = np.sin(wtau)
		cwtau = np.cos(wtau)
		swttau = swt * cwtau - cwt * swtau
		cwttau = cwt * cwtau + swt * swtau
		P[i] = (np.sum(ts_vector_not_nan * cwttau) ** 2.0) / np.sum(cwttau * cwttau) + (np.sum(ts_vector_not_nan * swttau) ** 2.0) / np.sum(swttau * swttau)
		P[i] /= (2.0 * sigma2)

	return P
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
### path_tot -> path to the main folder
### OUTPUT ARGUMENTS
### PNT -> generalised lomb spectrum
### freq -> frequencies vector
### pth -> threshold value
#############################################################################
def lomb_scargle(pn,t,ofac,Delta_t,path_tot):

    if (len(pn) * ofac) % 2 == 1:
        freq = np.linspace(1,np.floor(0.5 * len(pn) * ofac),np.floor(0.5 * len(pn) * ofac)) / (len(pn) * Delta_t * ofac)
    else:
        freq = np.linspace(1,np.floor(0.5 * len(pn) * ofac) - 1,np.floor(0.5 * len(pn) * ofac) - 1) / (len(pn) * Delta_t * ofac)
    PNT = gls(pn,t,freq)
    M = 2.0 * len(freq) / ofac
    peak_prob = 0.95
    pth = (-np.log(1.0 - peak_prob ** (1.0 / M)))

    path_file = os.path.join(path_tot,'spectrum_in.txt')
    f = open(path_file,'w')
    for i in range(len(PNT)):
        f.write('%f %f %f\n' % (freq[i],PNT[i],pth))
    f.close()

    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')
    plt.plot(freq,PNT,'b')
    plt.plot((freq[0],freq[-1]),(pth,pth),'r')
    plt.ylabel('$P(\\nu)$')
    plt.xlabel('$\\nu$')
    plt.xlim(0,freq[-1])
    plt.ylim(0,np.max(PNT)+10.0)
    plt.title('LS spectrum (initial)')
    plt.margins(0.2)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig(os.path.join(path_tot,'spectrum_in.pdf'))
    plt.close()

    return PNT,freq,pth
#############################################################################
#############################################################################
########################### PEAKS_FILTER_ALL ################################
###################### last modified 28/03/2017 #############################
### Finds and filters all peaks in the spectrum above a threshold, since
### no more peaks are found
### INPUT ARGUMENTS
### max_iter -> maximum number of iteration for peaks
### path_tot -> path to the main folder
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
#############################################################################
def peaks_filter_all(max_iter,path_tot,PNT,freq,pth,units_converter,year_in,year_fin,t,nan_pos_pn,Delta_t,pn,ofac):
    
    path_file = os.path.join(path_tot,'freq_peak_all.txt')
    f_peaks = open(path_file,'w')
    
    tot_spectrum = np.sum(PNT)
    PNT_all = PNT.copy()
    PNT_single = PNT.copy()
    freq_single = freq.copy()
    part_over_tot = 0.0
    iter_peaks = 0
    freq_fig = []
    
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
                    if iter_peaks > 1:
                        if freq_fig[iter_peaks - 1] == freq_fig[iter_peaks - 2]:
                            break
                    f_peaks.write('%2d\t%9.8f\t%9.4f\t%5.2f\n' % (iter_peaks,(1.0 / (units_converter * (year_fin - year_in + 1))),units_converter * (year_fin - year_in + 1),ratio))
                else:
                    freq_fig.append(1.0 / locs_new)
                    if iter_peaks > 1:
                        if freq_fig[iter_peaks - 1] == freq_fig[iter_peaks - 2]:
                            break
                    f_peaks.write('%2d\t%9.8f\t%9.4f\t%5.2f\n' % (iter_peaks,locs_new,(1.0 / locs_new),ratio))
        
            interval = np.array(interval,dtype = float)
            
            it_filt = 0
            while it_filt <= 10:
                it_filt += 1
                for i in nan_pos_pn:
                    t[i] = i * Delta_t + 1
                m = np.mean(pn[~np.isnan(pn)])
                dataContent = pn.copy()
                nandata=list(nan_pos_pn)
                if len(nandata) != 0:
                    dataContent[np.isnan(dataContent)] = np.interp(t[np.isnan(dataContent)],t[~np.isnan(dataContent)],dataContent[~np.isnan(dataContent)])
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
                    i_int = np.nonzero(np.logical_and(fdata >= interval[2 * i],fdata <= interval[2 * i + 1]))[0]
                    I[i_int] = 0.0
                I = np.concatenate([I,np.zeros(len(idata) - len(I))])
                idata = idata * I
                pn = np.fft.ifft(idata)
                pn = pn.real
                pn += m
                pn = pn[0:(len(pn) / ofac)]
                pn[nan_pos_pn] = np.nan
                t[nan_pos_pn] = np.nan
                PNT_single = gls(pn,t,freq_single)
        else:
            break

    f_peaks.close()
#############################################################################
#############################################################################
############################# PEAKS_FILTER ##################################
###################### last modified 28/03/2017 #############################
### Finds and filters peaks in the spectrum above a threshold
### INPUT ARGUMENTS
### max_iter -> maximum number of iteration for peaks
### path_tot -> path to the main folder
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
### freq_fig -> vector with time periods
### pn -> residual data after filtering
### sig_phi -> matrix with single filtered components
### sig_to_noise -> total percentage of filtered signal
#############################################################################
def peaks_filter(max_iter,path_tot,PNT,freq,pth,units_converter,year_in,year_fin,t,nan_pos_pn,Delta_t,pn,ofac):

    path_file = os.path.join(path_tot,'freq_peak_percentage.txt')
    f_peaks = open(path_file,'w')

    print('\nfreq        \ttime    \tpercentage\n\n')

    tot_spectrum = np.sum(PNT)
    PNT_single = PNT.copy()
    freq_single = freq.copy()
    pn_phi = pn.copy()
    t_phi = t.copy()
    for i in nan_pos_pn:
        t_phi[i] = i * Delta_t + 1
    part_over_tot = 0.0
    iter_peaks = 0
    freq_fig = []
    sig_phi = np.zeros((1,len(pn_phi)),dtype = float)

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

                sum_PNT = np.sum(PNT[x1:x2 + 1])
                ratio = (sum_PNT / tot_spectrum) * 100.0
                part_over_tot += ratio
                for i in range(len(PNT)):
                    if i >= x1 and i <= x2:
                        PNT[i] = 0.0

                if (1.0 / locs_new) > (units_converter * (year_fin - year_in + 1)):
                    freq_fig.append(units_converter * (year_fin - year_in + 1))
                    if iter_peaks > 1:
                        if freq_fig[iter_peaks - 1] == freq_fig[iter_peaks - 2]:
                            break
                    f_peaks.write('%2d\t%9.8f\t%9.4f\t%5.2f\n' % (iter_peaks,(1.0 / (units_converter * (year_fin - year_in + 1))),units_converter * (year_fin - year_in + 1),ratio))
                    print('%9.8f\t%9.4f\t%5.2f\n' % ((1.0 / (units_converter * (year_fin - year_in + 1))),units_converter * (year_fin - year_in + 1),ratio))
                else:
                    freq_fig.append(1.0 / locs_new)
                    if iter_peaks > 1:
                        if freq_fig[iter_peaks - 1] == freq_fig[iter_peaks - 2]:
                            break
                    f_peaks.write('%2d\t%9.8f\t%9.4f\t%5.2f\n' % (iter_peaks,locs_new,(1.0 / locs_new),ratio))
                    print('%9.8f\t%9.4f\t%5.2f\n' % (locs_new,(1.0 / locs_new),ratio))

            interval = np.array(interval,dtype = float)

            it_filt = 0
            while it_filt <= 10:
                it_filt += 1
                for i in nan_pos_pn:
                    t[i] = i * Delta_t + 1
                m = np.mean(pn[~np.isnan(pn)])
                dataContent = pn.copy()
                nandata=list(nan_pos_pn)
                if len(nandata) != 0:
                    dataContent[np.isnan(dataContent)] = np.interp(t[np.isnan(dataContent)],t[~np.isnan(dataContent)],dataContent[~np.isnan(dataContent)])
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
                    i_int = np.nonzero(np.logical_and(fdata >= interval[2 * i],fdata <= interval[2 * i + 1]))[0]
                    I[i_int] = 0.0
                I = np.concatenate([I,np.zeros(len(idata) - len(I))])
                idata = idata * I
                pn = np.fft.ifft(idata)
                pn = pn.real
                pn += m
                pn = pn[0:(len(pn) / ofac)]
                pn[nan_pos_pn] = np.nan
                t[nan_pos_pn] = np.nan

            it_filt_phi = 0
            freq_phi = freq_single.copy()
            while it_filt_phi <= 10:
                it_filt_phi += 1
                m_phi = np.mean(pn_phi[~np.isnan(pn_phi)])
                dataContent_phi = pn_phi.copy()
                nandata_phi = list(nan_pos_pn)
                if len(nandata_phi) != 0:
                    dataContent_phi[np.isnan(dataContent_phi)] = np.interp(t_phi[np.isnan(dataContent_phi)],t_phi[~np.isnan(dataContent_phi)],dataContent_phi[~np.isnan(dataContent_phi)])
                Ts_phi = Delta_t
                data_phi = dataContent_phi.copy()
                data_phi -= m_phi
                sz_phi = len(data_phi)
                idata_phi = np.fft.fft(data_phi,ofac * sz_phi)
                cont_freq_fft_phi = [0]
                if (len(pn_phi) * ofac) % 2 == 1:
                    fdata_phi = np.concatenate([cont_freq_fft_phi,freq_phi,freq_phi[::-1]])
                else:
                    freq_plus_phi = [np.floor(0.5 * len(pn_phi) * ofac) / (len(pn_phi) * Delta_t * ofac)]
                    fdata_phi = np.concatenate([cont_freq_fft_phi,freq_phi,freq_plus_phi,freq_phi[::-1]])
                I_phi = np.ones(len(fdata_phi))
                for i in range(0,it_peak,2):
                    i_int_phi = np.nonzero(np.logical_and(fdata_phi >= interval[i],fdata_phi <= interval[i + 1]))[0]
                    I_phi[i_int_phi] = 0.0
                I_phi = np.concatenate([I_phi,np.ones(len(idata_phi) - len(I_phi))])
                idata_phi = idata_phi * I_phi
                pn_phi_filt = np.fft.ifft(idata_phi)
                pn_phi_filt = pn_phi_filt.real
                pn_phi_filt += m_phi
                pn_phi_filt = pn_phi_filt[0:(len(pn_phi_filt) / ofac)]
                pn_phi_filt[nan_pos_pn] = np.nan

            sig_phi=np.vstack((sig_phi,pn_phi-pn_phi_filt))
        else:
            break

    sig_phi=sig_phi[1:,:]
    f_peaks.close()

    if part_over_tot != 0.0:
        sig_to_noise = part_over_tot
    else:
        freq_fig=0
        sig_phi=0
        sig_to_noise=0.0

    return freq_fig,pn,sig_phi,sig_to_noise
#############################################################################
#############################################################################
############################## PHASE_FIG ####################################
###################### last modified 28/03/2017 #############################
### Plots for filtered periods
### INPUT ARGUMENTS
### sig_phi -> matrix with single filtered components
### time_col -> number of times column
### time -> times vector as in the file
### time_label -> sub-vector of time for the x axis of plots
### ticks -> sub-vector of t with positions of time_label elements
### freq_fig -> vector with time periods
### path_tot -> path to the main folder
### t -> times vector (subsequent integers spaced by sampling time)
#############################################################################
def phase_fig(sig_phi,time_col,time,time_label,ticks,freq_fig,path_tot,t):

    if not np.isscalar(sig_phi):
        path_phi = os.path.join(path_tot,'phi_signal')
        os.makedirs(path_phi)
        sz_phi = np.shape(sig_phi)[0]
        for i in range(sz_phi):
            plt.rc('text',usetex = True)
            plt.rc('font',family = 'serif')
            plt.plot(t,sig_phi[i,:],'b')
            plt.ylabel('$S_t$')
            plt.xlabel('t')
            plt.title('signal ($T$ = ' + str(freq_fig[i]) + ')')
            plt.xlim((t[0],t[-1]))
            plt.xticks(ticks,time_label,rotation = 'vertical')
            plt.margins(0.2)
            plt.subplots_adjust(bottom = 0.2)
            plt.savefig(os.path.join(path_phi,'phi_signal_' + str(i + 1) + '.pdf'))
            plt.close()
            path_file = os.path.join(path_phi,'sig_' + str(i + 1) + '.txt')
            f = open(path_file,'w')
            for j in range(len(sig_phi[i,:])):
                f.write('%f %f %s\n' % (t[j],sig_phi[i,j],time[j]))
            f.close()
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
### path_tot -> path to the main folder
### OUTPUT ARGUMENTS
### pos_out -> position of outliers
#############################################################################
def res_lomb(pn,Delta_t,t,time,ticks,time_label,pth,path_tot):

    if len(pn) % 2 == 1:
        freq_single = np.linspace(1,np.floor(0.5 * len(pn)),np.floor(0.5 * len(pn))) / (len(pn) * Delta_t)
    else:
        freq_single = np.linspace(1,np.floor(0.5 * len(pn)) - 1,np.floor(0.5 * len(pn)) - 1) / (len(pn) * Delta_t)
    PNT_single = gls(pn,t,freq_single)

    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')
    plt.plot(freq_single,PNT_single,'b')
    plt.ylabel('$P(\\nu)$')
    plt.xlabel('$\\nu$')
    plt.xlim(0,freq_single[-1])
    plt.ylim(0,np.max(PNT_single) + 10.0)
    plt.title('LS spectrum of residuals')
    plt.margins(0.2)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig(os.path.join(path_tot,'spectrum_fin.pdf'))
    plt.close()

    m_res = np.nanmean(pn)
    s_res = np.nanstd(pn)
    pn_norm = (pn - m_res) / s_res
    pn_norm_notnan = pn_norm[~np.isnan(pn_norm)]
    outlier_lim = 3.0
    num_outliers_max = len(pn_norm_notnan[pn_norm_notnan > outlier_lim])
    num_outliers_min = len(pn_norm_notnan[pn_norm_notnan < -outlier_lim])
    num_outliers = num_outliers_max + num_outliers_min

    pos_out = []
    path_file = os.path.join(path_tot,'outliers.txt')
    f = open(path_file,'w')
    for i in range(len(t)):
        if pn_norm[i] < -outlier_lim or pn_norm[i] > outlier_lim:
            pos_out.append(i)
            f.write('%d %f %s\n' % (t[i],pn_norm[i],time[i]))
    f.close()

    plt.figure(figsize = (12,9))
    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')
    
    plt.subplot(2,1,1)
    plt.plot(t,pn)
    plt.xlim(t[0],t[-1])
    plt.xticks(ticks,'')
    plt.ylabel('$N_t$')
    plt.title('residuals / normalised residuals')
    plt.margins(0.2)
    plt.subplots_adjust(hspace = 0.0)

    plt.subplot(2,1,2)
    s_res = '%.2f' % s_res
    if int(matplotlib.__version__.split('.')[0]) == 2:
        plt.bar(t,pn_norm,width = 10,label = 'num outl = ' + str(num_outliers))
    else:
        plt.bar(t,pn_norm,width = 0.1,label = 'num outl = ' + str(num_outliers))
    plt.plot((t[0],t[-1]),(outlier_lim,outlier_lim),'r',label = '$\sigma$ = ' + str(s_res))
    plt.plot((t[0],t[-1]),(-outlier_lim,-outlier_lim),'r')
    plt.legend(loc = 0)
    plt.xlim(t[0],t[-1])
    plt.xticks(ticks,time_label,rotation = 'vertical')
    plt.ylabel('$N_t^{norm}$')
    plt.margins(0.2)
    plt.subplots_adjust(hspace = 0.0)

    plt.savefig(os.path.join(path_tot,'res.pdf'))
    plt.close()

    path_file = os.path.join(path_tot,'res.txt')
    f = open(path_file,'w')
    for i in range(len(t)):
        f.write('%f %f %f %s\n' % (t[i],pn[i],pn_norm[i],time[i]))
    f.close()

    return pos_out
#############################################################################
#############################################################################
########################### AUTOCORRELATION #################################
###################### last modified 05/07/2016 #############################
### Computes autocorrelation and plots it
### INPUT ARGUMENTS
### pn -> data vector
### path_tot -> path to the main folder
#############################################################################
def autocorrelation(pn,path_tot):

    max_lag = 20
    pn[np.isnan(pn)] = 0
    pn -= np.mean(pn)

    autocorr = np.correlate(pn,pn,mode = 'full')
    autocorr /= autocorr[np.argmax(autocorr)]
    half = autocorr[autocorr.size / 2:]
    bounds = [2.0 / np.sqrt(len(pn)),-2.0 / np.sqrt(len(pn))]

    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')
    markerline,stemlines,baseline = plt.stem(np.arange(0,max_lag + 1),half[0:max_lag + 1],'-')
    plt.setp(baseline,'color','r')
    plt.setp(stemlines,'color','r')
    plt.setp(markerline,'color','r',markeredgecolor = 'r')
    plt.plot((-1,max_lag + 1),(0,0),'r')
    plt.plot((-1,max_lag + 1),(bounds[0],bounds[0]),'b')
    plt.plot((-1,max_lag + 1),(bounds[1],bounds[1]),'b')
    plt.xlim(-1,max_lag + 1)
    plt.ylim(np.min(half) - 0.2,1.2)
    plt.ylabel('$R(\\tau)$')
    plt.xlabel('$\\tau$')
    plt.title('Autocorrelation of residuals')
    plt.savefig(os.path.join(path_tot,'autocorrelation_residuals.pdf'))
    plt.close()

    path_file = os.path.join(path_tot,'autocorrelation_residuals.txt')
    f = open(path_file,'w');
    for i in np.arange(0,max_lag + 1):
        f.write('%d %f\n' % (i,half[i]))
    f.close()
#############################################################################
#############################################################################
################################# DFA #######################################
###################### last modified 12/03/2017 #############################
### Computes detrended fluctuations analysis coefficient and plots it
### INPUT ARGUMENTS
### pn -> data vector
### min_win -> minimum scale
### rev_seg -> windows forward and backward
### path_tot -> path to the main folder
### OUTPUT ARGUMENTS
### H_mono -> Hurst exponent
#############################################################################
def dfa(pn,min_win,rev_seg,path_tot):

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
        F[i] = np.sqrt(np.nansum(F_nu) / float(len(F_nu[~np.isnan(F_nu)])))

    log_fit = np.polyfit(np.log(s),np.log(F),1)
    H_mono = '%.2f' % log_fit[0]
    log_fit = np.poly1d(log_fit)
    log_fit = log_fit(np.log(s))

    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')
    plt.plot(np.log(s),np.log(F),'o',label = '$H$ = ' + H_mono)
    plt.plot(np.log(s),log_fit,'r')
    plt.legend(loc = 0)
    plt.xlim(np.log(s[0]) - 0.3,np.log(s[-1]) + 0.3)
    plt.ylabel('log$(F(n))$')
    plt.xlabel('log$(n)$')
    plt.title('DFA fit')
    plt.margins(0.2)
    plt.subplots_adjust(bottom = 0.2)
    plt.savefig(os.path.join(path_tot,'dfa.pdf'))
    plt.close()

    path_file = os.path.join(path_tot,'dfa.txt')
    f = open(path_file,'w');
    for i in range(len(s)):
        f.write('%f %f %f\n' % (s[i],F[i],log_fit[i]))
    f.close()

    return H_mono
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
### path_tot -> path to the main folder
#############################################################################
def mdfa(H_mono,pn,min_win,q_max,rev_seg,path_tot):
    
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

    plt.figure(figsize = (11,11))
    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')

    plt.subplot(2,2,1)
    plt.plot(np.log(s),np.log(F[0,:]),'b.')
    plt.plot(np.log(s),MDFA_fit[:,0],'b',label = 'q = -3')
    plt.plot(np.log(s),np.log(F[50,:]),'r.')
    plt.plot(np.log(s),MDFA_fit[:,50],'r',label = 'q = 0')
    plt.plot(np.log(s),np.log(F[-1,:]),'g.')
    plt.plot(np.log(s),MDFA_fit[:,-1],'g',label = 'q = 3')
    plt.legend(loc = 0)
    plt.xlim(np.log(s[0]),np.log(s[-1]))
    plt.xlabel('log(n)')
    plt.ylabel('log(F(n))')
    plt.title('MDFA fit')
    plt.margins(0.2)
    plt.subplots_adjust(bottom = 0.2)

    plt.subplot(2,2,2)
    plt.plot(q,H,'b',label = 'h(q)')
    plt.plot((q[0],q[-1]),(H_mono,H_mono),'k',label = 'H')
    plt.legend(loc = 0)
    plt.xlim(q[0],q[-1])
    plt.xlabel('q')
    plt.ylabel('h(q)')
    plt.title('Generalised Hurst exponent')
    plt.margins(0.2)
    plt.subplots_adjust(bottom = 0.2)

    plt.subplot(2,2,3)
    plt.plot(alpha,sing_spec,'b')
    plt.ylim(min(sing_spec) - 0.2,1.2)
    plt.xlabel('$\\alpha$')
    plt.ylabel('$f(\\alpha)$')
    plt.title('Singularity spectrum')
    plt.margins(0.2)
    plt.subplots_adjust(bottom = 0.2)

    plt.savefig(os.path.join(path_tot,'mdfa.pdf'))
    plt.close()

    path_file = os.path.join(path_tot,'mdfa1.txt')
    f = open(path_file,'w');
    for i in range(len(s)):
        f.write('%f %f %f %f %f %f %f\n' % (s[i],F[0,i],MDFA_fit[i,0],F[50,i],MDFA_fit[i,50],F[-1,i],MDFA_fit[i,-1]))
    f.close()

    path_file = os.path.join(path_tot,'mdfa2.txt')
    f = open(path_file,'w');
    for i in range(len(q)):
        f.write('%f %f %f\n' % (q[i],H[i],tau[i]))
    f.close()

    path_file = os.path.join(path_tot,'mdfa3.txt')
    f = open(path_file,'w');
    for i in range(len(alpha)):
        f.write('%f %f\n' % (alpha[i],sing_spec[i]))
    f.close()
#############################################################################
#############################################################################
################################ MFDFA2 #####################################
###################### last modified 23/03/2017 #############################
### Computes Local Hurst exponent (python version of Ihlen's Matlab MFDFA2.m
### INPUT ARGUMENTS
### signal -> data vector
### scale -> scales at which compute H(t)
### m -> order for fluctuations polynomial fit
### path_tot -> path to the main folder
#############################################################################
def MFDFA2(signal,scale,m,path_tot):

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
        segments = len(X) / scale0[ns]
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

    plt.figure(figsize = (12,9))
    plt.rc('text',usetex = True)
    plt.rc('font',family = 'serif')

    Ht_plot = Ht[0,:]
    plt.subplot(2,1,1)
    ax = plt.gca()
    if int(matplotlib.__version__.split('.')[0]) == 2:
        ax.set_facecolor('black')
    else:
        ax.set_axis_bgcolor('black')
    plt.plot(Ht_plot,'y')
    plt.plot(0.5 * np.ones((len(Ht_plot),)),'w')
    plt.plot(np.ones((len(Ht_plot),)),'m')
    plt.plot(1.5 * np.ones((len(Ht_plot),)),'r')
    plt.xlim(0,len(Ht_plot))
    plt.ylim(0,3)
    plt.xlabel('time')
    plt.ylabel('$H_t$')
    plt.title('local Hurst exponent')
    plt.margins(0.2)
    plt.subplots_adjust(hspace = 0.3)

    plt.subplot(2,1,2)
    plt.plot(Htbin,Ph,'b',label = '$\mu$ = ' + mu)
    plt.plot(Htbin,fit_gauss,'r',linewidth = 2.0,label = '$\sigma$ = ' + sigma)
    plt.legend(loc = 0)
    plt.ylim(0,)
    plt.ylabel('Ph')
    plt.xlabel('$H_t$')
    plt.title('prob distr of $H_t$')
    plt.margins(0.2)
    plt.subplots_adjust(hspace = 0.3)

    plt.savefig(os.path.join(path_tot,'MFDFA2.pdf'))
    plt.close()

    path_file = os.path.join(path_tot,'Ht.txt')
    f = open(path_file,'w');
    for i in range(len(Ht_plot)):
        f.write('%f\n' % (Ht_plot[i]))
    f.close()

    path_file = os.path.join(path_tot,'Ht_distr.txt')
    f = open(path_file,'w');
    for i in range(len(Htbin)):
        f.write('%f %f %f %s %s\n' % (Htbin[i],Ph[i],fit_gauss[i],mu,sigma))
    f.close()
#############################################################################
#############################################################################
########################## DISTRIBUTIONS_FIT ################################
###################### last modified 05/07/2016 #############################
### Computes histogram of residuals and fits it with known distributions
### INPUT ARGUMENTS
### pn -> data vector
### pos_out -> positions of outliers (0 if not needed)
### path_tot -> path to the main folder
### type_hist -> initial, noOutliers, Outliers (for plot's labels)
#############################################################################
def distributions_fit(pn,pos_out,path_tot,type_hist):

    path_tot = os.path.join(path_tot,'distributions')
    if not os.path.exists(path_tot):
        os.makedirs(path_tot)

    if not np.isscalar(pos_out):
        pn[pos_out] = np.nan
    pn = pn[~np.isnan(pn)]
    if not np.any(pn < 0.0):
        f_hist = open(path_tot + '/distributions_fit_param_' + type_hist + '.txt','w')
        f_hist.write('Distribution fit_param1  fit_param2   fit_param3   fit_param4   test_ks  p_ks   shift\n')
        if int(np.__version__.split('.')[1]) >= 11:
            BinNumb = 'fd'
        else:
            BinNumb = int(np.sqrt(len(pn)))
        freq,Htbin = np.histogram(pn,bins = BinNumb)
        Htbin = (Htbin[:-1] + Htbin[1:]) / 2.0

        param = burr.fit(pn)
        fit_hist = burr.pdf(Htbin,param[0],param[1],loc = param[2],scale = param[3])
        test_value,p_value = kstest(pn,'burr',args = (param[0],param[1],param[2],param[3]))
        f_hist.write('Burr         %-11.3f %-12.3f %-12.3f %-12.3f %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],param[3],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Burr')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Burr_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = expon.fit(pn)
        fit_hist = expon.pdf(Htbin,loc = param[0],scale = param[1])
        test_value,p_value = kstest(pn,'expon',args = (param[0],param[1]))
        f_hist.write('Exponential  %-11.3f %-12.3f NaN          NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Exponential')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Exponential_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = gamma.fit(pn)
        fit_hist = gamma.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'gamma',args = (param[0],param[1],param[2]))
        f_hist.write('Gamma        %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Gamma')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Gamma_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = genpareto.fit(pn)
        fit_hist = genpareto.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'genpareto',args = (param[0],param[1],param[2]))
        f_hist.write('GenPareto    %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Generalised Pareto')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'GenPareto_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = invgauss.fit(pn)
        fit_hist = invgauss.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'invgauss',args = (param[0],param[1],param[2]))
        f_hist.write('InvGauss     %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Inverse Gaussian')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'InvGauss_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = lognorm.fit(pn)
        fit_hist = lognorm.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'lognorm',args = (param[0],param[1],param[2]))
        f_hist.write('Lognormal    %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Lognormal')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Lognormal_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = norm.fit(pn)
        fit_hist = norm.pdf(Htbin,loc = param[0],scale = param[1])
        test_value,p_value = kstest(pn,'norm',args = (param[0],param[1]))
        f_hist.write('Normal       %-11.3f %-12.3f NaN          NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Normal')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Normal_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = rayleigh.fit(pn)
        fit_hist = rayleigh.pdf(Htbin,loc = param[0],scale = param[1])
        test_value,p_value = kstest(pn,'rayleigh',args = (param[0],param[1]))
        f_hist.write('Rayleigh     %-11.3f %-12.3f NaN          NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Rayleigh')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Rayleigh_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = t.fit(pn)
        fit_hist = t.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'t',args = (param[0],param[1],param[2]))
        f_hist.write('T            %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : T')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'T_res_hist_' + type_hist + '.pdf'))
        plt.close()
    
        param = weibull_max.fit(pn)
        fit_hist = weibull_max.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'weibull_max',args = (param[0],param[1],param[2]))
        f_hist.write('Weib_max     %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Weibull max')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Weib_max_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = weibull_min.fit(pn)
        fit_hist = weibull_min.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'weibull_min',args = (param[0],param[1],param[2]))
        f_hist.write('Weib_min     %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Weibull min')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Weib_min_res_hist_' + type_hist + '.pdf'))
        plt.close()

        f_hist.close()
    else:
        f_hist = open(path_tot + '/distributions_fit_param_' + type_hist + '.txt','w')
        f_hist.write('Distribution fit_param1  fit_param2   fit_param3   fit_param4   test_ks  p_ks   shift\n')
        if int(np.__version__.split('.')[1]) >= 11:
            BinNumb = 'fd'
        else:
            BinNumb = int(np.sqrt(len(pn)))
        freq,Htbin = np.histogram(pn,bins = BinNumb)
        Htbin = (Htbin[:-1] + Htbin[1:]) / 2.0

        param = norm.fit(pn)
        fit_hist = norm.pdf(Htbin,loc = param[0],scale = param[1])
        test_value,p_value = kstest(pn,'norm',args = (param[0],param[1]))
        f_hist.write('Normal       %-11.3f %-12.3f NaN          NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Normal')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Normal_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = t.fit(pn)
        fit_hist = t.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn,'t',args = (param[0],param[1],param[2]))
        f_hist.write('T            %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f NaN\n' % (param[0],param[1],param[2],test_value,p_value))
        plt.hist(pn,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : T')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'T_res_hist_' + type_hist + '.pdf'))
        plt.close()

        shift = 2.0 * np.floor(np.min(pn))
        pn2 = pn - shift
        freq,Htbin = np.histogram(pn2,bins = BinNumb)
        Htbin = (Htbin[:-1] + Htbin[1:]) / 2.0

        param = burr.fit(pn2)
        fit_hist = burr.pdf(Htbin,param[0],param[1],loc = param[2],scale = param[3])
        test_value,p_value = kstest(pn2,'burr',args = (param[0],param[1],param[2],param[3]))
        f_hist.write('Burr         %-11.3f %-12.3f %-12.3f %-12.3f %-8.3f %-6.3f %d\n' % (param[0],param[1],param[2],param[3],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Burr')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Burr_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = expon.fit(pn2)
        fit_hist = expon.pdf(Htbin,loc = param[0],scale = param[1])
        test_value,p_value = kstest(pn2,'expon',args = (param[0],param[1]))
        f_hist.write('Exponential  %-11.3f %-12.3f NaN          NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Exponential')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Exponential_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = gamma.fit(pn2)
        fit_hist = gamma.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn2,'gamma',args = (param[0],param[1],param[2]))
        f_hist.write('Gamma        %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],param[2],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Gamma')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Gamma_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = genpareto.fit(pn2)
        fit_hist = genpareto.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn2,'genpareto',args = (param[0],param[1],param[2]))
        f_hist.write('GenPareto    %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],param[2],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Generalised Pareto')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'GenPareto_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = invgauss.fit(pn2)
        fit_hist = invgauss.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn2,'invgauss',args = (param[0],param[1],param[2]))
        f_hist.write('InvGauss     %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],param[2],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Inverse Gaussian')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'InvGauss_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = lognorm.fit(pn2)
        fit_hist = lognorm.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn2,'lognorm',args = (param[0],param[1],param[2]))
        f_hist.write('Lognormal    %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],param[2],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Lognormal')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Lognormal_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = rayleigh.fit(pn2)
        fit_hist = rayleigh.pdf(Htbin,loc = param[0],scale = param[1])
        test_value,p_value = kstest(pn2,'rayleigh',args = (param[0],param[1]))
        f_hist.write('Rayleigh     %-11.3f %-12.3f NaN          NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Rayleigh')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Rayleigh_res_hist_' + type_hist + '.pdf'))
        plt.close()

        param = weibull_max.fit(pn2)
        fit_hist = weibull_max.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn2,'weibull_max',args = (param[0],param[1],param[2]))
        f_hist.write('Weib_max     %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],param[2],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Weibull max')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Weib_max_res_hist_' + type_hist + '.pdf'))
        plt.close()
        
        param = weibull_min.fit(pn2)
        fit_hist = weibull_min.pdf(Htbin,param[0],loc = param[1],scale = param[2])
        test_value,p_value = kstest(pn2,'weibull_min',args = (param[0],param[1],param[2]))
        f_hist.write('Weib_min     %-11.3f %-12.3f %-12.3f NaN          %-8.3f %-6.3f %d\n' % (param[0],param[1],param[2],test_value,p_value,shift))
        plt.hist(pn2,BinNumb,normed = 1)
        plt.plot(Htbin,fit_hist,'r',linewidth = 2.0)
        plt.ylim(0,)
        plt.ylabel('$H(N_t)$')
        plt.xlabel('$N_t$')
        plt.title('Fitted distribution : Weibull min')
        plt.margins(0.2)
        plt.savefig(os.path.join(path_tot,'Weib_min_res_hist_' + type_hist + '.pdf'))
        plt.close()

        f_hist.close()
#############################################################################
