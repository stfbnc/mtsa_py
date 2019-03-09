# -*- coding: iso-8859-15 -*-
### last modified 01/04/2017

import MTsA_fun as MF
import warnings
warnings.filterwarnings("ignore")

file_name,time_col,data_col,Delta_t,time_units,ofac,typeoffit,year_in,year_fin,rev_seg,scale_min,scale_MFDFA = MF.args_parser()

print '\n\n        MTsA (v 2.13.1)        \n\n'

### check on inputs
units_converter = MF.input_checks(time_col,data_col,Delta_t,year_in,year_fin,ofac,time_units,typeoffit,rev_seg,scale_min,scale_MFDFA)

### main folder
path_tot = MF.main_folder(file_name,data_col)

### load file
pn,time,nan_pos_pn,t,time_label,ticks = MF.load_file(path_tot,file_name,data_col,time_col,Delta_t)

### detrending
pn = MF.trend_detrend(time_col,time,pn,t,typeoffit,path_tot,time_label,ticks)

### initial distributions
MF.distributions_fit(pn,0,path_tot,'initial')

### spectrum
PNT,freq,pth = MF.lomb_scargle(pn,t,ofac,Delta_t,path_tot)

### find all peaks iteratively
MF.peaks_filter_all(len(freq),path_tot,PNT,freq,pth,units_converter,year_in,year_fin,t,nan_pos_pn,Delta_t,pn,ofac)

### filter
freq_fig,pn,sig_phi,sig_to_noise = MF.peaks_filter(1,path_tot,PNT,freq,pth,units_converter,year_in,year_fin,t,nan_pos_pn,Delta_t,pn,ofac)

### residuals spectrum
pos_out = MF.res_lomb(pn,Delta_t,t,time,ticks,time_label,pth,path_tot)

### filtered signal figures
MF.phase_fig(sig_phi,time_col,time,time_label,ticks,freq_fig,path_tot,t)

### detrended fluctuation analysis
dfa_coeff = MF.dfa(pn,scale_min,rev_seg,path_tot)

### multifractal detrended fluctuation analysis
MF.mdfa(dfa_coeff,pn,scale_min,3.0,rev_seg,path_tot)

### local Hurst exponent
MF.MFDFA2(pn,scale_MFDFA,1,path_tot)

### distributions with outliers
MF.distributions_fit(pn,0,path_tot,'Outliers')

### distributions without outliers
MF.distributions_fit(pn,pos_out,path_tot,'noOutliers')

print '\n\n        END OF THE ANALYSIS        \n\n'