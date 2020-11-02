#!/usr/bin/env python3
# --*-- coding:utf-8 --*--

"""
This module reads a spike list form a csv file in the form of time-stamps and,
given time gates T1, .. Tn, it calculates the mean number of spikes per time gate 
and its variance. 

ToDo:
* Replace input parsing with argparse
* Output raw histograms


Things to note:
* Output file is showing Ts with 3 decimal points. Should it be more?
* For Axion's spike_list.csv analysis: the data file contains a 9 lines footer, which this script ignores. Update that if anything changes.
"""

# Owned
__author__ = "Or Nahmani"
__email__ = "ornahm@post.bgu.ac.il"
__company__ = "NIBN, Ben-Gurion University of the Negev"
__date__ = "20/03/2020"
__license__ = "MIT License"
__copyright__ = "Copyright 2020, NIBN, Israel"
__status__ = "Dev"
__credits__ = [__author__]
__version__ = "0.2.0"
   
import numpy as np
import sys
import csv
import re
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import math
import copy
import os
from collections import defaultdict, OrderedDict
from scipy import stats
import tkinter as tk
from tkinter import filedialog
import glob
import argparse
import xlsxwriter


#### AXION PARSING ####
# spikes list file fields: ",,Time (s),Electrode,Amplitude (mV)"
AMPLITUDE_IND = 4
ELECTRODE_IND = 3
TIMESTAMP_IND = 2
N_PLATE_TYPE_24 = 4 # 16 elecs
N_PLATE_TYPE_6 = 8# 64 elecs
REC_TIME_RE = '.* to ((?P<hr>\\d+)h)?((?P<min>\\d+)m)?(?P<sec>\\d+)s'
FOOTER_RE = '^Well Information.*'
LAST_DATA_LINE_RE = '.*,(?P<elec>[A-Z]\d_\d{2}),.*'
REC_TIME_LINE_RE = '.*Actual File Section Run*'
PLATE_TYPE_LINE_RE = '\s+Plate Type'
NUM_WELLS_RE = '.*MEA (?P<ptype>\\d+)'

def parse_axion_spikes_list(path):
    # Gets path of axion's _spike_list.csv file and returns a dict of WellxElectrode: spikes ts array
    try:
        with open(path) as spikesfile:
            # load spikes
            print("...reading file")
            lines = spikesfile.readlines();
            print("...done")
 
        # start parsing file
        print("...parsing file")

        # find header - Well Information
        iline =-1
        while re.match(FOOTER_RE, lines[iline]) == None:
            iline -= 1
            if abs(iline) == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND FOOTER (STRING: Well Information)\n')            
        data_wells = lines[iline+1].rstrip().split(",")[1:]
        data_colors = lines[iline+3].rstrip().replace("#","").split(",")[1:]

        # find last data line
        iline =-1
        while re.match(LAST_DATA_LINE_RE, lines[iline]) == None and abs(iline) < len(lines):
            iline -= 1
        data = np.loadtxt(lines[1:iline+1], dtype="str", delimiter=",")
        
        # parse the recording's length
        iline = 0
        while re.match(REC_TIME_LINE_RE, lines[iline]) == None:
            iline += 1
            if iline == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND RECORDING TIME (REQUESTED/ACTUAL DURATION OF MEASUREMENT)\n')                   
        time_match = re.match(REC_TIME_RE, lines[iline]).groupdict()
        
        rec_time = 0
        if time_match['hr']:
            rec_time += int(time_match['hr'])*3600
        if time_match['min']:
            rec_time += int(time_match['min'])*60
        if time_match['sec']:
            rec_time += int(time_match['sec']) 

        # parse the plate type - CytoView MEA 6/24
        iline = 0
        while re.match(PLATE_TYPE_LINE_RE, lines[iline]) == None:
            iline += 1
            if iline == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND PLATE TYPE (CytoView MEA 6/24)\n')            
        plate_type = int(re.match(NUM_WELLS_RE,lines[iline]).group('ptype'))

        # Parse footer (colors for sort)
        wells_color = defaultdict(list)   
        for i, j in zip(data_colors, data_wells):
            wells_color[i].append(j)
        wells_color = OrderedDict(sorted(wells_color.items())) # sort dict

        for i in wells_color.keys(): # sort list of values per key
            wells_color[i].sort(key=lambda s: s[::-1])

        electrodes = set(data[:, ELECTRODE_IND])        
        
        # Assuming all wells are the same (containing same electrodes alignment)
        elecs_per_well = {elec[-2:] for elec in electrodes}
        wells = {elec[:2] for elec in electrodes}
        
        # Initialize spikes dict for all electrodes
        elecs = {elec:[] for elec in elecs_per_well}
        #import pdb; pdb.set_trace()
        spikes_lists = {well: copy.deepcopy(elecs) for well in wells}
        amps_lists = {well: copy.deepcopy(elecs) for well in wells}
        for row in data:
            well, elec = row[ELECTRODE_IND].split('_')
            ts = float(row[TIMESTAMP_IND])
            amp = float(row[AMPLITUDE_IND])
            spikes_lists[well][elec].append(ts)
            amps_lists[well][elec].append(amp)

        # Convert to numpy arrays for further proccessing. this is not done initially because np.append copies the array.
        for well in wells:
            spikes_lists[well] = {elec:np.array(sp_list) for elec,sp_list in spikes_lists[well].items() if sp_list!=[]}
            amps_lists[well] = {elec:np.array(amp_list) for elec,amp_list in amps_lists[well].items() if amp_list!=[]}
        
        print("...done")
        
    except IOError:
        print("File not accessible")
    
    # pack returns
    spk_data = spikes_lists, wells, elecs_per_well,\
        rec_time, amps_lists, wells_color, plate_type
        
    return spk_data


def parse_mcs_spikes_list(path):
    # Gets path of multi channel system's spike list (exported from Multi Channel Analyzer) and returns a dict of Electrode: spikes ts array
    try:
        with open(path) as spikesfile:
            data = np.loadtxt(spikesfile.readlines()[3:], dtype="str", delimiter="\t")
        
        wells = [''] # A single well
        # Set the sample time to be the last spike's timestamp
        rec_time = 0
        spikes_lists = {}
        for col in data.T:
            elec = col[0][-2:]
            # timestamps are in µs. convert to seconds
            ts_list = [float(ts)/(10**6) for ts in col[1:] if ts!='']
            spikes_lists[elec] =  np.array(ts_list)
            rec_time = max(rec_time, max(ts_list))
        
    except IOError:
        print("File not accessible")
    
    return {'':spikes_lists}, wells, rec_time


# gets parsed spikeslist dict,returns a dict structure of TxWELLxELECTRODExSPT and plots the histograms
def spikes_per_T(elec_spikes_dict, wells, electrodes, Ts, sample_time):
    #
    # Ts - time durations (s)
    # sample_time - total sample time (s)
    print("...calculating spikes per T")
    data = {}
    # configure histograms figure 
    #ax_size = math.ceil(math.sqrt((len(electrodes)-1))) # num of rows/cols in the MEA grid

    for T in Ts:
        data[T] = {}
        bins = np.arange(0, sample_time+1, T)
        for well in wells:
            well_data = {elec: np.histogram(spikes, bins)[0] for elec, spikes in elec_spikes_dict[well].items()}

            # if T in plotting_Ts:
            #     # plot histograms
            #     fig, axs = plt.subplots(ax_size,ax_size, figsize=(25,15))
            #     fig.suptitle("Well %s, T = %.3f" %(well, T))
            #     for elec, spt in well_data.items():
            #         x, y = int(elec[0])-1, int(elec[1])-1
            #         axs[x, y].set_title(elec)
            #         axs[x, y].plot(bins[:-1], spt)
            #     #plt.savefig(os.path.join(output_dir,"hist_%s_%.3f.png" %(well,T)))
            #     plt.close(fig)

            # Add total well's spt as well
            well_data[""] = sum(well_data.values())
            data[T][well] = well_data
        
    print("...done")
    
    return data 


# gets parsed spikeslist dict,returns a dict structure of TxWELLxELECTRODExSPT and plots the histograms
def amps_per_T(elec_amps_dict, elec_spikes_dict, wells, electrodes, Ts, sample_time):
    #
    # Ts - time durations (s)
    # sample_time - total sample time (s)
    print("...calculating amplitudes per T")
    data = {}
    # configure histograms figure 
    #ax_size = math.ceil(math.sqrt((len(electrodes)-1))) # num of rows/cols in the MEA grid

    for T in Ts:
        data[T] = {}
        bins = np.arange(0, sample_time+1, T)
        for well in wells:
            well_data = {}
            for elec in elec_amps_dict[well].keys():
                spikes = elec_spikes_dict[well][elec] # time stamps of spikes
                amps   = elec_amps_dict[well][elec]
                well_data[elec] = np.nan_to_num(stats.binned_statistic(spikes, amps, statistic='mean', bins=bins)[0])                
                #well_data = {elec: np.histogram(spikes, bins)[0] for elec, spikes in elec_amps_dict[well].items()}

            # if T in plotting_Ts:
            #     # plot histograms
            #     fig, axs = plt.subplots(ax_size,ax_size, figsize=(25,15))
            #     fig.suptitle("Well %s, T = %.3f" %(well, T))
            #     for elec, spt in well_data.items():
            #         x, y = int(elec[0])-1, int(elec[1])-1
            #         axs[x, y].set_title(elec)
            #         axs[x, y].plot(bins[:-1], spt)
            #     #plt.savefig(os.path.join(output_dir,"hist_%s_%.3f.png" %(well,T)))
            #     plt.close(fig)

            # Add mean well's amp as well            
            ave_amp = np.zeros(len(bins)-1)
            for j in range(len(bins)-1):
                tot_amp = 0.0
                icount = 0
                for elec in well_data.keys():
                    if well_data[elec][j] > 0:
                        tot_amp += well_data[elec][j]
                        icount += 1
                if icount > 0: 
                    ave_amp[j] = tot_amp/icount                                     
                    
            well_data[""] = ave_amp
            
            data[T][well] = well_data
    
    print('...done')    
    
    return data 


def plot_heatmaps(data, well, T, labels, path):
    # data = 2D np array of single well's electrodes spikes per T
    fig, axs = plt.subplots(1, 2, figsize=(15,10))
    fig.suptitle("Well %s, T = %.3f" %(well, T))
    cov = np.cov(data)
    cov_g = sns.heatmap(cov, center=0, cmap=sns.diverging_palette(20, 220, n=200), ax=axs[0], square=True)
    cov_g.set_xticklabels(labels, rotation=45, horizontalalignment='right');
    cov_g.set_yticklabels(labels, rotation=0);
    cov_g.set_title('Electrodes Covariance')
    
    corrcoef = np.corrcoef(data)
    corrcoef_g = sns.heatmap(corrcoef, vmin=-1, vmax=1, center=0, cmap=sns.diverging_palette(20, 220, n=200), ax=axs[1], square=True)
    corrcoef_g.set_xticklabels(labels, rotation=45, horizontalalignment='right');
    corrcoef_g.set_yticklabels(labels, rotation=0);
    corrcoef_g.set_title('Electrodes Corrcoef')

    #plt.savefig(path)
    plt.close(fig)


def save_mean_vars(yfuncs_output_csv, Ts, output_path):
    output_title = ["Well/Electrode"] + ','.join(["T=%.3f mean,T=%.3f variance,T=%.3f var/mean" %(T,T,T) for T in Ts]).split(',')
    with open(output_path, 'w', newline='') as outfile:
        w = csv.writer(outfile, delimiter=',')
        w.writerow(output_title)
        # write ordered: wells first then electrodes.
        for item in sorted(sorted(yfuncs_output_csv.keys()), key=len):
            w.writerow([item] + yfuncs_output_csv[item])


def get_elecs_full_list(well, plate_type):
    n = N_PLATE_TYPE_24 if plate_type == 24 else N_PLATE_TYPE_6
    
    elec_list = []
    for i in range(1,n+1):
        for j in range(1,n+1):
            elec_list += [str(10*i+j)]
    
    return elec_list
    

def save_raw(data, Ts, ofname, workbook, wells_color, plate_type):
        
    for T in Ts:
        print("...saving file: %s" % ofname)
        with open(ofile, 'w', newline='') as outfile:
            w = csv.writer(outfile, delimiter=',')
            
            # span time axis
            exist_well = [*data[T].keys()][0]
            exist_elec = [*data[T][exist_well].keys()][0]
            nbins = len(data[T][exist_well][exist_elec])
            xtime = np.linspace(T,nbins*T,nbins)          
            w.writerow(['time'] + xtime.tolist())
            
            # write ordered: wells first then wells+electrodes.
            for color in wells_color:
                for well in wells_color[color]:     
                    if well in data[T]:
                        electrodes = data[T][well]
                            
                        # print only wells first
                        for elec in sorted(electrodes.keys()):  
                            if (elec == ''):
                                line = electrodes[elec].tolist()
                                ls = [x if x > 0 else "" for x in line]
                                w.writerow([well] + ls)
                    else:
                        ls = ['']*nbins
                        w.writerow([well] + ls)
                        
            w.writerow('') # empty line

            # write wells+electrodes.
            for color in wells_color:                                
                for well in wells_color[color]:
                    if well in data[T]:
                        electrodes = data[T][well]
                        
                        # print all        
                        elecs_full_list = get_elecs_full_list(well, plate_type)
                        for elec in elecs_full_list:
                            ttl = well if (elec == '') else (well + '_' + elec)
                            if elec in electrodes.keys():                                                           
                                line = electrodes[elec].tolist()
                                ls = [x if x > 0 else "" for x in line]                                
                            else:
                                ls = ['']*nbins
                            w.writerow([ttl] + ls)
                    else:
                        ls = ['']*nbins
                        w.writerow([well] + ls)
                        elecs_full_list = get_elecs_full_list(well, plate_type)
                        for elec in elecs_full_list:
                            ttl = well if (elec == '') else (well + '_' + elec)
                            ls = ['']*nbins
                            w.writerow([ttl] + ls)

                        # for elec in sorted(electrodes.keys()):  
                        #     ttl = well if (elec == '') else (well + '_' + elec)
                        #     line = electrodes[elec].tolist()
                        #     ls = [x if x > 0 else "" for x in line]
                        #     w.writerow([ttl] + ls)
                            
    print("...done")
    
# gets window size, and plots yfuncs and cov heatmaps for the rolling time window
def plot_time_windows(wells, Ts, plotting_Ts, data,output_dir, window_size=0):
    # Create yfuncs_output_csv - a dict of elec/well: list of mean,variance,mean,variance... for each given T in list, and in that order.
    yfuncs_output_csv = dict()
    wells_yfuncs = {well:[] for well in wells}
    for T in Ts:
        for well in wells:
            well_T_data = data[T][well]
            #update mean&var data - for the csv output
            for elec, spt in well_T_data.items():
                mean, var = spt.mean(), spt.var()
                print("%s  %s  %f  %f" % (well,elec,mean,var))
                key = well+"_" + elec if elec else well
                if key not in yfuncs_output_csv.keys():
                    yfuncs_output_csv[key] = []
                yfuncs_output_csv[key].extend([mean, var, var/mean])

            # Add well's yfunc (=var/mean) - for plotting
            wells_yfuncs[well].append(well_T_data[""].var()/well_T_data[""].mean())

            # plot covariance and corrcoef matrices
            if T in plotting_Ts and len(well_T_data.keys())>2:  # plot if we have at least 2 active electrodes. >2 because the well ("") is also a key
                electrodes = sorted(well_T_data.keys())
                electrodes.remove("")
                # Stack active electrodes in order (electrodes is a sorted list)
                spt_stack = np.stack([well_T_data[elec] for elec in electrodes], axis=0)
                plot_heatmaps(spt_stack, well, T, electrodes, os.path.join(output_dir,"covs_%s_%.3f.png" %(well,T)))

    # Save yfuncs_output_csv to file
    # save_mean_vars(yfuncs_output_csv, Ts, os.path.join(output_dir,"fano_factors.csv"))
    # Plot yfuncs
    plt.figure(figsize=(10,10))
    for well in sorted(wells_yfuncs.keys()):
        plt.plot(Ts, wells_yfuncs[well], label=well, marker='.')
    plt.legend()
    #plt.xlim(0, Ts[-1]+0.5)
    plt.xticks(Ts)
    plt.xlabel('T')
    plt.ylabel('var(T)/mean(T)')
    #plt.savefig(os.path.join(output_dir,"fano_factors.png"))


def analyze_file(input_type, Ts, input_file):

    print("...analyzing file: %s" % input_file)
    
    parsing_func = parse_axion_spikes_list if input_type=='axion' else parse_mcs_spikes_list
    spk_data = parsing_func(input_file)    
    
    spikes_dict, wells, electrodes, sample_time, amps_dict,\
        wells_color, plate_type = spk_data
        
    # check time gate
    print('...sample time is: %d s' % sample_time)
    if Ts[0] > sample_time:
        #sys.exit('\n>>>>> ERROR: TIME GATE (T) IS LARGER THAN THE SAMPLE LENGTH\n')
        print('...WARNING: time gate (T=%d s) is larger than the sample length (%d s)' % (Ts[0], sample_time))
        print('...skipping analysis of this file')
        return
    
    data1 = spikes_per_T(spikes_dict, wells, electrodes, Ts, sample_time)
    data2 = amps_per_T(amps_dict, spikes_dict, wells, electrodes, Ts, sample_time)
    
    #plot_time_windows(wells, Ts, plotting_Ts, data1, output_dir)
    
    print("...done analyzing\n\n")
    
    return data1, data2, wells_color, plate_type


def main():

    # arg parse    
    description = 'Process Axion data files.'
    epilog = 'Good luck Shira!'
    usage = 'script.py --type type --input-dir ipath --output-dir opath --time-gate gate --all_dir\n\
    \t  type: 1 - axion, 2 - multichannel systems\n\
    \t  input-dir: path for spikes list directory\n\
    \t  output-dir: path for output directory\n\
    \t  output-file: name of output file\n\
    \t  gate: time gate in seconds\n\
    \t  eg: script.py --type axion --input-dir=\"..\\spike_list.csv\" --output-dir=\"..\\output.csv\" --time-gate 60\n\
    \t      script.py -t axion -i \"..\\spike_list.csv\" -o \"..\\output.csv\" -T 60\n'
    parser = argparse.ArgumentParser(description=description,
                                     epilog=epilog, usage=usage)
    parser.add_argument("-t", "--type", help="type of input file. 1 - axion, 2 - multichannel systems", default='axion')
    parser.add_argument("-i", "--input-dir", help="path for spikes list directory", required=True)
    parser.add_argument("-o", "--output-dir", help="path for output directory", required=True)
    parser.add_argument("-f", "--output-file", help="output file name", required=True)
    parser.add_argument("-T", "--time-gate", help="time gate in seconds", required=True, nargs='*')
    parser.add_argument("-d", "--read-all-dir", action='store_true', help='process all csv files in input-dir automatically. No GUI.')
    args = parser.parse_args()
    
    input_type = args.type    
    input_dir = args.input_dir    
    output_dir = args.output_dir
    output_file = args.output_file
    Ts = args.time_gate
    
    # check if a specific file is given, not a directory
    if re.match(r'.*\.csv$',input_dir): 
        input_files = [input_dir]
    else:        
        use_gui = False if args.read_all_dir else True
    
        # read list of csv files for input
        if use_gui:        
            root = tk.Tk()
            root.filenames = filedialog.askopenfilenames(initialdir = input_dir, title = "Select file", filetypes = (("csv files","*.csv"),("all files","*.*")))    
            input_files = root.tk.splitlist(root.filenames)
        else:        
            input_files = [f for f in glob.glob(input_dir + "\\*.csv")]       
    
    # output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    Ts = [float(T.strip('*')) for T in Ts]
    Ts.sort()
    #plotting_Ts = [float(T.strip('*')) for T in sys.argv[3:-1] if '*' in T]

    # check extension of output file name
    root, ext = os.path.splitext(output_file)
    if not ext:
        ext = '.xlsx'
        output_file = root + ext
    ofile = os.path.join(output_dir,output_file)
    
    # create a workbook
    with xlsxwriter.Workbook(ofile,{'constant_memory': True}) as workbook:

        for input_file in sorted(input_files):
            data1, data2, wells_color, plate_type = analyze_file(input_type, Ts, input_file)
        
            # save raw spikes and amps count per timebin            
            ofname = ofile.replace(ext, '_SPK' + ext)
            save_raw(data1, Ts, ofname, workbook, wells_color, plate_type)
            ofname = ofile.replace(ext, '_AMP' + ext)
            save_raw(data2, Ts, ofname, workbook, wells_color, plate_type)

        


if __name__=='__main__':
    
    main()
    