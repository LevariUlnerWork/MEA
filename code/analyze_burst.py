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
__author__ = "Or Nahmani" + ", Lev-Ari Ulner"
__email__ = "ornahm@post.bgu.ac.il" + ", ulnerl@post.bgu.ac.il"
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
DURATION_IND = 5
SIZE_SPIKES_IND = 4
ELECTRODE_IND = 3
TIMESTAMP_IND = 2
N_PLATE_TYPE_24 = 4  # 16 elecs
N_PLATE_TYPE_6 = 8  # 64 elecs
MAX_SHEET_NAME_LENGTH = 31
SPK_RATE_TOL = 1  # spikes/sec, tollerance for a "living" electrod

REC_TIME_RE = '.* to ((?P<hr>\\d+)h)?((?P<min>\\d+)m)?(?P<sec>\\d+)s'
FOOTER_RE = '^Well Information.*'
LAST_DATA_LINE_RE = '.*,(?P<elec>[A-Z]\d_\d{2}),.*'
REC_TIME_LINE_RE = '.*Actual File Section Run*'
PLATE_TYPE_LINE_RE = '\s+Plate Type'
TREAT_LINE_RE = 'Treatment'
CONC_LINE_RE = 'Concentration'
NUM_WELLS_RE = '.*MEA (?P<ptype>\\d+)'


def parse_axion_bursts_list(path):
    # Gets path of axion's _bursts_list.csv file and returns a dict of WellxElectrode: spikes ts array
    try:
        with open(path, 'r', encoding='UTF-8') as burstsfile:
            # load spikes
            print("...reading file")
            lines = burstsfile.readlines();
            print("...done")

        # start parsing file
        print("...parsing file")

        # find header - Well Information
        iline = -1
        while re.match(FOOTER_RE, lines[iline]) == None:
            iline -= 1
            if abs(iline) == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND FOOTER (STRING: Well Information)\n')
        data_wells = lines[iline + 1].rstrip().split(",")[1:]
        data_colors = lines[iline + 3].rstrip().replace("#", "").split(",")[1:]

        # find Treatment line
        iline = -1
        while re.match(TREAT_LINE_RE, lines[iline]) == None:
            iline -= 1
            if abs(iline) == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND TREATMENT (STRING: Treatment)\n')
        data_treat = lines[iline].rstrip().split(",")[1:]
        if len(data_treat) != len(data_wells):
            sys.exit('\n>>>>> ERROR: LENGTH OF TREATMENT != WELLS\n')
            # define dict
        treatment = defaultdict(list)
        for i, j in zip(data_wells, data_treat):
            treatment[i].append(j)

        # find Concentration line
        iline = -1
        while re.match(CONC_LINE_RE, lines[iline]) == None:
            iline -= 1
            if abs(iline) == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND CONCENTRATION (STRING: Concentration)\n')
        data_conc = lines[iline].rstrip().split(",")[1:]
        if len(data_conc) != len(data_wells):
            sys.exit('\n>>>>> ERROR: LENGTH OF CONCENTRATION != WELLS\n')
            # define dict
        concentration = defaultdict(list)
        for i, j in zip(data_wells, data_conc):
            concentration[i].append(j)

        # find last data line
        iline = -1
        while re.match(LAST_DATA_LINE_RE, lines[iline]) == None and abs(iline) < len(lines):
            iline -= 1
        data = np.loadtxt(lines[1:iline + 1], dtype="str", delimiter=",")

        # parse the recording's length
        iline = 0
        while re.match(REC_TIME_LINE_RE, lines[iline]) == None:
            iline += 1
            if iline == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND RECORDING TIME (REQUESTED/ACTUAL DURATION OF MEASUREMENT)\n')
        time_match = re.match(REC_TIME_RE, lines[iline]).groupdict()

        rec_time = 0
        if time_match['hr']:
            rec_time += int(time_match['hr']) * 3600
        if time_match['min']:
            rec_time += int(time_match['min']) * 60
        if time_match['sec']:
            rec_time += int(time_match['sec'])

        # parse the plate type - CytoView MEA 6/24
        iline = 0
        while re.match(PLATE_TYPE_LINE_RE, lines[iline]) == None:
            iline += 1
            if iline == len(lines):
                sys.exit('\n>>>>> ERROR: COULD NOT FIND PLATE TYPE (CytoView MEA 6/24)\n')
        plate_type = int(re.match(NUM_WELLS_RE, lines[iline]).group('ptype'))

        # Parse footer (colors for sort)
        wells_color = defaultdict(list)
        for i, j in zip(data_colors, data_wells):
            wells_color[i].append(j)
        wells_color = OrderedDict(sorted(wells_color.items()))  # sort dict

        for i in wells_color.keys():  # sort list of values per key
            wells_color[i].sort(key=lambda s: s[::-1])

        electrodes = set(data[:, ELECTRODE_IND])

        # Assuming all wells are the same (containing same electrodes alignment)
        elecs_per_well = {elec[-2:] for elec in electrodes}
        wells = {elec[:2] for elec in electrodes}

        # Initialize spikes dict for all electrodes
        elecs = {elec: [] for elec in elecs_per_well}
        # import pdb; pdb.set_trace()
        timestamps_lists = {well: copy.deepcopy(elecs) for well in wells}
        bursts_lists = {well: copy.deepcopy(elecs) for well in wells}
        durations_lists = {well: copy.deepcopy(elecs) for well in wells}
        for row in data:
            well, elec = row[ELECTRODE_IND].split('_')
            ts = float(row[TIMESTAMP_IND])
            size = float(row[SIZE_SPIKES_IND])
            duration = float(row[DURATION_IND])
            timestamps_lists[well][elec].append(ts)
            bursts_lists[well][elec].append(size)
            durations_lists[well][elec].append(duration)

        # Convert to numpy arrays for further proccessing. this is not done initially because np.append copies the array.
        for well in wells:
            timestamps_lists[well] = {elec: np.array(br_list) for elec, br_list in timestamps_lists[well].items() if
                                      br_list != []}
            bursts_lists[well] = {elec: np.array(br_list) for elec, br_list in bursts_lists[well].items() if
                                  br_list != []}
            durations_lists[well] = {elec: np.array(dur_list) for elec, dur_list in durations_lists[well].items() if
                                     dur_list != []}

        print("...done")

    except IOError:
        print("File not accessible")

    # pack returns
    brst_data = timestamps_lists, bursts_lists, durations_lists, wells, elecs_per_well, \
                rec_time, wells_color, plate_type, \
                treatment, concentration

    return brst_data


# Not use:
def parse_mcs_spikes_list(path):
    # Gets path of multi channel system's spike list (exported from Multi Channel Analyzer) and returns a dict of Electrode: spikes ts array
    try:
        with open(path) as spikesfile:
            data = np.loadtxt(spikesfile.readlines()[3:], dtype="str", delimiter="\t")

        wells = ['']  # A single well
        # Set the sample time to be the last spike's timestamp
        rec_time = 0
        spikes_lists = {}
        for col in data.T:
            elec = col[0][-2:]
            # timestamps are in µs. convert to seconds
            ts_list = [float(ts) / (10 ** 6) for ts in col[1:] if ts != '']
            spikes_lists[elec] = np.array(ts_list)
            rec_time = max(rec_time, max(ts_list))

    except IOError:
        print("File not accessible")

    return {'': spikes_lists}, wells, rec_time


# gets parsed spikeslist dict,returns a dict structure of TxWELLxELECTRODExSPT and plots the histograms
def bursts_per_T(elec_timestamps_dict, elec_bursts_dict, wells, electrodes, Ts, sample_time):
    #
    # Ts - time durations (s)
    # sample_time - total sample time (s)
    print("...calculating bursts per T")
    data = {}
    # configure histograms figure
    # ax_size = math.ceil(math.sqrt((len(electrodes)-1))) # num of rows/cols in the MEA grid

    data = {}
    bins = np.arange(0, sample_time + 1, Ts)
    for well in wells:
        well_data = {}
        for elec in elec_bursts_dict[well].keys():
            timestamp = elec_timestamps_dict[well][elec]  # time stamps of spikes
            bursts = elec_bursts_dict[well][elec]
            well_data[elec] = np.nan_to_num(stats.binned_statistic(timestamp, bursts, statistic='mean', bins=bins)[0])
            # well_data = {elec: np.histogram(spikes, bins)[0] for elec, spikes in elec_amps_dict[well].items()}

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

            # Add mean well's durations as well
            tot_all_bur = np.zeros(len(bins) - 1)
            for j in range(len(bins) - 1):
                tot_bur = 0.0
                icount = 0
                for elec in well_data.keys():
                    if well_data[elec][j] > 0:
                        tot_bur += well_data[elec][j]
                        icount += 1
                if icount > 0:
                    tot_all_bur[j] = (tot_bur)  # to do

            well_data[""] = tot_all_bur

            data[well] = well_data
    print("...done")

    return data


# gets parsed spikeslist dict,returns a dict structure of TxWELLxELECTRODExSPT and plots the histograms
def durations_per_T(elec_durations_dict, elec_timestamps_dict, wells, electrodes, Ts, sample_time):
    #
    # Ts - time durations (s)
    # sample_time - total sample time (s)
    print("...calculating durations per Ts")
    data = {}
    # configure histograms figure
    # ax_size = math.ceil(math.sqrt((len(electrodes)-1))) # num of rows/cols in the MEA grid

    data = {}
    bins = np.arange(0, sample_time + 1, Ts)
    for well in wells:
        well_data = {}
        for elec in elec_durations_dict[well].keys():
            timestamp = elec_timestamps_dict[well][elec]  # time stamps of  bursts
            durations = elec_durations_dict[well][elec]
            well_data[elec] = np.nan_to_num(
                stats.binned_statistic(timestamp, durations, statistic='mean', bins=bins)[0])
            # well_data = {elec: np.histogram(spikes, bins)[0] for elec, spikes in elec_amps_dict[well].items()}

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

            # Add mean well's durations as well
            tot_all_dur = np.zeros(len(bins) - 1)
            for j in range(len(bins) - 1):
                tot_dur = 0.0
                icount = 0
                for elec in well_data.keys():
                    if well_data[elec][j] > 0:
                        tot_dur += well_data[elec][j]
                        icount += 1
                if icount > 0:
                    tot_all_dur[j] = tot_dur

            well_data[""] = tot_all_dur

            data[well] = well_data

    print('...done')

    return data


def plot_heatmaps(data, well, T, labels, path):
    # data = 2D np array of single well's electrodes bursts per T
    fig, axs = plt.subplots(1, 2, figsize=(15, 10))
    fig.suptitle("Well %s, T = %.3f" % (well, T))
    cov = np.cov(data)
    cov_g = sns.heatmap(cov, center=0, cmap=sns.diverging_palette(20, 220, n=200), ax=axs[0], square=True)
    cov_g.set_xticklabels(labels, rotation=45, horizontalalignment='right');
    cov_g.set_yticklabels(labels, rotation=0);
    cov_g.set_title('Electrodes Covariance')

    corrcoef = np.corrcoef(data)
    corrcoef_g = sns.heatmap(corrcoef, vmin=-1, vmax=1, center=0, cmap=sns.diverging_palette(20, 220, n=200), ax=axs[1],
                             square=True)
    corrcoef_g.set_xticklabels(labels, rotation=45, horizontalalignment='right');
    corrcoef_g.set_yticklabels(labels, rotation=0);
    corrcoef_g.set_title('Electrodes Corrcoef')

    # plt.savefig(path)
    plt.close(fig)


def save_mean_vars(yfuncs_output_csv, Ts, output_path):
    output_title = ["Well/Electrode"] + ','.join(
        ["T=%.3f mean,T=%.3f variance,T=%.3f var/mean" % (T, T, T) for T in Ts]).split(',')
    with open(output_path, 'w', newline='') as outfile:
        w = csv.writer(outfile, delimiter=',')
        w.writerow(output_title)
        # write ordered: wells first then electrodes.
        for item in sorted(sorted(yfuncs_output_csv.keys()), key=len):
            w.writerow([item] + yfuncs_output_csv[item])


def get_elecs_full_list(well, plate_type):
    n = N_PLATE_TYPE_24 if plate_type == 24 else N_PLATE_TYPE_6

    elec_list = []
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            elec_list += [str(10 * i + j)]

    return elec_list


def write_line(ws, wb, row, line, color):
    for col in range(len(line)):
        if col == 0:  # well
            cell_format = wb.add_format({'bold': True, 'bg_color': color})
            ws.write(row, col, line[col], cell_format)
        else:
            cell_format = wb.add_format({'bold': False, 'bg_color': color})
            ws.write(row, col, line[col], cell_format)


def write_sheet(data, data2, Ts, wb, ws, wells_color, plate_type, treatment, concentration):
    print("...writing sheet: %s" % ws.get_name())

    row = 0

    # span time axis
    exist_well = [*data.keys()][0]
    exist_elec = [*data[exist_well].keys()][0]
    nbins = len(data[exist_well][exist_elec])
    xtime = np.linspace(Ts, nbins * Ts, nbins)
    line = ['Time'] + ['Treatment (Concentration)'] + xtime.tolist()

    # write
    bold = wb.add_format({'bold': True})
    for col in range(len(line)):
        ws.write(row, col, line[col], bold)
    row += 1

    # write ordered: wells first then wells+electrodes.
    for color in wells_color:
        for well in wells_color[color]:
            if concentration[well] == ['']:
                tc = "".join(treatment[well])
            else:
                tc = "".join(treatment[well] + [' ('] + concentration[well] + [')'])
            if well in data:
                electrodes = data[well]

                # print only wells first
                for elec in sorted(electrodes.keys()):
                    if (elec == ''):
                        ln = electrodes[elec].tolist()
                        ls = [x if x > 0 else "" for x in ln]
                        line = [well] + [tc] + ls
                        # write
                        write_line(ws, wb, row, line, color)
                        row += 1
            else:
                ls = [''] * nbins
                line = [well] + [tc] + ls
                write_line(ws, wb, row, line, color)
                row += 1

    row += 1  # empty line

    row += 1
    ws.write(row, 0, 'Durations(s):', bold)
    row += 1

    # write ordered: durations wells.
    for color in wells_color:
        for well in wells_color[color]:
            if well in data2:
                electrodes = data2[well]

                # print only wells first
                for elec in sorted(electrodes.keys()):
                    if (elec == ''):
                        ln = electrodes[elec].tolist()
                        ls = [x if x > 0 else "" for x in ln]
                        line = [well] + [tc] + ls
                        # write
                        write_line(ws, wb, row, line, color)
                        row += 1
            else:
                ls = [''] * nbins
                line = [well] + [tc] + ls
                write_line(ws, wb, row, line, color)
                row += 1

    row += 1  # empty line

    # write ordered: wells+electrodes
    for color in wells_color:
        for well in wells_color[color]:
            if well in data:
                electrodes = data[well]

                # print all
                elecs_full_list = get_elecs_full_list(well, plate_type)
                for elec in elecs_full_list:
                    ttl = well if (elec == '') else (well + '_' + elec)
                    if elec in electrodes.keys():
                        line = electrodes[elec].tolist()
                        ls = [x if x > 0 else "" for x in line]
                    else:
                        ls = [''] * nbins
                    line = [ttl] + [''] + ls
                    write_line(ws, wb, row, line, color)
                    row += 1
                    # w.writerow([ttl] + ls)
            else:
                ls = [''] * nbins
                line = [well] + [''] + ls
                write_line(ws, wb, row, line, color)
                row += 1

                elecs_full_list = get_elecs_full_list(well, plate_type)
                for elec in elecs_full_list:
                    ttl = well if (elec == '') else (well + '_' + elec)
                    ls = [''] * (nbins + 1)
                    line = [ttl] + ls
                    write_line(ws, wb, row, line, color)
                    row += 1

                    print("...done")


    print("...done")


# Didnt use:
# gets window size, and plots yfuncs and cov heatmaps for the rolling time window
def plot_time_windows(wells, Ts, plotting_Ts, data, output_dir, window_size=0):
    # Create yfuncs_output_csv - a dict of elec/well: list of mean,variance,mean,variance... for each given T in list, and in that order.
    yfuncs_output_csv = dict()
    wells_yfuncs = {well: [] for well in wells}
    for T in Ts:
        for well in wells:
            well_T_data = data[T][well]
            # update mean&var data - for the csv output
            for elec, spt in well_T_data.items():
                mean, var = spt.mean(), spt.var()
                print("%s  %s  %f  %f" % (well, elec, mean, var))
                key = well + "_" + elec if elec else well
                if key not in yfuncs_output_csv.keys():
                    yfuncs_output_csv[key] = []
                yfuncs_output_csv[key].extend([mean, var, var / mean])

            # Add well's yfunc (=var/mean) - for plotting
            wells_yfuncs[well].append(well_T_data[""].var() / well_T_data[""].mean())

            # plot covariance and corrcoef matrices
            if T in plotting_Ts and len(
                    well_T_data.keys()) > 2:  # plot if we have at least 2 active electrodes. >2 because the well ("") is also a key
                electrodes = sorted(well_T_data.keys())
                electrodes.remove("")
                # Stack active electrodes in order (electrodes is a sorted list)
                spt_stack = np.stack([well_T_data[elec] for elec in electrodes], axis=0)
                plot_heatmaps(spt_stack, well, T, electrodes, os.path.join(output_dir, "covs_%s_%.3f.png" % (well, T)))

    # Save yfuncs_output_csv to file
    # save_mean_vars(yfuncs_output_csv, Ts, os.path.join(output_dir,"fano_factors.csv"))
    # Plot yfuncs
    plt.figure(figsize=(10, 10))
    for well in sorted(wells_yfuncs.keys()):
        plt.plot(Ts, wells_yfuncs[well], label=well, marker='.')
    plt.legend()
    # plt.xlim(0, Ts[-1]+0.5)
    plt.xticks(Ts)
    plt.xlabel('T')
    plt.ylabel('var(T)/mean(T)')
    # plt.savefig(os.path.join(output_dir,"fano_factors.png"))


def analyze_file(input_file, Ts, Tw, Tc):
    print("...analyzing file: %s" % input_file)

    parsing_func = parse_axion_bursts_list
    brst_data = parsing_func(input_file)

    timestamps_dict, bursts_dict, durations_dict, wells, electrodes, sample_time, \
    wells_color, plate_type, treatment, concentration = brst_data

    # check time gate
    print('...sample time is: %d s' % sample_time)
    if Ts > sample_time:
        # sys.exit('\n>>>>> ERROR: TIME GATE (T) IS LARGER THAN THE SAMPLE LENGTH\n')
        print('...WARNING: time gate (T=%d s) is larger than the sample length (%d s)' % (Ts, sample_time))
        print('...skipping analysis of this file')
        return

    data1 = bursts_per_T(timestamps_dict, bursts_dict, wells, electrodes, Ts, sample_time)
    data2 = durations_per_T(durations_dict, timestamps_dict, wells, electrodes, Ts, sample_time)

    # plot_time_windows(wells, Ts, plotting_Ts, data1, output_dir)

    print("...done analyzing\n\n")

    return data1, data2, wells_color, plate_type, treatment, concentration


def main():
    # arg parse
    description = 'Process Axion data files.'
    epilog = 'Good luck Shira!'

    usage = 'script.py --input-dir=ipath --output-dir=opath --output-file=fname --time-gate-spikes=Ts  --time-gate-correlation=Tw --time-gate-correlation-total=Tc --read-all_dir\n\
        \t  input-dir: path for spikes list directory\n\
        \t  output-dir: path for output directory\n\
        \t  output-file: name of output file\n\
        \t  time-gate-spikes: time gate for spikes [s]\n\
        \t  time-gate-correlation: time gate for correlations [s]\n\
        \t  time-gate-correlation-total: time gate for average correlations [s]\n\
        \t  read-all-dir: process all csv files in input-dir automatically. No GUI.\n\
        \t  e.g., from MEA/code, run:\n\
        \t  python3 analyze_brain.py -i "..\data\spike_list_1.csv" -o "..\output" -f "spk5" -Ts 120 -Tw 0.1 -Tc 10\n'
    parser = argparse.ArgumentParser(description=description,
                                     epilog=epilog, usage=usage)
    parser.add_argument("-i", "--input-dir", help="path for spikes list directory", required=True)
    parser.add_argument("-o", "--output-dir", help="path for output directory", required=True)
    parser.add_argument("-f", "--output-file", help="output file name", required=True)
    parser.add_argument("-Ts", "--time-gate-spikes", help="time gate in seconds for spikes", required=True, nargs='*', type=float, metavar='Ts', dest='Ts')
    parser.add_argument("-Tw", "--time-gate-correlation", help="time gate in seconds for spikes", required=True, nargs='*', type=float, metavar='Tw', dest='Tw')
    parser.add_argument("-Tc", "--time-gate-correlation-total", help="time gate in seconds for correlation", required=True, nargs='*', type=float, metavar='Tc', dest='Tc')
    parser.add_argument("-d", "--read-all-dir", action='store_true', help='process all csv files in input-dir automatically. No GUI.')
    args = parser.parse_args()

    # input_type = args.type
    input_dir = args.input_dir
    output_dir = args.output_dir
    output_file = args.output_file


    Ts = args.Ts[0]
    Tw = args.Tw[0]
    Tc = args.Tc[0]
    '''
    Ts = 120
    Tw = 0.1
    Tc = 10
    input_dir = "..\\data\\burst\\50k div12 cortex 20200804 After venom overnight(010)(000)_electrode_burst_list.csv"
    output_dir = "..\output"
    output_file = "brst55"
    '''
    # check if a specific file is given, not a directory
    if re.match(r'.*\.csv$', input_dir):
        input_files = [input_dir]
    else:
        use_gui = False if args.read_all_dir else True

        # read list of csv files for input
        if use_gui:
            root = tk.Tk()
            root.filenames = filedialog.askopenfilenames(initialdir=input_dir, title="Select file",
                                                         filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
            input_files = root.tk.splitlist(root.filenames)
        else:
            input_files = [f for f in glob.glob(input_dir + "\\*.csv")]

        # output dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    # check extension of output file name
    root, ext = os.path.splitext(output_file)
    if not ext:
        ext = '.xlsx'
        output_file = root + ext
    ofile = os.path.join(output_dir, output_file)

    # create workbooks - open xlsx files for writing
    print("...open output files for writing (%s)" % ofile)
    #ofnameBRST = ofile.replace(ext, '_BRST' + ext)
    ofnameBRST = ofile
    wbBRST = xlsxwriter.Workbook(ofnameBRST, {'constant_memory': True})
    print("...done")

    # loop over input files
    for input_file in sorted(input_files):
        data1, data2, wells_color, plate_type, treat, conc \
            = analyze_file(input_file, Ts, Tw, Tc)

        # save raw spikes and amps count per timebin
        ifile = os.path.basename(input_file)
        root, ext = os.path.splitext(ifile)
        # len(root) <= 31
        root = root.replace("_spike_list", "")
        root = root.replace("electrode_burst_list", "")
        root = root.replace("network_burst_list", "")
        if len(root) > MAX_SHEET_NAME_LENGTH:
            root = root[-MAX_SHEET_NAME_LENGTH:]  # leave last 31 chars
        # add s worksheets
        wsBRST = wbBRST.add_worksheet(root)

        write_sheet(data1, data2, Ts, wbBRST, wsBRST, wells_color, plate_type, treat, conc)
        # save_raw(data1, data2, Ts, ofnameBRST, workbook, wells_color, plate_type)
        # ofDur = ofile.replace(ext, '_AMP' + ext)
        # save_raw(data2, Ts, ofDur, workbook, wells_color, plate_type)

        print("...closing output files")
        wbBRST.close()
        print("...DONE\n\n")


if __name__ == '__main__':
    main()
