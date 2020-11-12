#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is a part of the analyze_brain file and is written to handle and find bursts in the sample.

ToDo:
* Output raw histograms
* Sort histograms by distance between elecs
* Remove multiple Ts from data
* Find bursts and save them
*


Things to note:
* Output file is showing Ts with 3 decimal points. Should it be more?
* For Axion's spike_list.csv analysis: the data file contains a 9 lines footer, which this script ignores. Update that if anything changes.
"""

# part_Owned
__author__ = "Or Nahmani, LevAri Ulner"
__email__ = "ulnerl@post.bgu.ac.il"
__company__ = "NIBN, Ben-Gurion University of the Negev"
__date__ = "05/11/2020"
__license__ = "MIT License"
__copyright__ = "Copyright 2020, NIBN, Israel"
__status__ = "Dev"
__credits__ = [__author__]
__version__ = "0.3.0"

import numpy as np
import sys
import re
import matplotlib

matplotlib.use('agg')
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
DURATION_IND = 5 # Lev change = The column of the duration of the burst
SIZE_SPIKES_IND = 4 # Lev change = The column of the number of spikes worked in that burst
ELECTRODE_IND = 3
TIMESTAMP_IND = 2
N_PLATE_TYPE_24 = 4  # 16 elecs
N_PLATE_TYPE_6 = 8  # 64 elecs
MAX_SHEET_NAME_LENGTH = 31
SPK_RATE_TOL = 1  # spikes/sec, tollerance for a "living" electrode

REC_TIME_RE = '.* to ((?P<hr>\\d+)h)?((?P<min>\\d+)m)?(?P<sec>\\d+)s'
FOOTER_RE = '^Well Information.*'
LAST_DATA_LINE_RE = '.*,(?P<elec>[A-Z]\d_\d{2}),.*'
REC_TIME_LINE_RE = '.*Actual File Section Run*'
PLATE_TYPE_LINE_RE = '\s+Plate Type'
TREAT_LINE_RE = 'Treatment'
CONC_LINE_RE = 'Concentration'
NUM_WELLS_RE = '.*MEA (?P<ptype>\\d+)'


def parse_axion_spikes_list(path):
    # Gets path of axion's _spike_list.csv file and returns a dict of WellxElectrode: spikes ts array
    try:
        with open(path, 'r', encoding='UTF-8') as spikesfile:
            # load spikes
            print("...reading file")
            lines = spikesfile.readlines();
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

        # Initialize bursts dict for all electrodes
        elecs = {elec: [] for elec in elecs_per_well}
        # import pdb; pdb.set_trace()
        spikes_lists = {well: copy.deepcopy(elecs) for well in wells}#Would change to spikes size
        duration_lists = {well: copy.deepcopy(elecs) for well in wells}
        for row in data:
            well, elec = row[ELECTRODE_IND].split('_')
            ts = float(row[TIMESTAMP_IND])
            dura = float(row[DURATION_IND])
            spikes_lists[well][elec].append(ts)
            duration_lists[well][elec].append(dura)

        # Convert to numpy arrays for further proccessing. this is not done initially because np.append copies the array.
        for well in wells:
            spikes_lists[well] = {elec: np.array(sp_list) for elec, sp_list in spikes_lists[well].items() if
                                  sp_list != []}
            duration_lists[well] = {elec: np.array(amp_list) for elec, amp_list in duration_lists[well].items() if
                                amp_list != []}

        print("...done")

    except IOError:
        print("File not accessible")

    # pack returns
    spk_data = spikes_lists, wells, elecs_per_well, \
               rec_time, duration_lists, wells_color, plate_type, \
               treatment, concentration

    return spk_data


# gets parsed burstslist dict,returns a dict structure of TxWELLxELECTRODExSPT and plots the histograms
def spikes_per_T(elec_burst_dict, wells, electrodes, T, sample_time):
    #
    # T - time duration (s)
    # sample_time - total sample time (s)
    print("...calculating burst per T")
    data = {}
    # configure histograms figure
    # ax_size = math.ceil(math.sqrt((len(electrodes)-1))) # num of rows/cols in the MEA grid

    bins = np.arange(0, sample_time + 1, T)
    for well in wells:
        well_data = {elec: np.histogram(bursts, bins)[0] for elec, bursts in elec_burst_dict[well].items()}

        # Add total well's spt as well
        well_data[""] = sum(well_data.values())
        data[well] = well_data

    print("...done")

    return data

"""
Lev: Don't change this function for now
"""

# gets parsed bursts_list dict,returns a dict structure of TxWELLxELECTRODExSPT and plots the histograms
def amps_per_T(elec_size_dict, elec_bursts_dict, wells, electrodes, T, sample_time):
    #
    # Ts - time durations (s)
    # sample_time - total sample time (s)
    print("...calculating amplitudes per T")
    data = {}
    # configure histograms figure
    # ax_size = math.ceil(math.sqrt((len(electrodes)-1))) # num of rows/cols in the MEA grid

    bins = np.arange(0, sample_time + 1, T)
    for well in wells:
        well_data = {}
        for elec in elec_amps_dict[well].keys():
            spikes = elec_spikes_dict[well][elec]  # time stamps of spikes
            amps = elec_amps_dict[well][elec]
            well_data[elec] = np.nan_to_num(stats.binned_statistic(spikes, amps, statistic='mean', bins=bins)[0])

            # Add mean well's amp as well
        ave_amp = np.zeros(len(bins) - 1)
        for j in range(len(bins) - 1):
            tot_amp = 0.0
            icount = 0
            for elec in well_data.keys():
                if well_data[elec][j] > 0:
                    tot_amp += well_data[elec][j]
                    icount += 1
            if icount > 0:
                ave_amp[j] = tot_amp / icount

        well_data[""] = ave_amp

        data[well] = well_data

    print('...done')

    return data


def get_elecs_full_list(well, plate_type):
    n = N_PLATE_TYPE_24 if plate_type == 24 else N_PLATE_TYPE_6

    elec_list = ['']
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


def write_sheet(data, T, wb, ws, wells_color, plate_type, \
                treatment, concentration):
    print("...writing sheet: %s" % ws.get_name())

    row = 0

    # span time axis
    exist_well = [*data.keys()][0]
    exist_elec = [*data[exist_well].keys()][0]
    nbins = len(data[exist_well][exist_elec])
    xtime = np.linspace(T, nbins * T, nbins)
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

    row += 1
    # w.writerow('') # empty line

    # write wells+electrodes.
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


def calc_cov_corr(data, wells, Tw, Tc):
    # check time gate
    print('...Tw = %f, Tc = %f' % (Tw, Tc))
    if Tw > Tc:
        sys.exit('\n>>>>> ERROR: TIME GATE (Tw) IS LARGER THAN Tc\n')
    if ((Tc / Tw) % 1) != 0.0:
        sys.exit('\n>>>>> ERROR: TIME GATE (Tw) IS NOT AN INTEGRAL DIVISION OF Tc\n')
    nw = round(Tc / Tw)

    # list of dicts
    c = [dict() for x in range(nw)]
    # cov = [dict() for x in range(nw)]

    # allocate
    for n in range(nw):
        c[n] = {well: [] for well in wells}
        # cov = {well:[] for well in wells}

    # find dead electrodes
    iwell = next(iter(data))  # find a well
    ielec = next(iter(data[iwell]))  # find an electrode
    nbins = len(data[iwell][ielec])
    Ttot = nbins * Tw
    # nnz = np.count_nonzero(data[iwell][ielec])
    # ave_spike_rate = np.sum(spk_stack,axis=1)/(nbins*T)

    # calc cov per bin Tc
    for n in range(nw):
        # stack data
        for well in wells:
            electrodes = data[well]
            if len(electrodes.keys()) > 2:
                elec_tags = sorted(electrodes.keys())
                elec_tags.remove("")
                spk_stack = np.stack([electrodes[elec] for elec in elec_tags], axis=0)

                # cov[well] = np.cov(spk_stack)
                c[n][well] = np.corrcoef(spk_stack)

    return c


def analyze_file(input_file, Ts, Tw, Tc):
    # unpack data
    spikes_dict, wells, electrodes, sample_time, amps_dict, \
    wells_color, plate_type, treatment, concentration \
        = parse_axion_spikes_list(input_file)

    print("...analyzing SPIKES in file: %s" % input_file)

    # check time gate
    print('...sample time is: %d s' % sample_time)
    if Ts > sample_time:
        # sys.exit('\n>>>>> ERROR: TIME GATE (T) IS LARGER THAN THE SAMPLE LENGTH\n')
        print('...WARNING: time gate (T=%d s) is larger than the sample length (%d s)' % (Ts, sample_time))
        print('...skipping analysis of this file')
        return

    # analyze spikes
    data_spk = spikes_per_T(spikes_dict, wells, electrodes, Ts, sample_time)
    data_amps = amps_per_T(amps_dict, spikes_dict, wells, electrodes, Ts, sample_time)

    # analyze correlations
    # data_corr = spikes_per_T(spikes_dict, wells, electrodes, Tw, sample_time)
    # cov, corr = calc_cov_corr(data_corr, wells, Tw, Tc)

    print("...done analyzing")

    return data_spk, data_amps, wells_color, plate_type, treatment, concentration


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
    # parser.add_argument("-t", "--type", help="type of input file. 1 - axion, 2 - multichannel systems", default='axion')
    parser.add_argument("-i", "--input-dir", help="path for spikes list directory", required=True)
    parser.add_argument("-o", "--output-dir", help="path for output directory", required=True)
    parser.add_argument("-f", "--output-file", help="output file name", required=True)
    parser.add_argument("-Ts", "--time-gate-spikes", help="time gate in seconds for spikes", required=True, nargs='*',
                        type=float, metavar='Ts', dest='Ts')
    parser.add_argument("-Tw", "--time-gate-correlation", help="time gate in seconds for spikes", required=True,
                        nargs='*', type=float, metavar='Tw', dest='Tw')
    parser.add_argument("-Tc", "--time-gate-correlation-total", help="time gate in seconds for correlation",
                        required=True, nargs='*', type=float, metavar='Tc', dest='Tc')
    parser.add_argument("-d", "--read-all-dir", action='store_true',
                        help='process all csv files in input-dir automatically. No GUI.')
    args = parser.parse_args()

    # input_type = args.type
    input_dir = args.input_dir
    output_dir = args.output_dir
    output_file = args.output_file
    Ts = args.Ts[0]
    Tw = args.Tw[0]
    Tc = args.Tc[0]

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
    ofnameSPK = ofile.replace(ext, '_SPK' + ext)
    wbSPK = xlsxwriter.Workbook(ofnameSPK, {'constant_memory': True})
    ofnameAMP = ofile.replace(ext, '_AMP' + ext)
    wbAMP = xlsxwriter.Workbook(ofnameAMP, {'constant_memory': True})
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
        wsSPK = wbSPK.add_worksheet(root)
        wsAMP = wbAMP.add_worksheet(root)

        write_sheet(data1, Ts, wbSPK, wsSPK, wells_color, plate_type, treat, conc)
        write_sheet(data2, Ts, wbAMP, wsAMP, wells_color, plate_type, treat, conc)

    print("...closing output files")
    wbSPK.close()
    wbAMP.close()
    print("...DONE\n\n")


if __name__ == '__main__':
    main()
