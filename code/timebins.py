import numpy as np
import sys
import csv

# spikes list file fields: "Investigator,,Time (s),Electrode,Amplitude (mV)"
ELECTRODE_IND = 3
TIMESTAMP_IND = 2


def parse_spikes_list(path):
    # Gets path of axion's _spike_list.csv file and returns a dict of electrode: list of spike timestamps
    with open(path) as spikesfile:
        # Theres a footer 9 lines long. first line is the fields.  <- not very generic...
        data = np.loadtxt(spikesfile.readlines()[1:-9], dtype="str", delimiter=",")
        # Get last spike's timestamp to use as the sample time.
        last_spike = float(data[-1][TIMESTAMP_IND])
        electrodes = set(data[:, ELECTRODE_IND])
        # Initialize spikes dict
        electrode_spikes = {elec:[] for elec in electrodes}
        for row in data:
            electrode_spikes[row[ELECTRODE_IND]].append(float(row[TIMESTAMP_IND]))
        return electrode_spikes, last_spike


# divide spikes_timestamps into timebins of size T and return the mean&std of spike count per bin
def get_mean_std(spikes_timestamps, T, sample_time):
    # spikes_timestamps - array of timestamps (s) of spike accurance
    # T - time duration (s)
    # sample_time - total sample time (s)
    spikes_timestamps = np.array(spikes_timestamps)
    bins = np.arange(0, int(sample_time), T)
    # ? should I ignore the leftover spikes? (if T divides with a remainder)
    bins_sum, b = np.histogram(spikes_timestamps, bins)
    return bins_sum.mean(), bins_sum.std()



def main():
    if len(sys.argv)!=4:
        print("usage: script <spikes_list_path> <T> <output_file>")
        print("\t  spikes_list_path - axion's _spike_list.csv file path")
        print("\t  T - timebin length in seconds")
        print("\t  output_file - full path to output file")
        return
    spikes_dict, sample_time = parse_spikes_list(sys.argv[1])
    T = int(sys.argv[2])
    means_stds = {elec : get_mean_std(spikes, T, sample_time) for elec, spikes in spikes_dict.items()}
    with open(sys.argv[3], 'w', newline='') as outfile:
        w = csv.writer(outfile, delimiter=',')
        w.writerow(["Electrode", "Mean", "STD"])
        for elec, (mean, std) in means_stds.items():
            w.writerow([elec, mean, std])


if __name__=='__main__':
    main()