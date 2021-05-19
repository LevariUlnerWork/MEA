'''
This file is the main file for the MEA project, it will active the specific script per excel
'''



import numpy as np
import sys
import csv
import re
import matplotlib

matplotlib.use('agg')
import os
import importlib
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
import time
import xlsxwriter



def main():
    os.system("cls")
    print ("      ____          ")
    print ("|     |     \      /")
    print ("|     ----   \    /")
    print ("|     |       \  /")
    print ("----  ----     \/")

    print("Welcome to MEA Analyze code")
    time.sleep(1)
    print("Please choose the csv files you want to analyze ")
    time.sleep(1)
    analyze_burst = __import__("analyze_burst")
    analyze_brain = __import__("analyze_brain")
    analyze_burst_network = __import__("analyze_burst_network")
    input_dir = "C:\\Users\\קוגניציה מולקולרית\\Desktop"
    root = tk.Tk()
    root.filenames = filedialog.askopenfilenames(initialdir=input_dir, title="Select file",
                                                 filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
    input_files = root.tk.splitlist(root.filenames)
    root.destroy()
    burst_files = []
    burst_network_files = []
    spikes_files = []
    error_files = []

    Ts = float(input("Please enter the Ts" + '\n'))

    print("Please choose the output folder:")
    time.sleep(2)

    root = tk.Tk()
    root.filenames = filedialog.askdirectory(initialdir=input_dir, title="Select folder")
    output_dir = root.filenames
    root.destroy()
    for file in input_files:
        filenameList = file.split('/')
        filename = filenameList[len(filenameList) - 1].split(".")[0]
        if len(filename) > 32:
            if (re.search("electrode_burst_list" , file)!= None) or (re.search("network_burst_list", file )!= None) or (re.search("network_burst_list", file )!= None):
                error_files.append(filename)
        elif (re.search("electrode_burst_list" , file)!= None):
            burst_files.append(file)

        elif (re.search("network_burst_list", file )!= None):
            burst_network_files.append(file)

        elif (re.search("spike_list", file )!= None):
            spikes_files.append(file)
    if(len(burst_files) > 0):
        analyze_burst.main(input_dir, burst_files, Ts,output_dir)
    if (len(burst_network_files) > 0):
        analyze_burst_network.main(input_dir,  burst_network_files, Ts, output_dir)
    if(len(spikes_files) > 0):
        analyze_brain.main(input_dir, spikes_files, Ts,output_dir)

    if (len(error_files) > 0):
        print ("The following files have more than 32 characters and didn't analyzed.")
        print("Please change their names:")
        for file in error_files:
            print(file)


    print ("Press any key to exit")

if __name__ == '__main__':
    main()