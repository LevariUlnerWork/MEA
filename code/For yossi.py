import csv
import numpy as np
import tkinter as tk
from tkinter import filedialog

# root = tk.Tk()
# root.filenames = filedialog.askopenfilenames(initialdir='C:/', title="Select file",
#                                              filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
# input_file = root.tk.splitlist(root.filenames)
# root.destroy()
path = 'C:/Users/levco/Downloads/V+D1(000)spike_list.csv'
with open(path, 'r', encoding='UTF-8') as burstsfile:
    # load spikes
    print("...reading file")
    lines = burstsfile.readlines();
    print("...done")
data_to_save = []
data_to_save.append(["Investigator","Yossi","Time (s)","Electrode","Amplitude (mV)"])
for line in range(1,37):
    if(lines[line]==data_to_save[0]):continue
    line_as_list = lines[line].split(',') [:2]
    data_to_save.append(line_as_list)
for line in range(-8,0):
    data_to_save.append(lines[line].split(','))

c = 0
ext = '.xlsx'
ofnameBRST = path.replace(ext, str(c) + ext)
writeFirstLine=True
filename = f'C:/Users/levco/Desktop/מחקר שירה/Yossi/spike_list_yossi{c}.csv'

for line in range (len(lines)):
    with open (filename, 'a' ,newline='', encoding='utf-8') as csvFile:
        w = csv.writer(csvFile, delimiter=',')
        if (writeFirstLine == True):
            data_to_write = data_to_save[0]
            w.writerow(data_to_write)
            writeFirstLine = False
        elif(line % 400000 == 0):
            w.writerow('')
            for i in range (8):
                data_to_save[37+i][len(data_to_save[37+i])-1] = data_to_save[37+i][len(data_to_save[37+i])-1].replace('\n','')
                w.writerow(data_to_save[37+i])
            c += 1
            filename = f'C:/Users/levco/Desktop/מחקר שירה/Yossi/spike_list_yossi{c}.csv'
            print("Changing file")
            writeFirstLine = True

        elif(line%400000 < 37):
            data_to_write = []
            for i in range (5):
                if i<2:
                    data_to_write.append(data_to_save[line%400000][i].replace('\n','')+",")
                else:
                    data_to_write.append(lines[line].split(',')[i].replace('\n','') + ",")
            data_to_write.append('\n')
            csvFile.writelines(data_to_write)
        else:
            data_to_write = []
            for i in range (5):
                data_to_write.append(lines[line].split(',')[i].replace('\n','') + ",")
            data_to_write.append('\n')
            csvFile.writelines(data_to_write)
        csvFile.close()

with open (filename, 'a' ,newline='', encoding='utf-8') as csvFile:
    w = csv.writer(csvFile, delimiter=',')
    w.writerow('')
    for i in range (8):
        data_to_save[37+i][len(data_to_save[37+i])-1] = data_to_save[37+i][len(data_to_save[37+i])-1].replace('\n','')
        w.writerow(data_to_save[37+i])
    csvFile.close()



