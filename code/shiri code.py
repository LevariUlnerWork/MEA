#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""

"""

#Owned:
__author__ = "Lev-Ari Ulner"
__email__ = "ulnerl@post.bgu.ac.il"
__company__ = "NIBN, Ben-Gurion University of the Negev"
__date__ = "17/02/2021"
__copyright__ = "Copyright 2021, NIBN, Israel"
__status__ = "Dev"
__credits__ = [__author__]
__version__ = "0.1.0"

#imports:

import os
import tkinter as tk
from tkinter import filedialog
import xlsxwriter



#Columns:
SIDE=1
COUNT = 7
STRENGTH = 8
rightData = {}
leftData = {}

def parse_data(path):
    # Gets path of axion's _spike_list.csv file and returns a dict of WellxElectrode: spikes ts array
    # try:
    with open(path, 'r', encoding='UTF-8') as mindfile:
        # load spikes
        print("...reading file")
        lines = mindfile.readlines();
        print("...done")


    # start parsing file
    print("...parsing file")

    for line in lines:

        line = line.split(',')

        if(line[0] == "Structure"):
            COUNT = line.index("Total")
            STRENGTH = COUNT+1

        #pass outside lines
        if(line[0] == "Outside" or line[0] == ""):
            continue

        # Taking data by side
        if(line[1] == "left"):
            if(line[0] not in leftData.keys()):
                leftData[line[0]] = {}
            leftData[line[0]]["count"] = line[COUNT]
            leftData[line[0]]["strength"] = line[STRENGTH]

        if (line[1] == "right"):
            if(line[0] not in rightData.keys()):
                rightData[line[0]] = {}
            rightData[line[0]]["count"] = line[COUNT]
            rightData[line[0]]["strength"] = line[STRENGTH]
# except:
#     print ("Apchi")


def writeSheet(wbMD,wsMD,leftData,rightData):
    mindParts = leftData.keys()
    lineCount = 0
    for part in mindParts:

        if(leftData[part]["count"] == '0' and rightData[part]["count"] == '0'):
            continue

        cell_format = wbMD.add_format({'bold': True})

        ans1 = str(float(leftData[part]["count"]) + float(rightData[part]["count"]) / 2)
        ans2 = str(float(leftData[part]["strength"]) + float(rightData[part]["strength"]) / 2)
        if not ((float(leftData["root"]["count"]) + float(rightData["root"]["count"])) == 0):
            ans3 = str((float(leftData[part]["count"]) + float(rightData[part]["count"])) / (float(leftData["root"]["count"]) + float(rightData["root"]["count"])))
        else:
            ans3 = "0"
        if not ((float(leftData["root"]["strength"]) + float(rightData["root"]["strength"])) == 0):
            ans4 = str((float(leftData[part]["strength"]) + float(rightData[part]["strength"])) / (float(leftData["root"]["strength"]) + float(rightData["root"]["strength"])))
        else:
            ans4 = "0"
        if not ((float(leftData[part]["strength"]) + float(rightData[part]["strength"])) == 0):
            ans5 = str(float(leftData[part]["strength"]) / (float(leftData[part]["strength"]) + float(rightData[part]["strength"])))
        else:
            ans5 = "0"
        if not ((float(leftData[part]["strength"]) + float(rightData[part]["strength"])) == 0):
            ans6 = str(float(rightData[part]["strength"]) / (float(leftData[part]["strength"]) + float(rightData[part]["strength"])))
        else:
            ans6 = "0"
        if not((float(leftData[part]["count"]) + float(rightData[part]["count"])) == 0):
            ans7 = str(float(leftData[part]["count"]) / (float(leftData[part]["count"]) + float(rightData[part]["count"])))
        else:
            ans7 = "0"
        if not ((float(leftData[part]["count"]) + float(rightData[part]["count"])) == 0):
            ans8 = str(float(rightData[part]["count"]) / (float(leftData[part]["count"]) + float(rightData[part]["count"])))
        else:
            ans8 = "0"
        line = ["count Avg (Left + Right/2)", ans1]
        line2 = ["strength Avg (Left + Right/2)", ans2]
        line3 = ["count prec (Left+ Right / all)", ans3]
        line4 = ["strength prec (Left+ Right / all)", ans4]
        line5 = ["left side strength prec (Left / all)", ans5]
        line6 = ["right side strength prec (Right / all)", ans6]
        line7 = ["left side count prec (Left / all)", ans7]
        line8 = ["right side count prec(Right / all)", ans8]
        wsMD._write(lineCount, 0, part + ":", cell_format)
        lineCount += 1
        wsMD._write(lineCount, 0, line[0], cell_format)
        wsMD._write(lineCount, 1, line[1], cell_format)
        lineCount+=1
        wsMD._write(lineCount, 0, line2[0], cell_format)
        wsMD._write(lineCount, 1, line2[1], cell_format)
        lineCount += 1
        wsMD._write(lineCount, 0, line3[0], cell_format)
        wsMD._write(lineCount, 1, line3[1], cell_format)
        lineCount += 1
        wsMD._write(lineCount, 0, line4[0], cell_format)
        wsMD._write(lineCount, 1, line4[1], cell_format)
        lineCount += 1
        wsMD._write(lineCount, 0, line5[0], cell_format)
        wsMD._write(lineCount, 1, line5[1], cell_format)
        lineCount += 1
        wsMD._write(lineCount, 0, line6[0], cell_format)
        wsMD._write(lineCount, 1, line6[1], cell_format)
        lineCount += 1
        wsMD._write(lineCount, 0, line7[0], cell_format)
        wsMD._write(lineCount, 1, line7[1], cell_format)
        lineCount += 1
        wsMD._write(lineCount, 0, line8[0], cell_format)
        wsMD._write(lineCount, 1, line8[1], cell_format)
        lineCount += 1
        lineCount += 1




#choose file
print("Please choose the csv files you want to analyze ")
root = tk.Tk()
root.filename = filedialog.askopenfilename(initialdir="C://", title="Select file",
                                             filetypes=(("csv files", "*.csv"), ("all files", "*.*")))
input_file = root.filename
root.destroy()

print("Please choose the output folder:")

root = tk.Tk()
root.filenames = filedialog.askdirectory(initialdir="C://", title="Select folder")
output_dir = root.filenames
root.destroy()



parse_data(input_file)

filenameList = input_file.split('/')
filename = filenameList[len(filenameList) - 1].split(".")[0]

MDName = output_dir + "//" + filename + "_output.xlsx"

ifile = os.path.basename(input_file)
root, ext = os.path.splitext(ifile)
root = root.replace("output","")

wbMD = xlsxwriter.Workbook(MDName,{'constant_memory': True})
wsMD = wbMD.add_worksheet(root)

writeSheet(wbMD,wsMD,leftData,rightData)

print("...closing output files")
wbMD.close()
print("...DONE\n\n")
