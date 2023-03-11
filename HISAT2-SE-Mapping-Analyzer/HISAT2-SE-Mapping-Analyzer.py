#!/usr/bin/python3

##############################################################################################
# Author: Mathias Verbeke
# Date of creation: 2022/02/18
# Summary: The HISAT2 Single-End Read Alignment Summary Analyzer is a Python command-line 
# application that takes a HISAT2 single-end read alignment summary file as input and 
# extracts all percentages.  The application then outputs these percentages to a CSV file and, 
# optionally, an Excel file.
##############################################################################################

##################
# Imported modules
##################

import argparse, re, sys, os
import pandas as pd

########################
# Command line arguments
########################

parser = argparse.ArgumentParser()

parser.add_argument('filename', help = 'Provide the name of the HISAT2 single-end read alignment summary file.')
parser.add_argument('-o', '--outdir', required =  False, default = '.', help = 'Provide the name of the output directory.')
parser.add_argument('-e', '--excel', action = 'store_true', help = 'Generate Excel file.')

args = parser.parse_args()

###########
# Variables
###########

filecomp = args.filename.split(".")

if len(filecomp) != 2:
    if len(filecomp) == 1:
        sys.exit("Error: invalid filename: absence of file extension")
    
    if len(filecomp) > 2:
        sys.exit("Error: invalid filename: file contains multiple '.'")

elif filecomp[1] == "csv":
    sys.exit("Error: invalid filename: invalid file extension (csv)")

outname = f"{args.outdir.rstrip('/')}/{os.path.basename(filecomp[0])}.csv"
excelname = f"{args.outdir.rstrip('/')}/{os.path.basename(filecomp[0])}.xlsx"

#############################
# Processing file information
#############################

infile = open(args.filename, 'r')
outfile = open(outname, 'w')
outfile.write("SampleNo,% unpaired,% aligned 0 times,% aligned exactly 1 time,% aligned >1 times,% overall alignment rate\n")

sample = 1
flag = "NOK"
for line in infile.readlines():

    line = line.strip()
    container = line.split()
    percentage = re.sub('[(,),%]', '', container[1])

    if "unpaired" in line or "aligned 0 times" in line or "aligned exactly 1 time" in line or "aligned >1 times" in line:
        if "unpaired" in line:
            outfile.write(f"{sample},")
        
        outfile.write(f"{percentage},")

    elif "overall alignment rate" in line:
        percentage = re.sub('[(,),%]', '', container[0])
        outfile.write(f"{percentage}\n")
        sample += 1

infile.close()
outfile.close()


if args.excel == True:
    # Opening the csvfile again
    with open(outname) as csvfile:

        # Reading the csv file
        csvDataframe = pd.read_csv(csvfile)

        # Creating output Excel file
        with pd.ExcelWriter(excelname) as excel:
        
            # Converting the csv file to an Excel file
            csvDataframe.to_excel(excel, index = False)