#!/usr/bin/env python
# coding: utf-8

# In[6]:


import csv
import numpy as np
import pandas as pd

# In[7]:


def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))


# In[8]:


def BedOverlap(bed1, bed2, overlapPercent):
    p0, p1 = 0, 0
    output = []
    
    ##### TODO change to raise error message
    if(not(all((isinstance(x[0], str)) and (isinstance(x[1], int)) and (isinstance(x[2], int))) for x in bed1)): return "Invalid Input in Bed1"
    if(not(all((isinstance(x[0], str)) and (isinstance(x[1], int)) and (isinstance(x[2], int))) for x in bed2)): return "Invalid Input in Bed2"
    
    while p0 < len(bed1) and p1 < len(bed2):
        a0, b0, a1, b1 = int(bed1[p0][1]), int(bed1[p0][2]), int(bed2[p1][1]), int(bed2[p1][2])
        if (a0<b0<a1<b1):
            p0+=1
        elif (a1<b1<a0<b0):
            p1+=1
        elif (a0<=a1<=b0<=b1):
            chrOverlap = overlap(a0, b0, a1, b1) # calculate the overlap between the current entry in bed1 and bed2
            bed1Perc = (chrOverlap / (b0 - a0) * 100) # determine the percentage of overlap with bed1
            if bed1Perc >= overlapPercent:
                tempList = []
                tempList.append(bed1[p0])
                tempList.append(bed2[p1])
                tempList.append(round(bed1Perc)) # append overlap percentage to the current line in bed1
                output.append(tempList) # append output
            p0 +=1
        elif a0<=a1<=b1<=b0:
            chrOverlap = overlap(a0, b0, a1, b1) # calculate the overlap between the current entry in bed1 and bed2
            bed1Perc = (chrOverlap / (b0 - a0) * 100) # determine the percentage of overlap with bed1
            if (bed1Perc >= overlapPercent):
                tempList = []
                tempList.append(bed1[p0])
                tempList.append(bed2[p1])
                tempList.append(round(bed1Perc)) # append overlap percentage to the current line in bed1
                output.append(tempList) # append output
            p1 +=1
        elif a1<=a0<=b0<=b1:
            tempList = []
            tempList.append(bed1[p0])
            tempList.append(bed2[p1])
            
            ### TODO confirm if this should be integer?
            tempList.append('100') 
            output.append(tempList)
            p0 +=1
        elif a0<=a1<=b1<=b0:
            chrOverlap = overlap(a0, b0, a1, b1) # calculate the overlap between the current entry in bed1 and bed2
            bed1Perc = (chrOverlap / (b0 - a0) * 100) # determine the percentage of overlap with bed1
            if (bed1Perc >= overlapPercent):
                tempList = []
                tempList.append(bed1[p0])
                tempList.append(bed2[p1])
                tempList.append(round(bed1Perc)) # append overlap percentage to the current line in bed1
                output.append(tempList) # append output
            p1 +=1
        elif a1<=a0<=b1<=b0:
            chrOverlap = overlap(a0, b0, a1, b1) # calculate the overlap between the current entry in bed1 and bed2
            bed1Perc = (chrOverlap / (b0 - a0) * 100) # determine the percentage of overlap with bed1
            if (bed1Perc >= overlapPercent):
                tempList = []
                tempList.append(bed1[p0])
                tempList.append(bed2[p1])
                tempList.append(round(bed1Perc)) # append overlap percentage to the current line in bed1
                output.append(tempList) # append output
            p1+=1
        #else:
            #break
            
    ## TODO: output list
    return np.array(output,dtype=object)
        
    


# In[9]:


def BedScan(fileName, delim):
    output = []
    BedFile = open(fileName)#open the file
    BedReader = csv.reader(BedFile, delimiter=delim) #intiate the reader
    for line in BedReader:
        output.append(line)
    BedFile.close
    output.pop(0)
    out = list((row[0], int(row[1]), int(row[2])) for row in output)
    return out


# In[10]:


def CompBeds():
    delim = "\t"
    print("Enter file name for bed 1: ")
    bed1Name = input()
    print("Enter file name for bed 2: ")
    bed2Name = input()
    bed1 = BedScan(bed1Name, delim)
    bed2 = BedScan(bed2Name, delim)
    print("enter minimum overlap: ")
    minOverlap = int(input())
    result = BedOverlap(bed1, bed2, minOverlap)
    print(result)
    #This assumes that the bed files are sorted, perhaps using this command sort -k2 -n file > out