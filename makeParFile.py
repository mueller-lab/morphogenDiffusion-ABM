"""Run this code to make a parameter file for nodal_ABM simulation. A .csv file is generated with each row being a parameter
combination for one simulation.

The name of the parFile needs to be updated in the morphogenDiffusionABM.py file to use the parameters from this file for the
simulation, and `-parScr Yes` should be provide as user argument while running morphogenDiffusionABM.py
"""
import numpy as np


# list of parameters
parNames = ["BS","nrw","recepDens", "rTm"]

# make lists of possible parameter values
parDict={}
parDict["BS"]=[2]
parDict["nrw"]=[0]
parDict["recepDens"]=[200]
parDict["rTm"]= [1,16]
kList=list(parDict)

# permute the parameters to generate an array of all parameter combinations
parArr= np.zeros((2,4))

parList= [a for a in parDict.values()]
ids=[len(a) for a in parList]

def makeListComb(list):
    listLen= [len(a) for a in list]
    colNum= len(listLen)
    rowNum= 1
    for a in listLen:
        rowNum=rowNum*a
    print(f"rowNum is {rowNum}")
    combArr= np.zeros((rowNum,colNum))
    #now go through each list within the list and modify column in the combArr
    for i, val in enumerate(list):
        #first repeat val (rowNum/(product of lengths uptil now)) times
        rptNum=1
        for j in range(colNum-1,i,-1):
            rptNum=rptNum*listLen[j]
        print(f"rptNum is {rptNum}")
        rptVal= np.repeat(val, rptNum)
        # then tile and reshape the val to equal to the rowNum

        tileVal = np.tile(rptVal, (int(rowNum/len(rptVal)), 1))
        colVal= np.reshape(tileVal, (rowNum))
        #modify i'th column in combArr
        combArr[:,i]=colVal

    return combArr

combArr=makeListComb(parList)

#save the array as a csv file
headStr = ",".join(kList)
np.savetxt('parFile_20220301.csv', combArr, fmt='%i',delimiter=',', header=headStr, comments='')

print("made parFile.csv")
