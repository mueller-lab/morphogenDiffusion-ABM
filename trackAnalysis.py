"""This is a script to load the saved tracks.csv files to -

1. make track movies
2. analyze jump sizes and estimate diffusion coefficients
3. make histogram of angle distributions
4. load image grid and select interface and cavities
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from datetime import datetime
import pandas as pd
import os
import time

# define a function to extract binding times and plot the binding time distribution for each ligand
def getBoundFrac(bindXY, trackXY, posArr=None):
    """
    A function to get agent bound-fraction at the end of simulation (overall, in cavity, in interface)
    """
    sp = bindXY.shape
    trackXY=trackXY.astype(int)

    if posArr is not None:
        bindXY2=np.copy(bindXY)
        posArr=posArr.astype(int)
        aSet = set([tuple(x) for x in posArr])
        agCnt = sp[0]
        for i in range(sp[0]):
            pos = tuple(trackXY[i,-1,:])
            if pos not in aSet:
                bindXY2[i,-1]=0
                agCnt-=1
        lastFrac= np.sum(bindXY2[:,-1])/agCnt if agCnt>0 else 0
    else:
        lastFrac= np.sum(bindXY[:,-1])/sp[0]

    return lastFrac

def plotBindTimes(bindingDf, dir, iRng, subG="all", xRng=[0,20], bins =40):
    bindDict={}
    # loop through all the rows to get binding times of all ligands
    bindingArr= np.array(bindingDf)
    agNames = bindingDf.index
    agNames = agNames[:iRng]

    # plot bound fraction over time
    morNames = [a for a in agNames if "morphogen" in a]
    morDf = np.array(bindingDf.loc[morNames])
    morFrac= np.sum(morDf, axis=0)/len(morNames) if len(morNames)>0 else np.sum(morDf, axis=0)
    timeAx = np.arange(0,morDf.shape[1]*0.01,0.01)

    plt.figure()
    plt.xlabel("Time (s)")
    plt.xlim(xRng)
    plt.ylim([0,1])
    plt.ylabel("Bound Fraction")
    plt.plot(timeAx, morFrac, label="morphogen")
    plt.legend()
    plt.savefig(dir/f"BoundFrac_vs_Time.tif")
    plt.close()

    for i in range(iRng):
        row=bindingArr[i]
        boundIntv = []
        counter = 0
        for j in range(len(row)):
            if row[j]==1:
                counter+=1
            else:
                if counter>0: boundIntv.append(counter*0.010)
                counter=0
        if row[-1]==1:
            boundIntv.append(counter*0.010)
        bindDict[agNames[i]]= boundIntv

    #subset binding times of nodal and lefty
    subDict={}
    subDict["morphogen"] = []
    for k,v in bindDict.items():
        k1 = k.index('_')
        subDict[k[:k1]].append(v)

    subDict["morphogen"]= np.array([a for arr in subDict["morphogen"] for a in arr])
    #plot binding time histogram
    for k,v in subDict.items():
        v=np.array(v)
        meanV = np.sum(v)/len(v)
        fig, axs= plt.subplots(2)
        fig.tight_layout()
        axs[0].set_xlim(xRng[0],xRng[1])
        axs[0].set(xlabel="Binding time (s)", ylabel="Count fraction")
        hist, bin = np.histogram(v, bins=bins)
        hist= hist/sum(hist) if sum(hist)>0 else hist
        binW= np.diff(bin)[0]
        axs[0].bar(bin[:-1], hist, width=binW, align='edge',label=k)
        ymx = np.max(hist)
        axs[0].text(0.6*xRng[1],0.2*ymx, f"mean rTm = {meanV:.2f} s")
        axs[0].legend()

        #plot dissociation rates on a log scale
        axs[1].set_xlim(-2,2)
        axs[1].set(xlabel="Dissociation rate log10(1/s)", ylabel="Count fraction")
        hist, bin = np.histogram(np.log10(1/v), bins=bins)
        hist= hist/sum(hist) if sum(hist)>0 else hist
        binW= np.diff(bin)[0]
        # axs[1].set_xscale('log')
        axs[1].bar(bin[:-1], hist, width=binW, align='edge',label=k)
        axs[1].legend()
        plt.subplots_adjust(left=0.12,bottom=0.1)

        plt.savefig(dir/f"bindTimeHist_{subG}_{k}.tif")
        plt.close()
    return meanV

#function to calculate jump distances
def getJumps(tracks, posArr=None, iRng = 6):
    jumpAr= [] #np.zeros((int(trkSp[0]/2), trkSp[1]-1))
    if posArr is not None: aSet = set([tuple(x) for x in posArr])
    for i in range(iRng):
        jumpAr_i = []
        xt =np.array(tracks.iloc[2*i,:])
        yt =np.array(tracks.iloc[2*i+1,:])
        xyAr = np.stack((xt,yt), axis=1)
        for j in range(len(xyAr)-1):
            if posArr is None:
                on1 = 1
            else:
                bSet = set([tuple(x) for x in [xyAr[j],xyAr[j+1]]])
                pos1 = [1 for x in bSet & aSet]
                on1 = 1 if sum(pos1)==2 else 0
            if on1==1:
                x1,y1=xyAr[j]
                x2,y2=xyAr[j+1]
                r = np.sqrt((x2-x1)**2 + (y2-y1)**2)
                if r<200:
                    jumpAr_i.append(r)
        jumpAr.append(np.array(jumpAr_i))
    out = np.array(jumpAr)
    return out

# function to get angles
def getAngles(tracks, posArr=None, iRng = 6):
    angAr=[]
    if posArr is not None: aSet = set([tuple(x) for x in posArr])
    for i in range(iRng):
        angAr_i = []
        xt =np.array(tracks.iloc[2*i,:])
        yt =np.array(tracks.iloc[2*i+1,:])
        xyAr = np.stack((xt,yt), axis=1)
        for j in range(len(xyAr)-2):
            if posArr is None:
                on1 = 1
            else:
                bSet = set([tuple(x) for x in [xyAr[j],xyAr[j+1],xyAr[j+2]]])
                pos1 = [1 for x in bSet & aSet]
                on1 = 1 if sum(pos1)==3 else 0
            if on1==1:
                x1,y1=xyAr[j]
                x2,y2=xyAr[j+1]
                x3,y3=xyAr[j+2]
                th1= np.arctan2((y2-y1),(x2-x1))
                th2= np.arctan2((y3-y2),(x3-x2))
                angTemp= th2-th1 if th2>th1 else (th2-th1+2*np.pi)
                angAr_i.append(angTemp)
        angAr.append(np.array(angAr_i))
    out = np.array(angAr)
    return out

# get average jump size
def getAvgJump(jumpDict):
    avgJmp={}
    for k,v in jumpDict.items():
        s=np.sum(v)
        n=v.size
        avgJmp[k]=0.01*(s/n) #average jump in um
    return avgJmp['morphogen']

#get spot density ratio (Cvt/Ifc)
def spotRto(trackXY, posArr):
    posArr=posArr.astype(int)
    trackXY=trackXY.astype(int)

    aSet = set([tuple(x) for x in posArr])
    sp = trackXY.shape

    cnt=0
    for i in range(sp[0]):
        for j in range(sp[1]):
            pos = tuple(trackXY[i,j,:])
            if pos in aSet:
                cnt+=1
    sptD = cnt/len(aSet)
    return sptD

# a function to plot jump distance histogram
def plotJumps(jumpDict, dir, subG = "all", xRng=[0,1.6], bins =100):
    for k,v in jumpDict.items():
        plt.figure()
        plt.xlim(xRng)
        if subG== "Cvt":
            plt.ylim([0,0.1])
            # plt.ylim([0,0.04])
        else:
            plt.ylim([0,0.1])
        plt.xlabel("Jump distance (\u03BCm)")
        plt.ylabel("Count fraction")
        hist, bin = np.histogram(0.01*v, bins=bins)
        hist= hist/sum(hist)  if sum(hist)>0 else hist
        binW= np.diff(bin)[0]
        plt.bar(bin[:-1], hist, width=binW, align='edge',label=k)
        plt.legend()
        plt.savefig(dir/f"jumpHist_{subG}_{k}.tif")
        plt.close()
# plt.show()

# a function to plot angle histogram
def plotAngles(angDict, dir, subG = "all"):
    N = 36
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=True) #
    bins = [(theta[i], theta[i+1]) for i in range(N-1)]
    theta= np.array([(theta[i]+theta[i+1])/2 for i in range(N-1)])
    width = 2 * np.pi/N #np.pi / 4 * np.random.rand(N)
    for k,v in angDict.items():
        aAr= v
        radii =np.array([ np.sum((a<=aAr) & (aAr<b)) for (a,b) in bins ]) #10 * np.random.rand(N)
        if max(radii)>0:
            radii = radii/max(radii)
            # colors = plt.cm.viridis(np.linspace(0,1,N-1))
            colors = plt.cm.viridis(radii)

            plt.figure()
            ax = plt.subplot(projection='polar')
            ax.bar(theta, radii, width=width, bottom=0.0, color= colors, alpha=0.5, label = k)
            plt.legend()
            plt.savefig(dir/f"angHist_{subG}_{k}.tif")
            plt.close()
        else:
            print("Max(radii)=1, cannot make the plot. ")

# a function to get neighbour positions at radius r
def drawCircle(pos,r):
    r= r+1
    vPosList=[]
    nodes=8
    p0,p1=pos
    for i in range(int(nodes)):
        i=float(i)
        a =p0+r*np.cos(2.*(i+1)*np.pi/nodes)
        b =p1+r*np.sin(2.*(i+1)*np.pi/nodes)
        a= int(round(a))
        b= int(round(b))
        vPosList.append([a,b])
    return vPosList

def getNeighbours(pos, r, gS, any=0):
    # neighbours = np.zeros(8*r)
    if any==1:
        p=pos
        iRan= range(p[0]-r,p[0]+r+1)
        jRan= range(p[1]-r,p[1]+r+1)
        neighbours = np.array([[i,j] for i in iRan for j in jRan])
    else:
        vPosList= drawCircle(pos,r) # [b for a in vPos for b in a]
        neighbours = np.array([a for a in vPosList if 0<a[0]<gS[0] and 0<a[1]<gS[1]])
    return neighbours

# a function to isolate cavity and interface subgrids
def getSubGrids(grid1, nonZeroPos, r=150, itrs=5000):
    print("getting the dictSubGrids[Cvt]...")
    grid = np.copy(grid1)
    gS = grid.shape
    # itrs= int(gS[0]*gS[1]*10**-3)
    dictSubGrids={}
    dictSubGrids["Cvt"] = []
    randIds = np.random.choice([i for i in range(len(nonZeroPos))], itrs)
    randPos= [nonZeroPos[i] for i in randIds]
    adPos = 0
    for p in randPos:
        nE = getNeighbours(p,r, gS)
        cnt = np.array([grid[i,j] for [i,j] in nE])
        goAhd=1
        if 0 in cnt:
            goAhd=0
        if goAhd:
            subGrid = getNeighbours(p,r, gS, any=1)
            dictSubGrids["Cvt"].append(subGrid)
            adPos+=1
    print(f"added {adPos} positions to cvt")
    dictSubGrids["Cvt"] = np.concatenate(dictSubGrids["Cvt"], axis=0)
    dictSubGrids["Cvt"] = np.unique(dictSubGrids["Cvt"], axis=0)
    print("getting the dictSubGrids[Ifc]...")
    dictSubGrids["Ifc"] = np.array(list(set(map(tuple, nonZeroPos))- set(map(tuple, dictSubGrids["Cvt"]))))
    return dictSubGrids

# a function to plot the subgrids
def plotSubGrids(grid1,dictSubGrids,dir):
    gridC= np.copy(grid1)
    for [i,j] in dictSubGrids["Cvt"]: gridC[i,j] = 20
    for [i,j] in dictSubGrids["Ifc"]: gridC[i,j] = 10
    plt.figure()
    plt.imshow(gridC.T, origin= "lower")
    plt.savefig(dir/"grid_cvt_if.tif")
    plt.close()

# modify the jumpAr so that all nodal jumps make one plot and all lefty jumps make another.
"""
jumpAr is an np.array. I'll subset that into three arrays - nodal, lefty, and oep.
I can make use of dictionary data type to do that.

First make a dictionary with keys being ag.name and values being arrays of jumpDist
then, merge keys that have same ag.name.
"""
def makeJumpDict(jumpAr, agList):
    jumpDict={}
    for arr, key1 in zip(jumpAr, agList):
        jumpDict[key1] = arr

    # print(f"jumpDict is {jumpDict}")
    jumpDt={}
    jumpDt["morphogen"]= []

    for k, v in jumpDict.items():
        # jumpDt[k[:-5]].append(v)
        # in case the agent names don't have the same pattern of indexing
        k1 = k.index('_')
        jumpDt[k[:k1]].append(v)

    jumpDt["morphogen"]= np.array([a for arr in jumpDt["morphogen"] for a in arr])
    return jumpDt

############################################
############################################
if __name__=="__main__":
    t1 = time.time()

    cwdPath=Path(os.path.abspath(os.getcwd()))
    # dataPath = cwdPath/"data"/"parScr1"
    dataPath = cwdPath/"data"
    dirList = [f for f in dataPath.iterdir() if f.is_dir()]
    # parFile = str(cwdPath/"parFile_20220213.csv")
    parFile = str(cwdPath/"parFile_20220301.csv")
    parDf = pd.read_csv(parFile, header=0, index_col = None)
    parDf.index = [i for i in range(parDf.shape[0])]

    # print(f"parDf.head is {parDf.head}")
    outList = []
    # inDir = "20220128_132858_2000steps_grid_s8192_4_scaled"
    for inDir in dirList:
        filePath = inDir/"analysis/angHist_Ifc_morphogen.tif"
        trackPath = inDir/"tracks.csv"

        if not filePath.is_file() and trackPath.is_file():
            parID = int(str(inDir)[-3:])
            print(f"analyzing tracks in {inDir}")
            inPath = cwdPath/"data"/inDir

            outPath = inPath/"analysis"
            outPath.mkdir(mode=0o777, parents=True, exist_ok=True)
            # print("made analysis dir")


            # load the tracks.csv
            trackFile = str(inPath/"tracks.csv")
            tracks = pd.read_csv(trackFile, header=0, index_col = 0)

            bindingFile = str(inPath/"binding.csv")
            bindingDf = pd.read_csv(bindingFile, header=0, index_col = 0)

            # remove "oep" tracks and binding data
            # use the .loc method remove rows that contain "oep"
            tracks =tracks.loc[~tracks.index.str.contains('oep'),:]
            bindingDf= bindingDf.loc[~bindingDf.index.str.contains('oep'),:]

            agList = [str(tracks.index[2*i])[:-2] for i in range(int(len(tracks.index)/2))]

            grid1 = np.loadtxt(inPath/'grid.csv',delimiter=',')

            iRng=bindingDf.shape[0]
            # make a numpy array by combining tracks positions with bindingDf
            # get even rows (x pos) and odd rows (y pos) from tracks
            evID = [2*i for i in range(iRng)]
            odID = [2*i+1 for i in range(iRng)]
            tracksX = tracks.iloc[evID,:-1].to_numpy()
            tracksY = tracks.iloc[odID,:-1].to_numpy()
            trkSp = tracksX.shape

            tracksXY = np.dstack([tracksX, tracksY]).reshape(trkSp[0], trkSp[1], 2)
            # print(f"trackXY.shape is {tracksXY.shape}")

            bindXY = bindingDf.to_numpy()


            mean_rTm = plotBindTimes(bindingDf, outPath, iRng, xRng=[0,100], bins= 50)
            nonZeroPos = np.argwhere(grid1)
            dictSubGrids = getSubGrids(grid1, nonZeroPos) # generates a dictionary with radius as key and np array of subgrid positions as values
            plotSubGrids(grid1,dictSubGrids,outPath)

            print("getting bound fractions ...")
            frac = getBoundFrac(bindXY, tracksXY)
            fracC= getBoundFrac(bindXY, tracksXY, dictSubGrids["Cvt"])
            fracI= getBoundFrac(bindXY, tracksXY, dictSubGrids["Ifc"])

            # # calculate jump distances
            print("getting jumps ...")

            jumpAr = getJumps(tracks, iRng = iRng)
            jumpAr_C = getJumps(tracks, dictSubGrids["Cvt"], iRng = iRng)
            jumpAr_I = getJumps(tracks, dictSubGrids["Ifc"], iRng = iRng)
            jumpDt = makeJumpDict(jumpAr, agList)
            jumpDt_C = makeJumpDict(jumpAr_C, agList)
            jumpDt_I = makeJumpDict(jumpAr_I, agList)

            # get the average jump sizes (all, in Cvt, in Ifc)
            jmp= getAvgJump(jumpDt)
            jmpC=getAvgJump(jumpDt_C)
            jmpI=getAvgJump(jumpDt_I)
            dt= 0.01
            avD = jmp**2/(4*dt)
            avD_C= jmpC**2/(4*dt)
            avD_I= jmpI**2/(4*dt)

            #get spot density ratio (cavity/Ifc)
            sptDC=spotRto(tracksXY, dictSubGrids["Cvt"])
            sptDI=spotRto(tracksXY, dictSubGrids["Ifc"])
            sptCbyI=sptDC/sptDI

            outList.append([parID, frac, fracC, fracI, jmp, jmpC, jmpI, sptDC, sptDI, sptCbyI, mean_rTm, avD, \
            avD_C, avD_I])
            # calculate angle distribution
            print("getting angles ...")
            angAr = getAngles(tracks, iRng = iRng)

            angAr_C = getAngles(tracks, dictSubGrids["Cvt"], iRng = iRng)
            #
            angAr_I = getAngles(tracks, dictSubGrids["Ifc"], iRng = iRng)

            angDt = makeJumpDict(angAr, agList)
            angDt_C = makeJumpDict(angAr_C, agList)
            angDt_I = makeJumpDict(angAr_I, agList)

            print("making plots ...")
            # plotRng = [i for i in range(iRng)]
            plotAngles(angDt, outPath)
            plotJumps(jumpDt, outPath)

            plotAngles(angDt_C, outPath,subG = "Cvt")
            plotJumps(jumpDt_C, outPath, subG = "Cvt")
            #
            plotAngles(angDt_I, outPath,subG = "Ifc" )
            plotJumps(jumpDt_I, outPath, subG = "Ifc")
        else:
            print(f"{inDir} is already analyzed.")

    bFracArr= np.array(outList)
    bFracArr=bFracArr[bFracArr[:,0].argsort()]
    # to make the bFracArr same size as parDf
    # bFracArr = np.concatenate((np.zeros((60,14)),bFracArr),  axis=0)

    parDf['bFrac']= bFracArr[:,1]
    parDf['bFracC']= bFracArr[:,2]
    parDf['bFracI']= bFracArr[:,3]
    parDf['jmpAv']=bFracArr[:,4]
    parDf['jmpAvC']=bFracArr[:,5]
    parDf['jmpAvI']=bFracArr[:,6]
    parDf['sptDC']=bFracArr[:,7]
    parDf['sptDI']=bFracArr[:,8]
    parDf['sptCbyI']=bFracArr[:,9]
    parDf['mean_rTm']=bFracArr[:,10]
    parDf['avD'] = bFracArr[:,11]
    parDf['avD_C']=bFracArr[:,12]
    parDf['avD_I']=bFracArr[:,13]


    # modify recepDens column -
    parDf['recepDens'] = np.log10(100/parDf['recepDens'])
    print(f"parDf.head is {parDf.head}")

    parDfName= str(dataPath/'forScatterPlot.csv')
    parDf.to_csv(parDfName, index_label="ID")

    t2= time.time()
    totalSec= t2-t1
    Sec=int(totalSec%60)
    Hrs=int(totalSec//3600)
    Min=int((totalSec%3600)//60)

    print ("Program completed in %sHr:%sMin:%ssec\n"%(Hrs,Min,Sec))
