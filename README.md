# MorphogenDiffusion-ABM

The program simulates an agent-based model of morphogen diffusion in extracellular space on a two-dimensional grid.

Developer- Amit Landge (amit.landge@uni-konstanz.de)  
Requires- Python 3.7.12  
Dependencies- listed in require.yml  

## Installation
 - Install conda

 Please follow the instructions given [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install Anaconda or Miniconda for your system.  

 - Create and activate a conda environment

 ```conda create --name envMD python=3.7 numpy scipy matplotlib pandas pathlib```

 ```conda activate envMD```

 - Clone the git repository  

To install git, please run  
 ```conda install -c anaconda git```  

 Then, clone the repository using -  
 ```git clone https://github.com/mueller-lab/morphogenDiffusion-ABM.git```

Typical install time ~ 5 min. 

## Usage

## morphogenDiffusionABM.py
```
morphogenDiffusionABM.py [-h] [-n NUM [NUM ...]]
                                [-gS GRIDSIZE [GRIDSIZE ...]] [-st STEPS]
                                [-nrw NARROWINTERFACE]
                                [-slr SAILOR [SAILOR ...]] [-lv LIVE]
                                [-uImg USEIMG] [-iName IMGNAME]
                                [-parScr PARAMETERSCREEN] [-parID PARAMETERID]
```
```
optional arguments:
  - h, --help            show this help message and exit
  - n NUM [NUM ...], --num NUM [NUM ...]
                        Number of agents. default=1
  - gS GRIDSIZE [GRIDSIZE ...], --gridSize GRIDSIZE [GRIDSIZE ...]
                        Size of the 2D grid. default=[100,100]
  - st STEPS, --steps STEPS
                        Number of simulation steps. default=500
  - nrw NARROWINTERFACE, --narrowInterface NARROWINTERFACE
                        Steps to reduce interface width. default=5
  - slr SAILOR [SAILOR ...], --sailor SAILOR [SAILOR ...]
                        Types of sailor (list). default=0
  - lv LIVE, --live LIVE
                        Whether to show plots during runtime or not.
                        default=No
  - uImg USEIMG, --useImg USEIMG
                        Whether to use given image to make the grid.
                        default=No
  - iName IMGNAME, --imgName IMGNAME
                        Name of the image file. default=s8192_4_scaled.tiff
  - parScr PARAMETERSCREEN, --parameterScreen PARAMETERSCREEN
                        Whether to perform a parameter screen. default=No
  - parID PARAMETERID, --parameterID PARAMETERID
                        Row number in the parameter file. default=0
```

- Example  
```python -u morphogenDiffusionABM.py -n 200 0 0 0 100 -st 100 -nrw 1 -slr 0 4 -uImg Yes -iName s8192_4_scaled```  
The output is stored in ./data/yyyymmdd_HHMMSS_*/. Note that for each simulation a separate output directory is created.

Runtime for this command = 0Hr:6Min:18sec, on a laptop with Intel(R) Core(TM) i7-8850H CPU @ 2.60GHz

Memory requirement - 2GB

## makeParFile.py
Generates a parameter file named 'parFile_20220301.csv'. This file contains parameter combinations that are used to run a parameter
screen for morphogen diffusion in extracellular space using morphogenDiffusionABM.py. Each row is a unique parameter combination. The
row numbers needs to be provided as input to morphogenDiffusionABM.py to run a simulation with that parameter combination.

- Example  
```python -u makeParFile.py```

After generating a parFile_20220301.csv, it can be used to set parameters for the simulations.

- Example  
```python -u morphogenDiffusionABM.py -n 200 0 0 0 100 -st 100 -nrw 1 -slr 0 4 -uImg Yes -iName s8192_4_scaled -parScr Yes -parID 0```


## trackAnalysis.py
Generates a directory named 'analysis' in each simulation output directory and stores plots and .csv files after analysis. Before running this
many simulations need to be completed. Each simulation with generate a separate output directory to store the results.

- Example  
```python -u trackAnalysis.py```

## plot3D_scatter.py
Generates a 3D scatterplot of outputs for all tested parameter conditions using forScatterPlot.csv generated by trackAnalysis.py as input.

- Example  
```python -u plot3D_scatter.py```
