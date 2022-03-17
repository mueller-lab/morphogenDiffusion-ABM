## morphogenDiffusion-ABM

The program simulates an agent-based model of morphogen diffusion in extracellular space on a two-dimensional grid.

Developer- Amit Landge (amit.landge@uni-konstanz.de)  
Requires- Python 3.7.12  
Dependencies- listed in require.yml  

# Usage
```
morphogenDiffusionABM.py [-h] [-n NUM [NUM ...]]
                                [-gS GRIDSIZE [GRIDSIZE ...]] [-st STEPS]
                                [-nrw NARROWINTERFACE]
                                [-slr SAILOR [SAILOR ...]] [-lv LIVE]
                                [-uImg USEIMG] [-iName IMGNAME]
                                [-parScr PARAMETERSCREEN] [-parID PARAMETERID]
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

- Example  
```python -u morphogenDiffusionABM.py -n 200 0 0 0 100 -st 100 -nrw 1 -slr 0 4 -uImg Yes -iName s8192_4_scaled 1> out.log 2>&1```

# trackAnalysis.py
Generates a directory named 'analysis' in each simulation output directory and stores plots and .csv files after analysis.

- Example  
```python -u trackAnalysis.py 1>out1.log 2>&1```

# plot3D_scatter.py
Generates a 3D scatterplot of outputs for all tested parameter conditions using forScatterPlot.csv generated by trackAnalysis.py as input.

- Example  
```python -u plot3D_scatter.py 1>out2.log 2>&1```
