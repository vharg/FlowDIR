## FlowDIR 2021
*Elly Tennant, Susanna Jenkins, Sebastien Biass*

FlowDIR is a MATLAB tool for forecasting the travel directions of topographically controlled hazardous flows. This page is meant as a quick start guide, for more information and an understanding of how FlowDIR works please refer to Tennant et al., (...). 
## Getting started


FlowDIR was written using MATLAB 2021a. Before starting ensure that you have the following MATLAB toolboxes installed: mapping; image processing, parallel computing.

<ol> 

<li> Download and unpack the zipfile into your MATLAB folder ensuring that the paths are set correctly. The folder contains all of the files needed to run FlowDIR along with several case studies. The content of the folder is as follows: 

```
FlowDIR
   ├── Code
   |	 └── FlowDIR.m
   ├── DEMs
   |	 └── Shinmoedake_2016_15m.tif
   ├── Dependancies
   |	 ├── topotoolbox-master
   |	 ├── invprctile.m
   |	 ├── polarwitherrorbar.m
   |	 └── polarPcolor.m
   └── Out
         └── Shinmoedake
         		├── 0
         		└── 1
```


<li> FlowDIR can be run from the command line or through dedicated graphical user interfaces. Type ` help FlowDIR` into the MATLAB command window to learn more about the inputs required for command line executable mode. To run with the GUIs simply type `FlowDIR` into the command window and the following GUI will pop up:
<br/>

<center><img src="https://github.com/EllyTennant/FlowDir/blob/main/images/GUI.png" width="250"></center>

<li> The input parameters are defined as:

|  FlowDIR input    |  |
| ----------- | ----------- |
| DEM file      | This should be a .tif format file and projected into the UTM coordinate system       |
| Default swath length   | This is the maximum distance away from the starting point that may be considerd in the calculation. FlowDIR generates swaths extending from the start point to this length. The highest elevation along each swath is identified and swaths are clipped to this and a buffer is applied. This value should extend well beyond the crater edge.        |
| Buffer (m)     | Buffer applied to the detected crater perimeter. The perimeter + this buffer will be the area included in the calculation. For dome summit topographies the maximum elevation along swath will be the summit, and swaths will extend from the start point to this buffer length.      |
|   Elevation threshold   | This is a threshold used for crater overtopping and is not relevant to dome type topographies. FlowDIR calculates the accumulated elevation change along the swath. Direction bins are identified as having along swath elevation gains above or below this threshold. Directions above this threshold cannot be overtopped.   |
| Maximum number of steps allowed   | Used in stage 2 of the calculation (POSD). FlowDIR calculates the least cost path through the elevation matrix. This is the maximum number of cells that can be traversed in this calculation.  |
|  Capture uncertainty in start? (0/1)   | Binary option to include uncertainty in the user defined start point. When capture uncertainty in start = 1, the initialisation polygon is enlarged and simulations are run using secondary start points within the polygon in addition to the initial starting coordinate. The size of the polygon is determined by the start uncertainty. |
|  Start uncertainty (m)  | Here we define the distance to which uncertainty is included, for example if the DEM resolution is 10 m and 20 m is selected here, simulations are initiated from the 2 cells adjacent to the selected start point in each direction in addition to the start point, such that 17 unique simulations are conducted. |



<li>FlowDIR will plot the DEM, and ask you to clip it to an area of interest by drawing a polygon. +/- at the top right hand side of the plot can be used to zoom in or out.
<br/>


<center><img src="https://github.com/EllyTennant/FlowDir/blob/main/images/clip_dem.png" width="400"></center>

The new zoomed in DEM is then plotted and you can click to define the start point.

<li> Congratulations FlowDIR is now running, its progress is stated as below:

```Running FlowDir, please wait...
Start point # [1/17] ...
Start point # [2/17] ...
Start point # [3/17] ...
Start point # [4/17] ...
...
Finished
```

<li> When FlowDIR is finished, all of the figures along with the workspace containing all variables, will be saved into the `out/VolcanoName/X` folder.

<br/>
<img src="https://github.com/EllyTennant/FlowDir/blob/main/images/Shinmoedake_ex.png" width="800">

## Citation
FlowDIR was published in ...

