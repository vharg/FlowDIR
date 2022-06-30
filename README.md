## FlowDIR 2022
*Elly Tennant, Susanna Jenkins, Sébastien Biass*

FlowDIR is a MATLAB tool for forecasting the travel directions of topographically controlled hazardous flows. This page is meant as a quick start guide, for more information and an understanding of how FlowDIR works please refer to Tennant et al., (in prep)... 
## Getting started


FlowDIR was written using MATLAB v9.12. Before starting ensure that you have the following MATLAB toolboxes installed: mapping; image processing, parallel computing.

<ol> 

<li> Download and unpack the zipfile into your MATLAB folder ensuring that the paths are set correctly. The folder contains all of the files needed to run FlowDIR along with several case studies. The content of the folder is as follows: 

```
FlowDIR
   ├── Code
   |	 └── FlowDIR.m
   ├── DEMs
   |	 └── Shinmoedake_2016_5m_clip.tif
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


<li> FlowDIR can be run from the command line as follows:

`FlowDIR('Shinmoedake', 'Shinmoedake_2016_5m_clip.tif',800, 678155, 3532081, 150,20, 500, 1, 10)`

Type <code>help FlowDIR</code> into the MATLAB command window to learn more about the inputs required for command line executable mode. 

Alternatively, to run with the dedicated graphical user interface (GUI) simply type  <code>FlowDIR</code> into the command window and the following will pop up:
<br/>


<center><img src="https://github.com/EllyTennant/FlowDir/blob/main/images/GUI.png" width="250"></center>

<li> The input parameters are defined as:

|  FlowDIR input    | Description | Suggested range|
| ----------- | ----------- | ----------- |
| DEM file      | This should be a .tif format file and projected into the UTM coordinate system       |-|
| Starting coordinate	| A single point in UTM coordinates (X,Y).|-|
| Swath length  (m) | Should be long enough to extend from the starting coordinate outside of the crater in all directions. | 500 – 1000 m (default 800 m)|
| Buffer (m)     | Swaths are clipped to the maximum elevation and a buffer is applied so that the swath extends outside of the crater. For breached craters a buffer towards the higher end of the range may be required.      | 50 – 150 m (no default)|
|   Elevation threshold (m)  | Total elevation change along swath above which the flow is not expected to overtop the crater. Should be based on knowledge of the current crater morphology (maximum height from the base to the rim) along with past erupted volumes, and flow rheologies.   | 20 – 50 m (no default)|
| Maximum number of steps allowed   | The maximum number of steps allowed in the path of steepest descent.  | 200-600 (default 500)|
|  Capture uncertainty in start? (0/1)   | Input 1 to run FlowDIR with uncertainty in the starting coordinate. When uncertainty = 0, the initialisation polygon edges are 1 DEM cell width from the starting coordinate, and the polygon vertices are used as additional initialisation points such that there are 5 in total. When uncertainty = 1, the user can set the size of the polygon to increase the number of initialisation points. | Default 0|
|  Start uncertainty (m)  | When uncertainty = 1, this is used to set the size of the polygon. | No default|



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

<li> When FlowDIR is finished, all of the figures along with the workspace containing all variables, will be saved into the <code>out/VolcanoName/X</code> folder.

<br/>
<img src="https://github.com/EllyTennant/FlowDir/blob/main/images/Shinmoedake_ex.png" width="800">

## Citation
FlowDIR was published in ....

