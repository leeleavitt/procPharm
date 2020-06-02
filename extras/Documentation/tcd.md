# traceClick.dev
This is an interactive visualization, which is constantly under development (thus the dev). In general to use the function follow,

## Loading experiment data
Loading the experiment is as easy as double clicking the `RD.experiment.Rdata` within the folder. Performing this action with automatically set the correct working directory.

If you would like to circumvent this, it is essential you set the correct working directory, then load the `RD.experiment.Rdata`. Below is an example, after opening the R console, or an R terminal.
```
# Load in the package
require(procPharm)

# Set the working directory (this is where the experiment is located)
# If you have downloaded the experiment from the paper, check the download folder.
setwd("Y:/Lee Leavitt/Presentations/Paper/Figures/TTA A2 Dose Response/190531.39.m.m3.p2 ttaa2 potassium dose response comparison/")

# Load the experiment. 
load("./RD.190531.39.m.m3.p2.Rdata")
```
For an overview of the composition of the experiment data, go [here](../../extras/Documentation/RDcomp.md).

## tcd
Now that both the software and the experiment data is loaded invoking the `tcd()` function. This function has flexible inputs and many features to use during visualization as well as analysis tasks.

Most consoles are interactive. To view what has populated the current work space us the `ls()` function. The output of this function will show the experiment name. 
```
ls()
```
prior to using tcd familiarize yourself with the controls
```R
?tcd

```



To enter the experiment for visualization, (the gif shows the initial view entering the visualization),
```
tcd(RD.190531.39.m.m3.p2)
```
![](../../extras/gifWT/openTCD.gif)


## Main controls
<kbd>up</kbd>: Moves through current selected list of cell

<kbd>down</kbd>: Moves down through selected list of cells

<kbd>1-9, -, =</kbd>: Collect cells into groups 1-12

<kbd>o</kbd> Order traces in single view based on **c.dat**, **bin**, **scp**. For information on these data frames please go [here](./RDcomp.md)

<kbd>F7</kbd>: Populate the groups with the cell_types

<kbd>shift</kbd><kbd>p</kbd>: Select group of cells to populate the groups

<kbd>s</kbd>: Stack all cells from the selected group

<kbd>shift</kbd><kbd>s</kbd>: Sample stack selected group. Usefule if selected group has more than 50 cells.

![](../../extras/gifWT/tcdIntroduction.gif)

## Cell View

Additionally the cells can be viewed using,

<kbd>v</kbd> View the cells from the selected group

<kbd>shift</kbd><kbd>i</kbd> Change the image for the multi-view and the single-view images. Only one image can be selected here.

<kbd>i</kbd> select the image to display next to the stacked traces. Can be multiple

![](../../extras/gifWT/tcdImage.gif)

## Advanced Trace view

Visualizing traces has many more options.

<kbd>t</kbd> select the type of trace to view

<kbd>h</kbd> select the color of the trace

<kbd>shift</kbd><kbd>d</kbd> Change the separation of the traces **< 1** closer together **> 1** further apart

<kbd>shift</kbd><kbd>o</kbd> Sort the traces on a continuous selected variable.

<kbd>r</kbd> Rename a selected group to what you would like.

<kbd>shift</kbd><kbd>v</kbd> Select the values to appear on the right side of the trace.

<kbd>u</kbd> Underline toggle for both stacked traces and single trace view

<kbd>l</kbd> Select the window regions to display. Any or all can be observed.

In the example below,
   * A new group is chosen using <kbd>shift</kbd><kbd>p</kbd>
   * The baseline corrected trace is selected <kbd>t</kbd>
   * The color of the trace is selected <kbd>h</kbd>
   * The underline is removed <kbd>u</kbd>
   * The traces are squished closer together <kbd>shift</kbd><kbd>d</kbd>
   * More values are added to the right of trace using <kbd>shift</kbd><kbd>v</kbd>
   
   ![](../../extras/gifWT/advancedTraceFunctioning.gif)


## Custom Statistics
Sorting traces based on non obvious characteristics is important. To capture these characteristics two functions have been created to make these statistics. These functions are encapsulated in the F1 and F2 key.

<kbd>F1</kbd> 
    
This function allows you to select various window regions to create a ratio. During the experiments we apply compounds which elicit amplification, block, or direct effects. Rapidly sorting the traces based on these effects is paramount for a rapid analysis. 

The function additionally has a "boxPlotSelector" functionality. This means. one can use a box plot to select cells based on the new statistic. 

The user has the option to save the statistic. This appends the newly created statistic to the end of the **scp** data frame with the user input name. Once this statistic is saved the user has the option of sorting cells based on this stat. 



