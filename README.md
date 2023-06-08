# CIRCOAST: a statistical hypothesis test for cellular colocalization with network structures

This software project provides a statistical hypothesis test (CRICOAST) for changes in cellular colocalization with the network structures such as vascular networks, and provide a suite of tools to study cellular colocalization with a variety of models (ARCAS).


--------------------------------------------------------------------------------
User Directions
--------------------------------------------------------------------------------

* Note: Example images found in the "test_images" folder in this directory.
* Note2: Experiment data download links in [ARCAS wiki](https://github.com/uva-peirce-cottler-lab/ARCAS/wiki)

## Startup
1. Open the file "USER_INITIALIZE.m" in Matlab.
2. Run this file within the matlab editor (Either hit Run icon in editor or hit F5).
3. Wait a couple of seconds for the program to pop up. 

## Basic GUI Operation
1. Click "Load Image" and choose an image to analyze.
2. Set the appropriate parameters for the image
    1. Number of Trials: total number of Monte Carlo trials to run in simulation.You can leave at 100,000, that is sufficient.
    2. Cells Per Trial: number of cells in the image that is being observed for Colocalization. Set this to the observed number of cells in this image.
    3. Microns per Pixel: resolution of the image. Used to determine cell size in pixel units.
    4. Cell diameter (um): diameter of the cell body being measured for colocalization. This parameter is determined by manually measuring a significant sample of cells and calculating average area of cell body.
3. If the image is a binary image, the program will display "Image is Binary" in the upper left portion of the window. You may then proceed to step 7.
4. If the image is not binary, the program will display "Image is Not Binary" in
   the upper left portion of the window.
    1. You need to threshold the image by activating a channel for thresholding. Thresholding controls are the buttons and sliders on the right side of the program window.
    2. Activate the channel that contains the features to be thresholded (such channel for vasculature/ EC cells) by selecting the R,G, or B checkbox at the bottom right hand side of the window.
    3. Set a threshold value for the channel by dragging the slider for the same channel that was activated until there is a nice threshold. The threshold will be displayed live as the slider is adjusted.
    4. You may hide the other channels by clicking on the R,G, or B boxes on the top right of the window.
5. Hit "Run Single Trial" button and the output of a single trial will be displayed.
6. Verify that the simulated cells (shown in red in main image) is roughly the size you would expect.
7. Hit "Run Simulation" for all of the requested trials to be run.
8. The mean colocalization fraction and standard deviation will be displayed at the bottom of the window for the Monte Carlo Model of Random Placement and Binomial Model of Random Placement.
    
## Advance Operation
The actual CIRCOAST statistics are conducted in datasets of images with two study groups. Input data required for these tests include (1) the thresholded (binary) images, and (2) a csv file that assigns each binary image to a group, along with the observed total number of COIs (cells of interest such as injected cells, etc.) and colocalizing cells.

From the ARCAS GUI, select Process>Process Experiment a from a CSV file that contains the required information (please see example data for syntax). 
The program will calculate:
1. CELLVCOAV: measure colcoalization versus random placement with a binomial hypothesis test for _each individual image_.
2. 1S-CIRCOAST: measure colocalization versus random placement for an entire _study group of images_.
3. 2S-CIRCOAST: measure coloczation differences between _two study groups_.

p values are printed to the matlab command line and results saved in the same directory as the input CSV file.
    
