
ARCAS: Automated Random Cell Association Simulator

--------------------------------------------------------------------------------
User Directions
--------------------------------------------------------------------------------

Note: Example images found in the "test_images" folder in this directory.

1. Open the file "USER_INITIALIZE.m" in Matlab.
2. Run this file within the matlab editor (Either hit Run icon in editor or hit
   F5).
3. Wait a couple of seconds for the program to pop up. click "Load Image" and 
   choose an image to analyze.
4. Set the appropriate parameters for the image
    a. Number of Trials: total number of Monte Carlo trials to run in 
       simulation.You can leave at 100,000, that is sufficient.
    b. Cells Per Trial: number of cells in the image that is being observed for
       Colocalization. Set this to the observed number of cells in this image.
    c. Microns per Pixel: resolution of the image. Used to determine cell size
       in pixel units.
    d. Cell diameter (um): diameter of the cell body being measured for 
       colocalization. This parameter is determined by manually measuring a 
       significant sample of cells and calculating average area of cell body.
5. If the image is a binary image, the program will display "Image is Binary" in
   the upper left portion of the window. You may then proceed to step 7.
6. If the image is not binary, the program will display "Image is Not Binary" in
   the upper left portion of the window.
    a. You need to threshold the image by activating a channel for thresholding.
       Thresholding controls are the buttons and sliders on the right side of
       the program window.
    b. ACtivate the channel that contains the features to be thresholded (such 
       channel for vasculature/ EC cells) by selecting the R,G, or B checkbox
       at the bottom right hand side of the window.
    c. Set a threshold value for the channel by dragging the slider for the same
       channel that was activated until there is a nice threshold. The threshold
       will be displayed live as the slider is adjusted.
    d. You may hide the other channels by clicking on the R,G, or B boxes on the
       top right of the window.
7. Hit "Run Single Trial" button and the output of a single trial will be
   displayed.
8. Verify that the simulated cells (shown in red in main image) is roughly the 
   size you would expect.
9. Hit "Run Simulation" for all of the requested trials to be run.
10. The mean colocalization fraction and standard deviation will be displayed at
    the bottom of the window.
    