# progerin-nuclear-import
Sample data and codes for the analysis of capsid nuclear import in progerin-expressing cells. 

First, open .ims file with FIJI, load the Transportassay-manually-select-region_get-mean-and-filename.ijm into FIJI, and manually select regions inside the nucleus. Play the timelapse each time to make sure the cells do not move, and the region remains inside the nucleus in all time frames (importantly, that nuclear envelope does not enter the region throughout the time lapse). Save the imaes as .csv file. 
After each region is measured, erase the selected region by pressing erase (it erases from the last frame, channel 3, where the code ends). then, duplicate this image (Frame 45, CH3) and save as tiff. to have a copy of the regions which you have measured. You can use this .tiff file for another analysis later. 

Next, run the Analysis_Matlab code into matlab, 
