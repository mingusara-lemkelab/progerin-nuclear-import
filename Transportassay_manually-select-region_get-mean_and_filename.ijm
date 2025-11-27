// Activate the first image in the stack
selectImage(1);

// Enable the "Polygon Selection" tool
//setTool("polygon");
// Measure the selected region in the current image
run("Set Measurements...", "area mean display redirect=None decimal=3 pixel");
     
    

// Ask the user to draw a polygon region in the first image
waitForUser("Draw a polygon region in the first image to select the region of interest. Press 'OK' when done.");

// Get the coordinates of the drawn polygon
getSelectionCoordinates(xPoints, yPoints);

// Loop through all images in the stack
for (i = 1; i <= nSlices; i++) {
    // Set the current slice
    setSlice(i);
        
    // Measure the selected region in the current image
    run("Measure");
    
    // Wait for a short moment to allow manual adjustments before moving to the next slice
    //wait(100);
}

// Close the selection
//run("Select None");

// Activate the first image in the stack
selectImage(1);