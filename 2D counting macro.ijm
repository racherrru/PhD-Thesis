path = getDirectory("Choose a Directory "); 
outputFolder = getDirectory("Choose_a_Directory "); 
filelist = getFileList(path); // Load array of all files inside input directory

for (i = 0; i < filelist.length; i++) {
    // Process .oir files only
    if (endsWith(filelist[i], ".oir")) {
        // Open file with Bio-Formats
        run("Bio-Formats Importer", "open=[" + path + filelist[i] + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
        name = getTitle();
        run("Duplicate...", "duplicate");
        stack = getImageID(); // Original
        close("\\Others");
        
        // Max Intensity projection
        run("Z Project...", "projection=[Max Intensity]");
        run("Duplicate...", "title=max duplicate");
        replacedFileName = replace(filelist[i], " ", "_");

        c1 = replacedFileName + "-C1";
        c2 = replacedFileName + "-C2";
        c3 = replacedFileName + "-C3";
        
        selectImage("max");
        run("Duplicate...", "title="+ c1 +" duplicate channels=1"); // DAPI
        selectImage("max");
        run("Duplicate...", "title="+ c2 +" duplicate channels=2"); // K14    
        selectImage("max");
        run("Duplicate...", "title="+ c3 +" duplicate channels=3"); // Melanosomes

        // Count all cells present
        selectWindow(c1); // Nuclei (C1)
        run("Threshold");
        setAutoThreshold("Default dark no-reset"); // Adjust thresholds as needed
        run("Analyze Particles...", "size=20-Infinity show=Overlay include summarize");
        Table.save(path + "Summary.csv", "Summary")

        // Create binary masks for C2 and C3 to find common cells
        selectWindow(c2);
        run("Convert to Mask"); // Convert C2 to binary mask
        run("Duplicate...", "title=C2_mask");
        
        selectWindow(c3);
        run("Convert to Mask"); // Convert C3 to binary mask
        run("Duplicate...", "title=C3_mask");

        imageCalculator("AND create", "C2_mask","C3_mask");
        
        rename(replacedFileName + "C3-C2-Common");
        selectImage(replacedFileName + "C3-C2-Common");
        run("Analyze Particles...", "clear include summarize");
        Table.save(path + "Summary.csv", "Summary")
        
        imageCalculator("Subtract create", "C3_mask","C2_mask");
        rename(replacedFileName + "C3-C2-Subtract");
        selectImage(replacedFileName + "C3-C2-Subtract");
        run("Analyze Particles...", "clear include summarize");
                Table.save(path + "Summary.csv", "Summary")
        
       close("*");	
    }

}
