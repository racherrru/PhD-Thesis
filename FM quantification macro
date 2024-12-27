path = getDirectory("Choose a Directory "); 
outputFolder = getDirectory("Choose_a_Directory "); 
filelist = getFileList(path); 

 for (i=0; i< filelist.length; i++) {
 	
     // process vsi files only
     if (endsWith(filelist[i], ".vsi")) {
     	run("Bio-Formats Importer", "open=["+ path + filelist[i] +"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
         name=getTitle();
         run("Duplicate...", "duplicate");
         stack = getImageID(); //original
         close("\\Others");
     
     	//run("Channels Tool...");
		run("Stack to RGB");
		run("Colour Deconvolution2", "vectors=[FastRed FastBlue DAB] output=8bit_Transmittance simulated cross hide");

		
		selectImage(filelist[i] + " - macro image-1 (RGB)-(Colour_1)");
		//run("Threshold...");
		setThreshold(0, 150, "raw");
		run("Analyze Particles...", "size=8000-Infinity circularity=0.00-0.05 show=Overlay display clear summarize add composite");
	    roiManager("Select",0);
		roiManager("Add");
		roiManager("Measure");
		roiManager("Save", outputFolder + filelist[i]+ "roi_selection_counter"+".zip");
	
	   
	    selectImage(filelist[i] + " - macro image-1 (RGB)-(Colour_2)");
	    //run("Threshold...");
		roiManager("Select",0);
		roiManager("Add");
		
		setThreshold(0, 175, "raw");
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1 show=Overlay display clear summarize add composite");
		roiManager("Measure");
		roiManager("Save", outputFolder + filelist[i]+ "roi_selection_mel"+".zip");

		roiManager("Select",0)
		roiManager("Delete");
		close("*");
		    
     }
 Table.save(path + "Summary.csv", "Summary");
 }
 

