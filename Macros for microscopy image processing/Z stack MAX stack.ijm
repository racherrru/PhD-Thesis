path = getDirectory("Choose a Directory "); 
outputFolder = getDirectory("Choose_a_Directory "); 
filelist = getFileList(path); //load array of all files inside input directory
 for (i=0; i< filelist.length; i++) {
     // process oir files only
     if (endsWith(filelist[i], ".oir")) {
         
         // open file, requires LOCI tools (aka Bio-Formats)
     run("Bio-Formats Importer", "open=["+ path + filelist[i] +"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		
		//run z project analysis, to project based on average intensity
		run("Z Project...", "projection=[Max Intensity]");

//run("Channels Tool...");
		Stack.setDisplayMode("composite");

		run("Stack to RGB");
		run("Scale Bar...", "width=100 height=5 font=18 color=White background=None location=[Lower Right] bold");
		saveAs("TIF", path + filelist[i] + "RGB.tif");
		close();
		close("*");
		     }
 }
