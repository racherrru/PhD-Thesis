path = getDirectory("Choose a Directory "); 
outputFolder = getDirectory("Choose_a_Directory "); 
filelist = getFileList(path); //load array of all files inside input directory

run("Bio-Formats Macro Extensions");
processBioFormatFiles(path);

function processBioFormatFiles(path) {
	for (i = 0; i < filelist.length; i++){
	if(endsWith(filelist[i], "lif")){
		Ext.setId(path + filelist[i])'
		Ext.getSeriesCount(seriesCount);
		
			run("Bio-Formats Importer", "open=["+ path + filelist[i] +"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+seriesCount);
			name=getTitle();
         	run("Duplicate...", "duplicate");
         	stack = getImageID(); //original
         	close("\\Others");
         	
         	selectImage(stack);
         run("Duplicate...", "title=melanosomes duplicate channels=2"); //NKI
         selectImage(stack);
         run("Duplicate...", "title=basal duplicate channels=3"); // basal    
         selectImage(stack);
         run("Duplicate...", "title=suprabasal duplicate channels=4"); //suprabasal
         
         imageCalculator("AND create", "melanosomes", "basal");
         run("Duplicate...", "title=basal_melanin");
         run("Auto Threshold", "method=Yen white");
         setThreshold(1, 255);
         setOption("BlackBackground", true);
         run("Convert to Mask");
         run("Divide...", "value=255");
         run("Measure");
         
         imageCalculator("AND create", "suprabasal", "basal");
         run("Duplicate...", "title=basal_suprabasal");
         imageCalculator("AND create", "basal_suprabasal", "melanosomes");
         run("Duplicate...", "title=suprabasal_melanin");
         run("Auto Threshold", "method=Yen white");
         setThreshold(1, 255);
         setOption("BlackBackground", true);
         run("Convert to Mask");
         run("Divide...", "value=255");
         run("Measure");


        close("*"); //closing all images in loop cycle
	
		}
	}
}
saveAs("Results", outputFolder+"results"+".csv");
close("Results");       

        
