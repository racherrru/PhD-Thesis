path = getDirectory("Choose a Directory "); 
outputFolder = getDirectory("Choose_a_Directory "); 
filelist = getFileList(path); 
 for (i=0; i< filelist.length; i++) {
     // process vsi files only
     if (endsWith(filelist[i], ".vsi")) {
     run("Viewer", "open=["+ path + filelist[i] +"]");
	run("Scale Bar...", "width=100 height=147 thickness=4 font=30 color=Black background=None location=[Lower Right] horizontal bold overlay");
	saveAs("PNG", path + filelist[i] + ".png");
	close();
}
 }
