/*
 * author: Theresa Weltermann and Ujjwal Mukund Mahajan
 *
 * **************************************************
 * **            artifacts substraction            **
 * **************************************************
 *
 */

#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", value = ".tif") suffix

var outputDir
var inputDir

outputDir = File.getName(output);
inputDir = File.getName(input);

setBatchMode(true);
processFolder(input);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix) && (indexOf(list[i], "_R1_") != -1))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
		   // import artifacts
           run("Image Sequence...", "open=["+ input + File.separator + list[i] +"] type=RGB starting=2 sort");
           name = getTitle();
           print("----------------------------------------------");
           print("Processing: ",name);
           print("----------------------------------------------");
           titleStack = name + "_" + "Stack";
           titleMask = name + "_" + "Mask";
           rename(titleStack);
		   // import masks
           run("Image Sequence...", "open=["+ input + list[i] +"] type=RGB number=1 sort");
           run("8-bit");
           run("Invert LUT");
           setOption("BlackBackground", true);
           run("Convert to Mask");
           makeRectangle(70, 60, 250, 150);
           run("Measure");
           val = getResult("Mean");
           run("Select None");
           if (val < 10) {
           run("Invert LUT");
           }
           run("Erode");
           rename(titleMask);
           // image calculator addition
           imageCalculator("Max create stack", titleStack, titleMask);
           run("Image Sequence... ", "format=TIFF start=1 digits=4 use save=["+ output + File.separator +"]");
           run("Close All");
}

close("Results");
