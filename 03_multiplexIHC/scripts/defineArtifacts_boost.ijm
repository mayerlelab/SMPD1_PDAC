/*
* author: Theresa Weltermann and Ujjwal Mukund Mahajan
*
* **************************************************
* ** Batch processing code to define artifacts    **
* **************************************************
*
*/


#@File(label = "Input directory", style = "directory") input
#@File(label = "InputLut directory", style = "directory") inputLut
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", description = ".tif") suffix
#@String(label = "mask") mask

var outputDir
var inputDir
var inputLutDir

outputDir = File.getName(output);
inputDir = File.getName(input);
inputLutDir = File.getName(inputLut);

setBatchMode(true);
processFiles(input);

function processFiles(input) {
list = getFileList(input);
for (i = 0; i < list.length; i++) {
if(endsWith(list[i], suffix))
processFile(input, inputLut, output, list[i]);
}
}

// Artifact extraction
function processFile(input, inputLut, output, file) {
open(input + File.separator + file);
title = getTitle();

stringColor1 = "CD-0001";
stringColor2 = "CD-0002";
stringColor3 = "CD-0003";

print("----------------------------------------------");
print("Processing: ",title);
print("----------------------------------------------");

// Define Tissue boundaries
if (title.matches(".*"+mask+".*")) {
print(mask);    
run("8-bit");
run("Auto Threshold", "method=Default setthreshold");
//run("Threshold...");
setThreshold(245,255);
//setOption("BlackBackground", false);
run("Convert to Mask");
run("Gaussian Blur...", "sigma=15");
setMinAndMax(80, 125);
run("Apply LUT");
saveAs("Tiff", output + File.separator + file);
run("Close All");
} else if (title.matches(".*_HE.*")) {
open(inputLut + File.separator + "lutHE.tif");
lutName = getTitle();
run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
run("Stack to Images", "");
selectImage(stringColor3);
run("Auto Threshold", "method=Default setthreshold");
//run("Threshold...");
setThreshold(0, 140);
setOption("BlackBackground", false);
run("Convert to Mask");
saveAs("Tiff", output + File.separator + file);
run("Close All");
} else if (mask.matches(".*_SR.*")==false) {
open(inputLut + File.separator + "lutHE.tif");
lutName = getTitle();
run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
run("Stack to Images", "");
selectImage(stringColor3);
run("Auto Threshold", "method=Default setthreshold");
//run("Threshold...");
setThreshold(0, 100);
setOption("BlackBackground", false);
run("Convert to Mask");
saveAs("Tiff", output + File.separator + file);
run("Close All");
} else if (mask.matches(".*AB.*")==false) {
open(inputLut + File.separator + "lutAB_PAS.tif");
lutName = getTitle();
run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
run("Stack to Images", "");
selectImage(stringColor3);
run("Auto Threshold", "method=Default setthreshold");
//run("Threshold...");
setThreshold(0, 100);
setOption("BlackBackground", false);
run("Convert to Mask");
saveAs("Tiff", output + File.separator + file);
run("Close All");
} else {
open(inputLut + File.separator + "lutAEC.tif");
lutName = getTitle();
run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
run("Stack to Images", "");
selectImage(stringColor3);
run("Auto Threshold", "method=Default setthreshold");
//run("Threshold...");
setThreshold(0, 120);
setOption("BlackBackground", false);
run("Convert to Mask");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
}
