/*
* author: Ujjwal Mukund Mahajan and Theresa Weltermann
*
* **************************************************
* **         Deconvolution of images              **
* **************************************************
*
*/

#@File(label = "Input directory", style = "directory") input
#@File(label = "InputLut directory", style = "directory") inputLut
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", value = ".tif") suffix
#@String(label = "Co staining", value = "costain") costain
#@String(label = "Co staining1", value = "costain1") costain1
#@String(label = "Co staining2", value = "costain2") costain2
#@String(label = "HE Co staining", value = "heCostain") heCostain

var outputDir
var inputDir
var inputLutDir

outputDir = File.getName(output);
inputDir = File.getName(input);
inputLutDir = File.getName(inputLut);

setBatchMode(true);
processFiles(input, inputLut);

function processFiles(input, inputLut) {
list = getFileList(input);
for (i = 0; i < list.length; i++)  {
if(endsWith(list[i], suffix))
processFile(input, inputLut, output, list[i]);
}
}

function processFile(input, inputLut, output, file) {
open(input + File.separator + file);
title = getTitle();

stringColor1 = "CD-0001";
stringColor2 = "CD-0002";
stringColor3 = "CD-0003";

print("----------------------------------------------");
print("Processing: ",title);
print("----------------------------------------------");
//Hematoxylin & Eosin
if (matches(title, ".*_HE.*")) {
	nameHE = replace(title, "HE", "Eosin");
	open(inputLut + File.separator + "lutHE.tif");
	lutName = getTitle();
	run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
	run("Stack to Images", "");
	selectImage(stringColor1);
	run("Brightness/Contrast...");
	setMinAndMax(50, 130);
	// run("Log");
	run("Apply LUT");
	run("Smooth");
	run("Subtract Background...", "rolling=50 light disable");
	// run("Kuwahara Filter", "sampling=5");
	run("RGB Color");
	saveAs("Tiff", output + File.separator + title);
	selectImage(stringColor2);
	run("Brightness/Contrast...");
    setMinAndMax(65, 175);
    // run("Log");
    run("Apply LUT");
    run("Smooth");
    run("Subtract Background...", "rolling=50 light disable");
    // run("Kuwahara Filter", "sampling=5");
    rename(nameHE);
    run("RGB Color");
    saveAs("Tiff", output + File.separator + nameHE);
	run("Close All");
}
else if (matches(title, ".*_"+heCostain+".*")) {
	nameHE = replace(title, heCostain, "HE");
	open(inputLut + File.separator + "lutAEC.tif");
	lutName = getTitle();
	run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
	run("Stack to Images", "");
	selectImage(stringColor1);
	run("Brightness/Contrast...");
	setMinAndMax(100, 200);
	// run("Log");
	run("Apply LUT");
	run("Smooth");
	run("Subtract Background...", "rolling=50 light disable");
	// run("Kuwahara Filter", "sampling=5");
	rename(nameHE);
	run("RGB Color");
	saveAs("Tiff", output + File.separator + nameHE);
	selectImage(stringColor2);
	run("Brightness/Contrast...");
	setMinAndMax(150, 200);
	// run("Log");
	run("Apply LUT");
	run("Smooth");
	run("Subtract Background...", "rolling=50 light disable");
	// run("Kuwahara Filter", "sampling=5");
	rename(title);
	run("RGB Color");
	saveAs("Tiff", output + File.separator + title);
	run("Close All");
}
else if (matches(title, ".*_AB.*")) {
  nameAlcB = replace(title, "AB-PAS", "AB");
  open(inputLut + File.separator + "lutAB_PAS.tif");
  lutName = getTitle();
  run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
  run("Stack to Images", "");
  selectImage(stringColor3);
  run("Brightness/Contrast...");
  setMinAndMax(80, 150);
  run("Smooth");
  // run("Log");
  // run("Apply LUT");
  // run("Subtract Background...", "rolling=50 light disable");
  // run("Kuwahara Filter", "sampling=5");
  rename(nameAlcB);
  run("RGB Color");
  saveAs("Tiff", output + File.separator + nameAlcB);
  selectImage(stringColor1);
  namePAS = replace(title, "AB-PAS", "PAS");
  run("Brightness/Contrast...");
  setMinAndMax(100, 200);
  // run("Log");
  run("Apply LUT");
  run("Smooth");
  run("Subtract Background...", "rolling=50 light disable");
  // run("Kuwahara Filter", "sampling=5");
  rename(namePAS);
  run("RGB Color");
  saveAs("Tiff", output + File.separator + namePAS);
  run("Close All");
}
else if (matches(title, ".*_"+costain+"_.*")) {
  //nameAMEC = replace(title, costain, costain+".AMEC");
  open(inputLut + File.separator + "lutCostain.tif");
  lutName = getTitle();
  run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
  run("Stack to Images", "");
  // selectImage(stringColor3);
  // run("Brightness/Contrast...");
  // setMinAndMax(80, 150);
  // run("Log");
  // run("Apply LUT");
  // run("Smooth");
  // run("Subtract Background...", "rolling=50 light disable");
  // run("Kuwahara Filter", "sampling=5");
  // rename(nameAMEC);
  // run("RGB Color");
  // saveAs("Tiff", output + File.separator + nameAMEC);
  selectImage(stringColor2);
  // nameAP = replace(title, costain, costain+".AP");
  run("Brightness/Contrast...");
  setMinAndMax(80, 180);
  // run("Log");
  run("Apply LUT");
  run("Smooth");
  run("Subtract Background...", "rolling=50 light disable");
  // run("Kuwahara Filter", "sampling=5");
  // rename(nameAP);
  run("RGB Color");
  saveAs("Tiff", output + File.separator + title);
  run("Close All");
}
else if (matches(title, ".*_"+costain1+"_.*")) {
  // nameAMEC = replace(title, costain1, costain1+".AMEC");
  open(inputLut + File.separator + "lutCostain.tif");
  lutName = getTitle();
  run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
  run("Stack to Images", "");
  selectImage(stringColor3);
  run("Brightness/Contrast...");
  setMinAndMax(20, 140);
  // run("Log");
  run("Apply LUT");
  run("Smooth");
  run("Subtract Background...", "rolling=50 light disable");
  // run("Kuwahara Filter", "sampling=5");
  // rename(nameAMEC);
  run("RGB Color");
  saveAs("Tiff", output + File.separator + title);
  run("Close All");
}
else if (matches(title, ".*_"+costain2+"_.*")) {
  // nameAMEC = replace(title, costain2, costain2+".AMEC");
  open(inputLut + File.separator + "lutCostain.tif");
  lutName = getTitle();
  run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
  run("Stack to Images", "");
  selectImage(stringColor3);
  run("Brightness/Contrast...");
  setMinAndMax(20, 150);
  // run("Log");
  run("Apply LUT");
  run("Smooth");
  run("Subtract Background...", "rolling=50 light disable");
  // run("Kuwahara Filter", "sampling=5");
  // rename(nameAMEC);
  run("RGB Color");
  saveAs("Tiff", output + File.separator + title);
  run("Close All");
}
// picrosirius red fast green
else if (matches(title, ".*_SR.*")){
	saveAs("Tiff", output + File.separator + title);
	run("Close All");
}
else {
	open(inputLut + File.separator + "lutAEC.tif");
	lutName = getTitle();
	run("Apply CD-LUT", "image="+title+" lut="+lutName+"");
	run("Stack to Images", "");
	selectImage(stringColor2);
	run("Brightness/Contrast...");
	setMinAndMax(50, 160);
	// run("Log");
	run("Apply LUT");
	run("Smooth");
	run("Subtract Background...", "rolling=50 light disable");
	// run("Kuwahara Filter", "sampling=5");
	rename(title);
	run("RGB Color");
	saveAs("Tiff", output + File.separator + title);
	run("Close All");
}
}
