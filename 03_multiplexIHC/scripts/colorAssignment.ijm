/*
* author: Ujjwal Mukund Mahajan and Theresa Weltermann
*
* **************************************************
* **               Colour assingment              **
* **************************************************
*
*/

#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", value = ".tif") suffix
#@String(label = "Costaining name") costain
#@String(label = "Costaining name1") costain1
#@String(label = "Costaining name2") costain2
#@String(label = "stain1costain") stain1costain
#@String(label = "stain2costain") stain2costain
#@String(label = "stain1costain1") stain1costain1
#@String(label = "stain2costain1") stain2costain1
#@String(label = "stain1costain2") stain1costain2
#@String(label = "stain2costain2") stain2costain2
#@String(label = "stain1 Name costain") stain1Namecostain
#@String(label = "stain2 Name costain") stain2Namecostain
#@String(label = "stain1 Name costain1") stain1Namecostain1
#@String(label = "stain2 Name costain1") stain2Namecostain1
#@String(label = "stain1 Name costain2") stain1Namecostain2
#@String(label = "stain2 Name costain2") stain2Namecostain2
#@String(label = "mask name") mask


var outputDir
var inputDir

outputDir = File.getName(output);
inputDir = File.getName(input);

setBatchMode(true);
processFiles(input);

function processFiles(input) {
list = getFileList(input);
for (i = 0; i < list.length; i++) {
if(endsWith(list[i], suffix))
processFile(input, output, list[i]);
}
}

function setLUTByWavelength(wavelength) {
	// These values will be between 0 and 1
	red = green = blue = 0;

	if (wavelength < 380 || wavelength > 780)
		abort("Only wavelengths between 380 and 780 are supported");
	else if (wavelength <= 390) {
		red = (440 - wavelength) / (440 - 380);
		blue = 1;
	}
	else if (wavelength <= 490) {
		green = (wavelength - 440) / (490 - 440);
		blue = 1;
	}
	else if (wavelength <= 510) {
		green = 1;
		blue = (510 - wavelength) / (510 - 490);
	}
	else if (wavelength <= 580) {
		red = (wavelength - 510) / (580 - 510);
		green = 1;
	}
	else if (wavelength <= 645) {
		red = 1;
		green = (645 - wavelength) / (645 - 580);
	}
	else
		red = 1;

	intensity = 1;
	if (wavelength > 700)
		intensity = 0.3 + 0.7 * (780 - wavelength) / (780 - 700);
	else if (wavelength < 420)
		intensity = 0.3 + 0.7 * (wavelength - 380) / (420 - 380);

	red *= intensity;
	green *= intensity;
	blue *= intensity;

	print(red+" "+green+" "+blue);
	// assuming gamma == 1
	// setForegroundColor(red * 255, green * 255, blue * 255);
	reds = newArray(256); 
    greens = newArray(256); 
    blues = newArray(256);
    
    for (i=0; i<256; i++) {
        reds[i] = round(i*red);
        greens[i] = round(i*green);
        blues[i] = round(i*blue);
    }
    setLut(reds, greens, blues);
}


function processFile(input, output, file) {
open(input + File.separator + file);
title = getTitle();
print("----------------------------------------------");
print("Processing: ",title);
print("----------------------------------------------");
if (matches(title, ".*HE.*")) {
run("8-bit");
run("Invert");
run("Grays");
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R1_.*")==true && matches(title, ".*HE.*")==false) {
run("8-bit");
run("Invert");
setLUTByWavelength(450);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R2_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(460);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R3_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(465);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R4_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(475);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R5_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(510);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R6_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(550);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R7_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(570);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R8_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(580);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R9_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(590);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R10_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(610);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_R11_.*")) {
run("8-bit");
run("Invert");
setLUTByWavelength(630);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*_"+costain+"_.*")==true && matches(title,".*_"+stain1costain+"_.*")==true) {
stainName = title.replace(costain+"."+stain1costain,stain1Namecostain);
run("8-bit");
run("Invert");
setLUTByWavelength(650);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + stainName);
run("Close All");
}
else if (matches(title, ".*_"+costain+"_.*")==true && matches(title, ".*_"+stain2costain+"_.*")==true) {
stainName = title.replace(costain+"."+stain2costain,stain2Namecostain);
run("8-bit");
run("Invert");
setLUTByWavelength(700);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + stainName);
run("Close All");
}
else if (matches(title, ".*_"+costain1+"_.*")==true && matches(title,".*_"+stain1costain1+"_.*")==true) {
stainName = title.replace(costain1+"."+stain1costain1,stain1Namecostain1);
run("8-bit");
run("Invert");
setLUTByWavelength(640);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + stainName);
run("Close All");
}
else if (matches(title, ".*_"+costain1+"_.*")==true && matches(title, ".*_"+stain2costain1+"_.*")==true) {
stainName = title.replace(costain1+"."+stain2costain1,stain2Namecostain1);
run("8-bit");
run("Invert");
setLUTByWavelength(680);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + stainName);
run("Close All");
}
else if (matches(title, ".*_"+costain2+"_.*")==true && matches(title,".*_"+stain1costain2+"_.*")==true) {
stainName = title.replace(costain2+"."+stain1costain2,stain1Namecostain2);
run("8-bit");
run("Invert");
setLUTByWavelength(670);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + stainName);
run("Close All");
}
else if (matches(title, ".*_"+costain2+"_.*")==true && matches(title, ".*_"+stain2costain2+"_.*")==true) {
stainName = title.replace(costain2+"."+stain2costain2,stain2Namecostain2);
run("8-bit");
run("Invert");
setLUTByWavelength(690);
run("RGB Color");
run("Subtract Background...", "rolling=200");
run("Smooth", "");
//run("Enhance Contrast...", "saturated=0.2");
saveAs("Tiff", output + File.separator + stainName);
run("Close All");
}

else if (matches(title, ".*"+mask+".*")) {
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*AB.*") && matches(title, mask)==false) {
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
else if (matches(title, ".*SR.*") && matches(title, mask)==false) {
saveAs("Tiff", output + File.separator + file);
run("Close All");
}
}