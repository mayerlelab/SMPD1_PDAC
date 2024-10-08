"""
 author: Ujjwal Mukund Mahajan
 **************************************************
 **                image z stacks                **
 **************************************************
"""

#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", value = ".tif") suffix
#@String(label = "Number of Images", value = "nuImages") nuPseudo

import os
import sys
from glob import glob
from ij import IJ, ImagePlus, WindowManager
from ij.plugin import ZProjector, ImageCalculator

inputDir = input.getAbsolutePath()
outputDir = output.getAbsolutePath()

def process_sequence(folder, output, files):
	name = os.path.dirname(folder)
	subDirName = name.split("/")[-1]

	print "----------------------------------------------"
	print "Processing folder: ", subDirName
	print "----------------------------------------------"

        image = IJ.run("Image Sequence...", "open=[" + folder + "] number=["+ str(nuPseudo) +"] sort use");
        IJ.run(image, "Subtract Background...", "rolling=150 disable stack");
        #IJ.run(image, "Kuwahara Filter", "sampling=5 filter stack");
        IJ.run(image, "Remove Outliers...", "radius=5 threshold=50 which=Bright stack");
        IJ.run(image, "Gamma...", "value=1.30 stack");
        image2 = IJ.run(image, "Grouped Z Project...", "projection=[Max Intensity] group=["+str(nuPseudo)+"]");
        IJ.saveAs(image2, "Tiff", os.path.join(outputDir, subDirName));
        IJ.run("Close All");

def process_folder(folder, output, suffix):

	for subdir in glob(os.path.join(folder, '*/')):
		process_folder(subdir, output, suffix)

	files = glob(os.path.join(folder, '*' + suffix))
	if len(files) > 0:
		process_sequence(folder, output, files)



process_folder(inputDir, outputDir, suffix)

#os._exit(0)
