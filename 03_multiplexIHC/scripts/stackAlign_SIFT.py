"""
 author: Ujjwal Mukund Mahajan

 **************************************************
 **  scale-invariant feature transform (SIFT)    **
 **             for image alignment              **
 **************************************************

"""

#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output
#@String(label = "Number of Images") nuImages
#@String(label = "File suffix", value = ".tiff") suffix


import os
import sys
from glob import glob
from ij import IJ, ImagePlus, WindowManager

inputDir = input.getAbsolutePath()
outputDir = output.getAbsolutePath()

def process_sequence(folder, output, files):
	name = os.path.dirname(folder)
	subDirName = name.split("/")[-1]
	newName = subDirName + '_'

	string = ["Aligned", str(nuImages), "of", str(nuImages)]
	stack = " ".join(string)


	print "----------------------------------------------"
	print "Processing folder: ", subDirName
	print "----------------------------------------------"

        IJ.run("Image Sequence...", "open=" + folder + " sort use");
        WindowManager.getCurrentImage();
        IJ.run("Linear Stack Alignment with SIFT", "initial_gaussian_blur=1 steps_per_scale_octave=3 minimum_image_size=64 maximum_image_size=1024 feature_descriptor_size=4 feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.92 maximal_alignment_error=25 inlier_ratio=0.05 expected_transformation=Rigid interpolate");
        WindowManager.getImage(stack);
        IJ.run("Image Sequence... ", "format=TIFF name=[" + newName + "] start=1 digits=4 use save=["+ output +"]");
        IJ.run("Close All");

def process_folder(folder, output, suffix):

	for subdir in glob(os.path.join(folder, '*/')):
		process_folder(subdir, output, suffix)

	files = glob(os.path.join(folder, '*' + suffix))
	if len(files) > 0:
		process_sequence(folder, output, files)

process_folder(inputDir, outputDir, suffix)
os._exit(0)
