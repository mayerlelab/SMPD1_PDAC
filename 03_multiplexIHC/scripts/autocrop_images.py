"""
author: Ujjwal Mukund Mahajan

 **************************************************
 ** Batch processing code for autocroping images **
 **************************************************

"""


#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", value = ".tif") suffix


import os
import ij
from ij import IJ, ImagePlus
from ij import WindowManager as WM

def run():
# define path
  inputDir = input.getAbsolutePath()
  outputDir = output.getAbsolutePath()

  for root, directories, filenames in os.walk(inputDir):
    filenames.sort();
    for fileName in filenames:

      # check for file extension
      if not fileName.endswith(suffix):
        continue
      process(inputDir, outputDir,fileName)

def process(inputDir, outputDir, fileName):
	print "----------------------------------------------"
	print "Processing: ", fileName
	print "----------------------------------------------"
        # open image

        image = IJ.openImage(os.path.join(inputDir, fileName));
        # autocrop
        IJ.run(image, "Auto Crop (guess background color)", "");
        # resize
        IJ.run(image, "Size...", "width=2685 height=2685 depth=1 constrain average interpolation=Bicubic");
        # enhance contrast
        IJ.run(image,"Enhance Contrast...", "saturated=0.2 process_all");
        # Save image
        IJ.saveAs(image, "Tiff", os.path.join(outputDir, fileName));
        image.close();

run()
