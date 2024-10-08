"""

author: Ujjwal Mukund Mahajan

**************************************************
**              Deconvolution of images         **
**************************************************

"""
#@File(label = "Output directory", style = "directory") output


import os, json
from ij import IJ, ImagePlus
from ij import WindowManager


outputDir = output.getAbsolutePath()

## HE staining
boost = IJ.run("Create RGB Template", "");
IJ.run(boost, "Colour Deconvolution", "vectors=H&E hide");
IJ.run("Create CD-LUT", "channel1=RGB24bit-(Colour_1) channel2=RGB24bit-(Colour_2) channel3=RGB24bit-(Colour_3)");
lut = WindowManager.getCurrentImage();
IJ.saveAs(lut, "Tiff", os.path.join(outputDir, "lut"));
IJ.run("Close All");
## H-AEC
boost = IJ.run("Create RGB Template", "");
IJ.run(boost, "Colour Deconvolution", "vectors=[H AEC] hide");
IJ.run("Create CD-LUT", "channel1=RGB24bit-(Colour_1) channel2=RGB24bit-(Colour_2) channel3=RGB24bit-(Colour_3)");
lut = WindowManager.getCurrentImage();
IJ.saveAs(lut, "Tiff", os.path.join(outputDir, "lutAEC"));
IJ.run("Close All");
## PAS-AB
boost = IJ.run("Create RGB Template", "");
IJ.run(boost, "Colour Deconvolution", "vectors=[User values] [r1]=0.43347 [g1]=0.85503 [b1]=0.284647 [r2]=0.7086714 [g2]=0.20601025 [b2]=0.6747923 [r3]=0.66848 [g3]=0.703754 [b3]=0.240529 hide");
IJ.run("Create CD-LUT", "channel1=RGB24bit-(Colour_1) channel2=RGB24bit-(Colour_2) channel3=RGB24bit-(Colour_3)");
lut = WindowManager.getCurrentImage();
IJ.saveAs(lut, "Tiff", os.path.join(outputDir, "lutAB_PAS"));
IJ.run("Close All");
## AMEC-AP
boost = IJ.run("Create RGB Template", "");
IJ.run(boost, "Colour Deconvolution", "vectors=[FastRed FastBlue DAB] hide");
IJ.run("Create CD-LUT", "channel1=RGB24bit-(Colour_1) channel2=RGB24bit-(Colour_2) channel3=RGB24bit-(Colour_3)");
lut = WindowManager.getCurrentImage();
IJ.saveAs(lut, "Tiff", os.path.join(outputDir, "lutCostain"));
IJ.run("Close All");
## PSFG
boost = IJ.run("Create RGB Template", "");
IJ.run(boost, "Colour Deconvolution", "vectors=[User values] [r1]=0.59227383 [g1]=0.3264422 [b1]=0.7366459 [r2]=0.099971585 [g2]=0.73738605 [b2]=0.6680326 [r3]=0.7995107 [g3]=0.5913521 [b3]=0.10528667 hide");
IJ.run("Create CD-LUT", "channel1=RGB24bit-(Colour_1) channel2=RGB24bit-(Colour_2) channel3=RGB24bit-(Colour_3)");
lut = WindowManager.getCurrentImage();
IJ.saveAs(lut, "Tiff", os.path.join(outputDir, "lutPS"));
IJ.run("Close All");
