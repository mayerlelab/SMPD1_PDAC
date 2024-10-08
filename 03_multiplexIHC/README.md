# Multiplex IHC Analysis

## QuPath
   tested for QuPath-0.2.3-m

1. **Preprocessing**

    1. Create QuPath Project

    1. Import all the images per TMA / staining to a project. Slide name should be "slideID_staining.mrxs". eg. TMA1_CD45.mrxs

    1. Create TMA grid using TMA dearrayer on refenrece Image.

    1. Import TMA labels from .txt files (TMA labels will be unique patient IDs)

1. **TMA grid to Annotations**

    run [tmaGrid2Annotations.groovy](03_multiplexIHC/scripts/quPath/tmaGrid2annotations.groovy) script.

1. **Align all Images in project using affine transformation**

    1. Open reference image

    1. Open Interactive image alignement and arrange moving Image using **shift + E** and running afine transformation followed by update

    1. Run [interactiveAlignmentMatrix.groovy](03_multiplexIHC/scripts/quPath/interactiveAlignmentMatrix.groovy) script adding name of moving image and save transformed Matrix.

    1. Repeat step 2 and 3 for all moving images.

    1. Do not close Interactive image alignement box for next step.

1. **Export all aligned annotations to all images in a project**

    1. Add String refStain = "referenece Image staining" (eg. CD45) to [alignmentAfterTransformation.groovy](03_multiplexIHC/scripts/quPath/alignmentAfterTransformation.groovy).
    
    1. Select slide other that reference slide without closing alignment box.

    1. Run [alignmentAfterTransformation.groovy](03_multiplexIHC/scripts/quPath/alignmentAfterTransformation.groovy) for a project.

    1. Close Interactive image alignement box.

    1. Notes:  Check the code if you need it. : # import qupath.ext.align.gui.ImageServerOverlay -> for QuPAth version 3.0.0. 
    2. If you have Qupath 3.0.0 or more, you download the extended version. Otherwise you have error. So deactivate the upper code with // , run the lower one with extension.

1. **Save all annotations as Tiff**

    1. run [annotations2Tiff.groovy](03_multiplexIHC/scripts/quPath/annotations2tiff.groovy for aproject

    1. all Tiff Images will be saved in raw_cores folder within project. Images name siwll be in format of unique "PatientID_Staining.tff".

## Fiji and cell profiler
  tested for 2.1.0/1.53c 
  
Run all fiji analysis using [multiplexIHC.sh](03_multiplexIHC/multiplexIHC.sh) then processed with for [Image Cytometry](03_multiplexIHC/image_cytometry/01_Image_Cytometry.Rmd), [data merging](03_multiplexIHC/image_cytometry/02_data_merging.Rmd) and anndata preparation.


Cleaned Image cytometry analyzed results then processed for [multiplex_IHC analysis](03_multiplexIHC/image_cytometry/04_multiplex_IHC_analysis.Rmd). 

Patients stratification for short- and long-term survivors was performed using [auton-survival](https://github.com/autonlab/auton-survival.git). 
