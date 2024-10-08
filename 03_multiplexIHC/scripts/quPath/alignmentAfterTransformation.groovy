/**
 * *******************************************************************************************
 * *********         Script to export pixels & annotations for TMA images        *************   
 * *******************************************************************************************
 *
 * Adapted from Mike Nelson's post (https://forum.image.sc/t/qupath-multiple-image-alignment-and-object-transfer/35521/2) to work on transformation matrices created from the alignment of multiple target slides onto reference slides in a single QuPath project.
 * This script assumes WSI filenames are in the format: slideID_stain.fileExt
 *
 * Requires creating each affine transformation from the target images so that there are multiple transform files with different names.
 *
 * The labelled image can also optionally use indexed colors to depict the colors of the
 * original classifications within QuPath for easier visualization & comparison.
 *
 * @author Ujjwal Mukund Mahajan (modified from @author Yau Mun Lim and Mike Nelson) 
 *
 * (https://forum.image.sc/t/qupath-multiple-image-alignment-and-object-transfer/35521/2)
 * (https://forum.image.sc/t/calling-image-alignment-function-from-a-script/35617/10)
 */

 
// Delete existing objects
def deleteExisting = true

// Change this if things end up in the wrong place
def createInverse = false

// Specify reference stain
String refStain = "Cers6"

import qupath.lib.objects.PathCellObject
import qupath.lib.objects.PathDetectionObject
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.objects.PathTileObject
import qupath.lib.roi.RoiTools
import qupath.lib.roi.interfaces.ROI


import java.awt.geom.AffineTransform

import static qupath.lib.gui.scripting.QPEx.*

// Affine folder path
path = buildFilePath(PROJECT_BASE_DIR, 'Affine')

// Get list of all images in project
def projectImageList = getProject().getImageList()

new File(path).eachFile{ f->
    f.withObjectInputStream {
        matrix = it.readObject()

        def targetFileName = f.getName()
        def (targetImageName, imageExt) = targetFileName.split('\\.')
        def (slideID,targetStain) = targetImageName.split('_')

        def targetImage = projectImageList.find {it.getImageName() == targetFileName}
        if (targetImage == null) {
            print 'Could not find image with name ' + f.getName()
            return
        }
        def targetImageData = targetImage.readImageData()
        def targetHierarchy = targetImageData.getHierarchy()

        refFileName = slideID + "_"  + refStain + "." + imageExt
        def refImage = projectImageList.find {it.getImageName() == refFileName}
        def refImageData = refImage.readImageData()
        def refHierarchy = refImageData.getHierarchy()

        def pathObjects = refHierarchy.getAnnotationObjects()
        
        print 'Aligning objects from reference slide ' + refFileName + ' onto target slide ' + targetFileName

        // Define the transformation matrix
        def transform = new AffineTransform(
                matrix[0], matrix[3], matrix[1],
                matrix[4], matrix[2], matrix[5]
        )
        if (createInverse)
            transform = transform.createInverse()
            
        if (deleteExisting)
            targetHierarchy.clearAll()
            
        def newObjects = []
        for (pathObject in pathObjects) {
            newObjects << transformObject(pathObject, transform)
            
        }
        
        targetHierarchy.addPathObjects(newObjects)
        targetImage.saveImageData(targetImageData)
    }
}
print 'Done!'

/**
 * Transform object, recursively transforming all child objects
 *
 * @param pathObject
 * @param transform
 * @return
 */
PathObject transformObject(PathObject pathObject, AffineTransform transform) {
    // Create a new object with the converted ROI
    def roi = pathObject.getROI()
    def roi2 = transformROI(roi, transform)
    def newObject = null
    if (pathObject instanceof PathCellObject) {
        def nucleusROI = pathObject.getNucleusROI()
        if (nucleusROI == null)
            newObject = PathObjects.createCellObject(roi2, pathObject.getPathClass(), pathObject.getMeasurementList())
        else
            newObject = PathObjects.createCellObject(roi2, transformROI(nucleusROI, transform), pathObject.getPathClass(), pathObject.getMeasurementList())
    } else if (pathObject instanceof PathTileObject) {
        newObject = PathObjects.createTileObject(roi2, pathObject.getPathClass(), pathObject.getMeasurementList())
    } else if (pathObject instanceof PathDetectionObject) {
        newObject = PathObjects.createDetectionObject(roi2, pathObject.getPathClass(), pathObject.getMeasurementList())
    } else {
        newObject = PathObjects.createAnnotationObject(roi2, pathObject.getPathClass(), pathObject.getMeasurementList())
    }
    // Handle child objects
    if (pathObject.hasChildren()) {
        newObject.addPathObjects(pathObject.getChildObjects().collect({transformObject(it, transform)}))
    }
    newObject.setName(pathObject.getName())
    return newObject
}

/**
 * Transform ROI (via conversion to Java AWT shape)
 *
 * @param roi
 * @param transform
 * @return
 */
ROI transformROI(ROI roi, AffineTransform transform) {
    def shape = RoiTools.getShape(roi) // Should be able to use roi.getShape() - but there's currently a bug in it for rectangles/ellipses!
    shape2 = transform.createTransformedShape(shape)
    return RoiTools.getShapeROI(shape2, roi.getImagePlane(), 0.5)
}

