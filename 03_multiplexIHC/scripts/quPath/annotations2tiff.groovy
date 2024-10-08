/**
 * *******************************************************************************************
 * *********        Script to export pixels & annotations for TMA images             *********
 * *******************************************************************************************
 *
 * The downsample value and coordinates are encoded in each image file name.
 *
 * The annotations are exported as 8-bit labelled images.
 * These labels depend upon annotation classifications; a text file giving the key is written for reference.
 *
 * The labelled image can also optionally use indexed colors to depict the colors of the
 * original classifications within QuPath for easier visualization & comparison.
 *
 * @author Ujjwal Mukund Mahajan (modified from @author Pete Bankhead)
 * (https://petebankhead.github.io/qupath/scripting/2018/03/14/script-export-labelled-images.html)
 */

import qupath.lib.common.ColorTools
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.RegionRequest
import qupath.lib.roi.RoiTools
import qupath.lib.gui.scripting.QPEx

import javax.imageio.ImageIO
import java.awt.Color
import java.awt.image.BufferedImage
import java.awt.image.DataBufferByte
import java.awt.image.IndexColorModel
import qupath.lib.images.servers.LabeledImageServer

// Requested pixel size - used to define output resolution
// Set <= 0 to use the full resolution (whatever that may be)
double requestedPixelSizeMicrons = 0

// Maximum size of an image tile when exporting
int maxTileSize = 1000

// Export the original pixels (assumed ot be RGB) for each tile
boolean exportOriginalPixels = true

// Export a labelled image for each tile containing annotations
boolean exportAnnotationLabelledImage = true

// Don't export any cores marked as 'missing'
boolean skipMissingCores = true

// NOTE: The following parameters only matter if exportAnnotationLabelledImage == true
// Ignore annotations that don't have a classification set
boolean skipUnclassifiedAnnotations = true
// Skip tiles without annotations (only applies to label exports - all image tiles will be written)
boolean skipUnannotatedTiles = true
// Create an 8-bit indexed image
// This is very useful for display/previewing - although need to be careful when importing into other software,
// which can prefer to replaces labels with the RGB colors they refer to
boolean createIndexedImageLabels = true

// NOTE: The following parameter only matters if exportOriginalPixels == true
// Define the format for the export image
def imageFormat = 'TIFF'

// Output directory for storing the tiles
def pathOutput = QPEx.buildFilePath(QPEx.PROJECT_BASE_DIR, 'raw_cores')
QPEx.mkdirs(pathOutput)

//------------------------------------------------------------------------------------

// Get the main QuPath data structures
def imageData = QPEx.getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()

// Get the annotations that have ROIs & are have classifications (if required)
def annotations = hierarchy.getFlattenedObjectList(null).findAll {
    it.isAnnotation() && it.hasROI() && (!skipUnclassifiedAnnotations || it.getPathClass() != null) }

// Get all the represented classifications
def pathClasses = annotations.collect({it.getPathClass()}) as Set

// We can't handle more than 255 classes (because of 8-bit representation)
if (pathClasses.size() > 255) {
    print 'Sorry! Cannot handle > 255 classications - number here is ' + pathClasses.size()
    return
}

// For each classification, create an integer label >= 1 & cache a color for drawing it
// Also, create a LUT for visualizing more easily
def labelKey = new StringBuilder()
def pathClassColors = [:]
int n = pathClasses.size() + 1
def r = ([(byte)0] * n) as byte[]
def g = r.clone()
def b = r.clone()
def a = r.clone()
pathClasses.eachWithIndex{ PathClass entry, int i ->
    // Record integer label for key
    int label = i+1
    String name = entry == null ? 'None' : entry.toString()
    labelKey << name << '\t' << label << System.lineSeparator()
    // Store Color based on label - needed for painting with Graphics2D
    pathClassColors.put(entry, new Color(label, label, label))
    // Update RGB values for IndexColorModel
    // Use gray as the default color indicating no classification
    int rgb = entry == null ? ColorTools.makeRGB(127, 127, 127) : entry.getColor()
    r[label] = ColorTools.red(rgb)
    g[label] = ColorTools.green(rgb)
    b[label] = ColorTools.blue(rgb)
    a[label] = 255
}

// Check if we've anything to do
if (!exportAnnotationLabelledImage && !exportOriginalPixels) {
    print 'Nothing to export!'
    return
}

// Calculate the downsample value
double downsample = 1
if (requestedPixelSizeMicrons > 0)
    downsample = requestedPixelSizeMicrons / server.getPixelCalibration().getAveragedPixelSizeMicrons()
    
// Calculate the tile spacing in full resolution pixels
int spacing = (int)(maxTileSize * downsample)

// Get the cores
def cores = hierarchy.getAnnotationObjects()


// Write the label 'key' // ## not necessary
// println labelKey
// def keyName = String.format('%s_(downsample=%.3f,tiles=%d)-key.txt', getProjectEntry().getImageName(), downsample, maxTileSize)
// def fileLabels = new File(pathOutput, keyName) 
// fileLabels.text = labelKey.toString()

def stain_name = getProjectEntry().getImageName()
String stain_name1 = stain_name.substring(stain_name.indexOf("_")+1)
String stain_name2 =  stain_name1.substring(0,  stain_name1.indexOf("."));



// Handle the cores in parallel
cores.parallelStream().forEach { core ->

    def request = RegionRequest.createInstance(server.getPath(), downsample, core.getROI())
    
    // Create a suitable base image name
    String name = String.format('%s_%s',
            core.getName(),
            stain_name2.trim()
    )

    // Export the raw image pixels if necessary
    // If we do this, store the width & height - to make sure we have an exact match
    int width = -1
    int height = -1
    if (exportOriginalPixels) {
        def img = server.readBufferedImage(request)
        width = img.getWidth()
        height = img.getHeight()
        def fileOutput = new File(pathOutput, name + '.' + imageFormat.toLowerCase())
        ImageIO.write(img, imageFormat, fileOutput)
    }

    // Export the labelled tiles if necessary
    if (exportAnnotationLabelledImage) {
        // Calculate dimensions if we don't know them already
        if (width < 0 || height < 0) {
            width = Math.round(request.getWidth() / downsample)
            height = Math.round(request.getHeight() / downsample)
        }
        // Fill the annotations with the appropriate label
        def imgMask = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY)
        def g2d = imgMask.createGraphics()
        g2d.setClip(0, 0, width, height)
        g2d.scale(1.0/downsample, 1.0/downsample)
        g2d.translate(-request.getX(), -request.getY())
        int count = 0
        for (annotation in annotations) {
            def roi = annotation.getROI()
            if (!request.intersects(roi.getBoundsX(), roi.getBoundsY(), roi.getBoundsWidth(), roi.getBoundsHeight()))
                continue
            def shape = PathROIToolsAwt.getShape(roi)
            def color = pathClassColors.get(annotation.getPathClass())
            g2d.setColor(color)
            g2d.fill(shape)
            count++
        }
        g2d.dispose()
        if (count > 0 || !skipUnannotatedTiles) {
            // Extract the bytes from the image
            def buf = imgMask.getRaster().getDataBuffer() as DataBufferByte
            def bytes = buf.getData()
            // Check if we actually have any non-zero pixels, if necessary -
            // we might not if the annotation bounding box intersected the region, but the annotation itself does not
            if (skipUnannotatedTiles && !bytes.any { it != (byte)0 })
                return
            // If we want an indexed color image, create one using the data buffer & LUT
            if (createIndexedImageLabels) {
                def colorModel = new IndexColorModel(8, n, r, g, b, a)
                def imgMaskColor = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_INDEXED, colorModel)
                System.arraycopy(bytes, 0, imgMaskColor.getRaster().getDataBuffer().getData(), 0, width*height)
                imgMask = imgMaskColor
            }
            // Write the mask
            def fileOutput = new File(pathOutput, name + '-labels.png')
            ImageIO.write(imgMask, 'PNG', fileOutput)
        }
    }
}
print 'Done!'
