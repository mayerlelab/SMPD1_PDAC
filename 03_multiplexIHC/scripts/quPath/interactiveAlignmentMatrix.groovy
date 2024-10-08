/**
 * ******************************************************************************************* 
 * *********         Alignemnt of all Images using interactive alignemt             **********
 * *******************************************************************************************
 *
 * Use this script if you want to pull a transformation matrix from the Interactive Image Alignment GUI
 *
 * 1. Open reference image in viewer
 * 2. Open the Interactive Image Alignment overlay, align an image
 * 3. While the overlay is still open, set 'name' to the name of the current moving image, and run script
 *
 * @author Ujjwal Mukund Mahajan (modified from @author Mike Nelson and Mark Zaidi)
 * (https://petebankhead.github.io/qupath/scripting/2018/03/14/script-export-labelled-images.html)
 * (https://github.com/MarkZaidi/QuPath-scripts)
 */
 
 
def name='TMA1_CK19.mrxs' 

path = buildFilePath(PROJECT_BASE_DIR, 'Affine')
mkdirs(path)
path = buildFilePath(PROJECT_BASE_DIR, 'Affine', name)



import qupath.lib.gui.align.ImageServerOverlay
// import qupath.ext.align.gui.ImageServerOverlay -> for QuPAth version 3.0.0

def overlay = getCurrentViewer().getCustomOverlayLayers().find {it instanceof ImageServerOverlay}

affine = overlay.getAffine()

print affine
afString = affine.toString()
afString = afString.minus('Affine [').minus(']').trim().split('\n')
cleanAffine =[]
afString.each{
    temp = it.split(',')
    temp.each{cleanAffine << Double.parseDouble(it)}
}

def matrix = []
affineList = [0,1,3,4,5,7]
for (i=0;i<12; i++){
if (affineList.contains(i))
    matrix << cleanAffine[i]
}

new File(path).withObjectOutputStream {
    it.writeObject(matrix)
}
print 'Done!'
