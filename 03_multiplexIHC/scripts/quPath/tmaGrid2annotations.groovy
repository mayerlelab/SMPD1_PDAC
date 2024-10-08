/**
 * ******************************************************************************************* 
 * *********       Script to export TMA grid & convert it to annotations            **********
 * *******************************************************************************************
 *
 * @author Ujjwal Mukund Mahajan (modified from @ Pete Bankhead) (https://forum.image.sc/t/transform-and-copy-for-tma-grid/25024/4)
 */
 
 
// Grid to annotations
def annotations = getTMACoreList().collect {
    def annotation = new qupath.lib.objects.PathAnnotationObject(it.getROI())
    annotation.setName(it.getUniqueID())
    return annotation
}
clearAllObjects()
addObjects(annotations)
fireHierarchyUpdate()
