import qupath.lib.roi.ROIs
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane

// Parameters - adjust these as needed
def numSquares = 10           // Number of random squares to generate
def squareSizeUM = 100       // Size of each square (in micrometers/um)
def minDistanceUM = 50        // Minimum distance between squares (in um, optional)

// Get the current image data
def imageData = getCurrentImageData()
def server = imageData.getServer()

// Get pixel calibration
def cal = server.getPixelCalibration()
def pixelWidth = cal.getPixelWidthMicrons()
def pixelHeight = cal.getPixelHeightMicrons()

// Check if calibration is available
if (!cal.hasPixelSizeMicrons()) {
    println "Warning: No pixel calibration found. Using pixels instead of micrometers."
    pixelWidth = 1.0
    pixelHeight = 1.0
}

// Convert um to pixels
def squareSize = Math.round(squareSizeUM / pixelWidth) as int
def minDistance = Math.round(minDistanceUM / pixelWidth) as int

println "Square size: ${squareSizeUM} um = ${squareSize} pixels"
println "Pixel size: ${pixelWidth} x ${pixelHeight} um"

// Get image dimensions
def width = server.getWidth()
def height = server.getHeight()

// Create a random number generator
def random = new Random()

// Store created squares to check for overlaps (if needed)
def squares = []

// Generate random squares
for (int i = 0; i < numSquares; i++) {
    def validLocation = false
    def x, y
    def attempts = 0
    def maxAttempts = 1000
    
    // Try to find a valid location
    while (!validLocation && attempts < maxAttempts) {
        // Generate random coordinates
        x = random.nextInt(width - squareSize)
        y = random.nextInt(height - squareSize)
        
        // Check if square is within bounds
        if (x >= 0 && y >= 0 && x + squareSize <= width && y + squareSize <= height) {
            validLocation = true
            
            // Optional: Check minimum distance from other squares
            if (minDistance > 0) {
                for (square in squares) {
                    def dx = Math.abs(square[0] - x)
                    def dy = Math.abs(square[1] - y)
                    if (dx < squareSize + minDistance && dy < squareSize + minDistance) {
                        validLocation = false
                        break
                    }
                }
            }
        }
        attempts++
    }
    
    if (validLocation) {
        // Create square ROI
        def roi = ROIs.createRectangleROI(x, y, squareSize, squareSize, ImagePlane.getDefaultPlane())
        
        // Create annotation with name
        def annotation = PathObjects.createAnnotationObject(roi)
        annotation.setName("R${i + 1}")  // Set name as R1, R2, R3, etc.
        
        // Add to hierarchy
        imageData.getHierarchy().addObject(annotation)
        
        // Store coordinates
        squares.add([x, y])
        
        println "Created annotation R${i + 1} at position (${x}, ${y})"
    } else {
        println "Could not find valid location for square R${i + 1} after ${maxAttempts} attempts"
    }
}

// Resolve hierarchy and update display
imageData.getHierarchy().resolveHierarchy()
fireHierarchyUpdate()

println "Successfully created ${squares.size()} random square annotations (${squareSizeUM} um each)"