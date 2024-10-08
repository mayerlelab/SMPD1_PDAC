/*
 * author: Ujjwal Mukund Mahajan
 *
 * **************************************************
 * **                 image montage                **
 * **************************************************
 *
 */

#@File(label = "Input directory", style = "directory") input
#@File(label = "Output directory", style = "directory") output
#@String(label = "File suffix", value = ".tif") suffix

var outputDir
var inputDir

outputDir = File.getName(output);
inputDir = File.getName(input);

setBatchMode(true);
processFolder(input);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix) && (indexOf(list[i], "_R1_") != -1))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
           run("Image Sequence...", "open=["+input + File.separator + list[i]+"] sort use");
					 print("----------------------------------------------");
					 print("Processing: ",getTitle());
					 print("----------------------------------------------");
           run("Subtract Background...", "rolling=150 disable stack");
           //run("Kuwahara Filter", "sampling=5 filter stack");
           run("Remove Outliers...", "radius=5 threshold=50 which=Bright stack");
           run("Gamma...", "value=1.30 stack");
					 stack = getTitle();
           num_channels = nSlices;
           // Get each slice label.
           labels = newArray(num_channels);
           for (i = 1; i <= num_channels; i++) {
                setSlice(i);
                labels[i-1] = getInfo("slice.label");
            }
           // Extract channel from label by removing common strings.
           channels = newArray(num_channels);
           for (j = 0; j < num_channels; j++) {
           	    label = labels[j];
           	    channels[j] = substring(label, lastIndexOf(label, "_"));
           	    channels[j] = replace(channels[j], "_", "");
           }
           // Combine into final stack.
           run("Stack to Images");
           run("Images to Stack", "name=Stack title=[] use");
           new_stack = getTitle();
           // Label with channel names.
           setSlice(1);
           run("Set Label...", "label=Composite");
           for (i = 0; i < num_channels; i++) {
           	     setSlice(i+1);
           	     label = channels[i];
           	     run("Set Label...", "label="+label+"");
           	}
           run("Make Montage...",  "columns=3 rows=5 scale=0.5 border=10 font=100 label use");
           saveAs("JPEG", output + File.separator + stack);
           run("Close All");
}
