input = "C:/Users/80027908/Desktop/K00_GemcitabineExposure_033023/F6_1/HALO Markup";
output = input + File.separator + "output";
File.makeDirectory(output);

list = getFileList(input);

setBatchMode(true);
for (i = 0; i < list.length; i++){
	saveSelected(input, output, list[i]);
}

function saveSelected(input, output, filename) {
	open(input + File.separator + filename);
	filename_pure = File.nameWithoutExtension;
	run("RGB Color", "slices");
	run("Stack to Images");
	selectImage(filename_pure + "-0001");
	close("\\Others");
	run("Scale...", "x=.2854 y=.2854 width=1468 height=1100 interpolation=Bilinear average create");
    saveAs("tiff", output + File.separator + filename);
    print(filename);
	close("*");
}


// Overlay
input_1 = "C:/Users/80027908/Desktop/K00_GemcitabineExposure_033023/F6_1/HALO Markup/output";
input_2 = "C:/Users/80027908/Desktop/K00_GemcitabineExposure_033023/F6_1";
output = input_1 + File.separator + "overlay_output";
File.makeDirectory(output);

list_1 = getFileList(input_1);
list_2 = getFileList(input_2);

setBatchMode(true);
for (i = 0; i < list_1.length; i++){
	makeOverlay(input_1, input_2, output, list_1[i], list_2[i]);
}

function makeOverlay(input_1, input_2, output, HALO_name, bright_name) {
	open(input_1 + File.separator + HALO_name);
	rename("HALO");
	open(input_2 + File.separator + bright_name);
	rename("brightfield");
	run("Add Image...", "image=HALO x=0 y=0 opacity=60 zero");
	run("Flatten");
    saveAs("tiff", output + File.separator + bright_name);
    print(bright_name);
	close("*");
}












