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