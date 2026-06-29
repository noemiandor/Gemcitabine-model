input = "C:/Users/80027908/Documents/";
output = input + File.separator + "F6_1_corr_output";
File.makeDirectory(output);

list = newArray("176", "154", "186", "181", "1082", "2147", "114", "2187", "184", "15");

setBatchMode(true);
for (i = 0; i < list.length; i++){
	makeChange(input, output, list[i]);
}

function makeChange(input, output, track) {
	File.openSequence(input, " filter=(corr_F6_1_" + track + "_)");
	rename("brightfield");
	
	File.openSequence(input, " filter=(corr_" + track + "_)");
	rename("plots");
		
	run("Combine...", "stack1=brightfield stack2=plots");
	
	saveAs("tiff", output + File.separator + track);
	print(track);
	close("*");
}


// inter-division

input = "C:/Users/80027908/Documents/";
output = input + File.separator + "F6_1_corr_div_output";
File.makeDirectory(output);

list = newArray("322", "940", "2910", "2051", "1128", "599", "3336", "1306", "539", "1336");

setBatchMode(true);
for (i = 0; i < list.length; i++){
	makeChange(input, output, list[i]);
}

function makeChange(input, output, track) {
	File.openSequence(input, " filter=(corr_div_F6_1_" + track + "_)");
	rename("brightfield");
	
	File.openSequence(input, " filter=(corr_div_" + track + "_)");
	rename("plots");
		
	run("Combine...", "stack1=brightfield stack2=plots");
	
	saveAs("tiff", output + File.separator + track);
	print(track);
	close("*");
}