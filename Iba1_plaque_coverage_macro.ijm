//definition of variables
dir = getDirectory("image")
name=getTitle()
path = dir + name
// splitting images
name = getTitle();
run("Split Channels");
selectWindow("C2-"+name)
close()
selectWindow("C3-"+name)
close()
selectWindow("C4-"+name)
close()
// define ROIs 
selectWindow("C1-"+name);
run("Set Measurements...", "area mean standard min perimeter shape feret's display redirect=None decimal=3");
raw = getTitle();
setTool("polygon");
waitForUser("Define Brain, cortex, hippo, fimbria, alveus, thalamus ROI")
roiManager("Select", 0);
roiManager("Rename", "Brain");
roiManager("Select", 1);
roiManager("Rename", "Cortex");
roiManager("Select", 2);
roiManager("Rename", "Hippocampus");
roiManager("Select", 3);
roiManager("Rename", "Fimbria");
roiManager("Select", 4);
roiManager("Rename", "Alveus");
roiManager("Select", 5);
roiManager("Rename", "Thalamus");
roiManager("Select", newArray(0,1,2,3,4,5));
roiManager("Measure");
selectWindow("Results");
// thresholding
roiManager("Select", 0);
run("Threshold...");
waitForUser("set threshold method");
run("Convert to Mask");
run("Create Selection");
roiManager("Add");
roiManager("Select", 6);
roiManager("Rename", "Threshold Brain");
roiManager("Measure");
// duplication of images
selectWindow(raw)
roiManager("Show All");
roiManager("Show None");
selectWindow(raw)
run("Duplicate...", " ");
rename("cortex_duplicate1")
run("Duplicate...", " ");
rename("hippo_duplicate2")
run("Duplicate...", " ");
rename("fimbria_duplicate3");
run("Duplicate...", " ");
rename("alveus_duplicate4");
run("Duplicate...", " ");
rename("thalamus_duplicate5");
// cortex roi
selectWindow("cortex_duplicate1")
roiManager("Select", 1);
run("Clear Outside");
roiManager("Show All");
roiManager("Show None");
run("Create Selection");
roiManager("Add");
roiManager("Select", 7);
roiManager("Rename", "Threshold Cortex");
roiManager("Measure");
// hippo roi
selectWindow("hippo_duplicate2");
roiManager("Select", 2);
run("Clear Outside");
roiManager("Show All");
roiManager("Show None");
run("Create Selection");
roiManager("Add");
roiManager("Select", 8);
roiManager("Rename", "Threshold Hippo");
roiManager("Measure");
// fimbria roi
selectWindow("fimbria_duplicate3");
roiManager("Select", 3);
run("Clear Outside");
roiManager("Show All");
roiManager("Show None");
run("Create Selection");
roiManager("Add");
roiManager("Select", 9);
roiManager("Rename", "Threshold Fimbria");
roiManager("Measure");
// alveus roi
selectWindow("alveus_duplicate4");
roiManager("Select", 4);
run("Clear Outside");
roiManager("Show All");
roiManager("Show None");
run("Create Selection");
roiManager("Add");
roiManager("Select", 10);
roiManager("Rename", "Threshold Alveus");
roiManager("Measure");
// thalamus roi
selectWindow("thalamus_duplicate5");
roiManager("Select", 5);
run("Clear Outside");
roiManager("Show All");
roiManager("Show None");
run("Create Selection");
roiManager("Add");
roiManager("Select", 11);
roiManager("Rename", "Threshold Thalamus");
roiManager("Measure");
// particle analyzer brain
selectWindow(raw)
run("Analyze Particles...", "size=20-2000 circularity=0.06-1.00 show=[Overlay Masks] display exclude include summarize add in_situ");
// particle analyzer cortex
selectWindow("cortex_duplicate1");
roiManager("Show All");
roiManager("Show None");
selectWindow("cortex_duplicate1");
run("Analyze Particles...", "size=20-2000 circularity=0.06-1.00 show=[Overlay Masks] display exclude include summarize add in_situ");
// particle analyzer hippo
selectWindow("hippo_duplicate2");
roiManager("Show All");
roiManager("Show None");
selectWindow("hippo_duplicate2");
run("Analyze Particles...", "size=20-2000 circularity=0.06-1.00 show=[Overlay Masks] display exclude include summarize add in_situ");
// particle analyzer fim
selectWindow("fimbria_duplicate3");
roiManager("Show All");
roiManager("Show None");
selectWindow("fimbria_duplicate3");
run("Analyze Particles...", "size=20-2000 circularity=0.06-1.00 show=[Overlay Masks] display exclude include summarize add in_situ");
// particle analyzer alveus
selectWindow("alveus_duplicate4");
roiManager("Show All");
roiManager("Show None");
selectWindow("alveus_duplicate4");
run("Analyze Particles...", "size=20-2000 circularity=0.06-1.00 show=[Overlay Masks] display exclude include summarize add in_situ");
// particle analyzer thalamus
selectWindow("thalamus_duplicate5");
roiManager("Show All");
roiManager("Show None");
selectWindow("thalamus_duplicate5");
run("Analyze Particles...", "size=20-2000 circularity=0.06-1.00 show=[Overlay Masks] display exclude include summarize add in_situ");
// safe ROIs, results and summary
selectWindow("Results")
saveAs("Results", path + "_results.csv")
selectWindow("Summary")
saveAs("Results", path + "_summary.csv")
selectWindow("ROI Manager")
roiManager("Save",  path + "_ROIset.zip")
// close everything
run("Close All")
close("Results")
close("ROI Manager")
close(name + "_summary.csv")
waitForUser("All done")

