Pipeline to analysed the individual images from 2_SectionSegmentation.

1. First, copy the files from the folder 'Example_save' and place it in the 'AutuCUTS_Pipeline'.
Rename the files in the folder to, e.g. 1 and all files will be chronologically ordered from 1 (1)...1 (N).
The folder should thereafter be renamed to the desired name and add '_1_order' e.g. 'Example_1_order'.

2. Run 'aAlignAuto' and images will be aligned.

3. Run 'bCropToDeepLearning' to crop areas that only include tissue

4. After the images from 'bCropToDeepLearning' have been segmented with the UNetDense architecture, then run 'cCropLayer3' region of interest for analysis.

5. Run 'dPredicSegmentation' to do the 3D-reconstruction of neurons and get morphological measurements out.

6. Run 'eHistogramDataPlot' to read the output from 'dPredicSegmentationFinal' to estimate mean, standard deviation, coefficient of variation and plot histograms.
