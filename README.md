# Introduction
Demonstration to run the AutoCUTS-LM-Analysis code with pictures.
To download the code with some example images, please go to [google drive](https://drive.google.com/drive/folders/1PFH3ZDogbBHG3NMdwgxhIrmVJBk39LKw?usp=sharing)

Cite the code:
[![DOI](https://zenodo.org/badge/313559780.svg)](https://zenodo.org/badge/latestdoi/313559780)

## Requirements 
Before starting, you will need the following:

Install [Image processing toolbox](https://www.mathworks.com/products/image.html) and [Computer vision toolbox](https://www.mathworks.com/products/computer-vision.html) before running the code. 


## Code Guideline

There are four folders in this repository.

1. BiopsyExtraction contains a userfriendly script to random sample biopsies from a block of tissue.  

![01_combineImg](https://user-images.githubusercontent.com/70948370/100009936-cdd99100-2dcf-11eb-893a-8a4dd9c1ef1d.png)



2. Section_segmentation script that is used to detect sections in focus and export them as individual files.

![01_segmentationCombined](https://user-images.githubusercontent.com/70948370/100010000-ea75c900-2dcf-11eb-9ce9-e4a37dcf3f7b.png)



3. AutoCUTSAnalysis contain scripts that show the following steps: 
- Align sections
- Crop aligned stack where only tissue appears
- Using the UNetDense architecture to segment cells, the region of interest could then be defined based on a density map of centroids.
- Filtering small cells away with k-means and reconstruct neurons in 3D. 

![01_alignCombined](https://user-images.githubusercontent.com/70948370/100010106-0d07e200-2dd0-11eb-8599-7741753e5d92.png)

![02_cromCombined](https://user-images.githubusercontent.com/70948370/100010154-20b34880-2dd0-11eb-8523-ddbdca6de5d5.png)

![03_cropArea2](https://user-images.githubusercontent.com/70948370/100010304-5b1ce580-2dd0-11eb-894a-373cf6171c68.png)

![04_reconstruction_combined](https://user-images.githubusercontent.com/70948370/100068174-75d67500-2e37-11eb-89a6-39aa7ab141be.png)



4. ValidationDeeplearning shows the pixelwise and objectwise validation of the UNetDense model output.

![01_pixelValidation](https://user-images.githubusercontent.com/70948370/100068549-fe551580-2e37-11eb-98ab-5f6a0f963764.png)

![01_objectviseValidation](https://user-images.githubusercontent.com/70948370/100068564-01500600-2e38-11eb-829a-c3b91dce2898.png)


## Contact
If you have any questions or suggestions, you can reach me via e-mail at nick.yin.larsen@gmail.com

#### Author: Nick Yin Larsen

