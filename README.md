# Introduction
Demonstration to run the AutoCUTS-LM-Analysis code with pictures.
The full setup with scripts and some image examples will be available after publication.

Cite the code: https://github.com/Nick7900/AutoCUTS-LM-Analysis or article after publication.
## Requirements 
Before starting, you will need the following:

- Install [MATLAB](https://www.mathworks.com/downloads/).
- Install [Image processing toolbox](https://www.mathworks.com/products/image.html) and [Computer vision toolbox](https://www.mathworks.com/products/computer-vision.html) before running the code. 
- Install [R](https://cran.r-project.org/).

## Code Guideline

There are five folders in this repository.

1. BiopsyExtraction contains a userfriendly script to random sample biopsies from a block of tissue.  

![combineImg_demoUpload](https://user-images.githubusercontent.com/70948370/100726005-ab91d580-33c4-11eb-8aa9-3a10eab9d7fc.PNG)



2. Section_segmentation script was used to detect focused sections and export them as individual files. The script was appled on a test image.

![01_segmentationCombined](https://user-images.githubusercontent.com/70948370/100010000-ea75c900-2dcf-11eb-9ce9-e4a37dcf3f7b.png)



3. AutoCUTSAnalysis includes scripts implemented on a test collection of 50 images that involve the following steps:
- Align sections
- Crop aligned stack where only tissue appears
- Using the UNetDense architecture to segment cells, the region of interest could then be defined based on a density map of centroids.
- Filtering small cells away with k-means and reconstruct neurons in 3D. 

![01_alignCombined](https://user-images.githubusercontent.com/70948370/100010106-0d07e200-2dd0-11eb-8599-7741753e5d92.png)

![02_cromCombined](https://user-images.githubusercontent.com/70948370/100010154-20b34880-2dd0-11eb-8523-ddbdca6de5d5.png)

![03_cropArea2](https://user-images.githubusercontent.com/70948370/100010304-5b1ce580-2dd0-11eb-894a-373cf6171c68.png)

![04_reconstruction_combined](https://user-images.githubusercontent.com/70948370/100068174-75d67500-2e37-11eb-89a6-39aa7ab141be.png)



4. ValidationDeeplearning will estimate the pixelwise and objectwise validation of the UNetDense model output.
- Pictures of objectwise validation will be attach after publication.
- Look at False-Negative plot in the article for further information.
![01_objectviseValidation_2](https://user-images.githubusercontent.com/70948370/100723354-f3fbc400-33c1-11eb-8324-dc380dd80791.PNG)

5. CylindricalKfunction
Demonstration to run the cylindrical K-function code with pictures. The function was introduced in '[The cylindrical K-function and Poisson line cluster point processe](https://arxiv.org/abs/1503.07423)' by Jesper Møller et al. 
All credits go to the script writer, Ninna Vihrs.
- The code analyses the columnar structure in 3D point patterns for one test subject.
- Pictures will be attach after publication.



## Contact
If you have any questions or suggestions, you can reach Nick via e-mail at nick.yin.larsen@gmail.com for the MATLAB scripts or Ninna via e-mail at nvihrs@math.aau.dk for the R-script.

#### Author: Nick Yin Larsen

