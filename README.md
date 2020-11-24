# AutoCTUS-LM
Demostration to run the AutoCUTS-pipeline code with pictures.
To download the code with some example images, please go to [google drive](https://pages.github.com/)
Cite the code:
[![DOI](https://zenodo.org/badge/313559780.svg)](https://zenodo.org/badge/latestdoi/313559780)


## Code Guideline

There are four folders in this repository.

1. BiopsyExtraction contains a script to random sample the  

2. Section_segmentation scripts that train the three models described in the paper. There are a few data samples that the code can run on them smoothly.

3. AutoCUTSAnalysis  shows how we run the trained models during testing. The three types of trained models can be downloaded in the google drive. Please put them under "./demo_step003_training/models" to run them.

4. ValidationDeeplearning shows how we run the trained models during testing. The three types of trained models can be downloaded in the google drive. Please put them under "./demo_step003_training/models" to run them

![01_combineImg](https://user-images.githubusercontent.com/70948370/100009936-cdd99100-2dcf-11eb-893a-8a4dd9c1ef1d.png)

![01_segmentationCombined](https://user-images.githubusercontent.com/70948370/100010000-ea75c900-2dcf-11eb-9ce9-e4a37dcf3f7b.png)

![01_alignCombined](https://user-images.githubusercontent.com/70948370/100010106-0d07e200-2dd0-11eb-8599-7741753e5d92.png)

![02_cromCombined](https://user-images.githubusercontent.com/70948370/100010154-20b34880-2dd0-11eb-8523-ddbdca6de5d5.png)

![03_cropArea2](https://user-images.githubusercontent.com/70948370/100010304-5b1ce580-2dd0-11eb-894a-373cf6171c68.png)

![04_reconstruction_combined](https://user-images.githubusercontent.com/70948370/100068174-75d67500-2e37-11eb-89a6-39aa7ab141be.png)

![01_pixelValidation](https://user-images.githubusercontent.com/70948370/100068549-fe551580-2e37-11eb-98ab-5f6a0f963764.png)

![01_objectviseValidation](https://user-images.githubusercontent.com/70948370/100068564-01500600-2e38-11eb-829a-c3b91dce2898.png)



####
Author: Nick Yin Larsen

