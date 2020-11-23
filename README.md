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

![03_cropArea](https://user-images.githubusercontent.com/70948370/100010204-3294eb80-2dd0-11eb-9347-7f7d134086b0.png)

![04_reconstruction](https://user-images.githubusercontent.com/70948370/100010227-39bbf980-2dd0-11eb-96d4-b26f9f004f3c.png)


####
Author: Nick Y. Larsen

