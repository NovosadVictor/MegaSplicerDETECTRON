# SRPseq

Splicing Regulation Prediction from RNA-seq data. Please cite this paper if you are using SRPseq in your work:

Nersisyan S, Novosad V, Tonevitsky A. SRPseq: Splicing Regulation Prediction from RNA-seq data. bioRxiv. 2021 Aug 04. doi: [10.1101/2021.08.03.454798](https://doi.org/10.1101/2021.08.03.454798).


## Table of Contents  
[Introduction](#introduction)  
[Installation](#installation)  
[Tutorial](#tutorial)  
[Running SRPseq](#running-srpseq)  
[Functions and classes](#functions-and-classes)  

# Introduction
<img align="right" width="400px" src="https:https://github.com/NovosadVictor/SRPseq/static/flowchart.png?raw=true">
<div>
<p>In general words, the core idea of our approach consists in the sequentiality of the splicing process. Specifically, by focusing on a single pre-RNA of a certain gene, our aim is to model the splicing process as a set of specific steps where each step is dependent on the previous steps and the “environment” of the current step. In particular, we assign each step to the spliceosome decision on whether to include or to exclude a considered exon. Importantly, each step has its probability which depends on the previously included/excluded exons and the expression levels of splicing factors related to the considered exon. The latter means that it is believed that if the splicing factor is an enhancer for a specific exon in a specific tissue, then the higher its expression the higher the probability of the exon inclusion, and otherwise for the silencing splicing factors. Briefly, a pipeline is implemented as follows:</p>
<ol>
  <li><i>Construction of isoform-exon tree:</i> Divide groups of isoforms by sequence of common exons.</li>
  <li><i>List of RBPs with matching motifs:</i> For each node of the tree set a list of RBPs which motifs are present on the node's exon.</li>
  <li><i>Model training:</i>Train linear model based on the specified RNA-seq data for each node of the tree.</li>
  <li><i>Evaluation:</i> evaluate each model and make summary of all steps.</li>
</ol>

Input data can consist from different tissue-samples and datasets, and each dataset should be labeled by one of the following types:
<ol>
<li><i>Training set:</i> samples from training datasets will be used for tuning the regression models. At least one such dataset is required.</li>
<li><i>Validation (test) set:</i> performance of models evaluated on the validation sets. At least one such dataset is required.</li>
</ol>

If not specified, the whole data randomly split into training and validation datasets with 3 to 1 ratio.

</div>

# Installation

### Prerequisites:
Make sure you have installed all of the following prerequisites on your development machine:
  - python3.6+  
  - pip3


### SRPseq installation:  
`pip3 install srpseq`

# Tutorial


# Running SRPseq

## Step 1: data preparation

Before running the tool, you should prepare optional two tsv tables containing expression data of RNA-binding proteins (or all genes expression levels), expression levels of considered gene isoforms and a table with RNA-binding proteins motifs. RNA-seq tables should contain log2 transformed FPKM values associated with samples (rows) and genes / isoforms (columns):

<details>
  <summary>Example</summary>
  
  |            | RBP_1 | RBP_2 |
  | ---------- | --------- | --------- |
  | Sample 1   | 17.17     | 365.1     |
  | Sample 2   | 56.99     | 123.9     |
  | ...        |           |           |
  | Sample 98  | 22.22     | 123.4     |
  | Sample 99  | 23.23     | 567.8     |
  | ...        |           |           |
  | Sample 511 | 10.82     | 665.8     |
  | Sample 512 | 11.11     | 200.2     |
</details>


Annotation table format is different for classification and survival analysis. For classification it should contain binary indicator of sample class (e.g., 1 for recurrent tumor and 0 for non-recurrent), dataset (batch) label and dataset type (Training/Filtration/Validation).  
It is important that `Class = 1` represents "Positives" and `Class = 0` are "negatives", otherwise accuracy scores may be calculated incorrectly.   
Note that annotation should be present for each sample listed in the data table in the same order:

<details>
  <summary>Example</summary>
  
  |            | Class | Dataset  | Dataset type |
  | ---------- | ----- | -------- | ------------ |
  | Sample 1   | 1     | GSE3494  | Training     |
  | Sample 2   | 0     | GSE3494  | Training     |
  | ...        |       |          |              |
  | Sample 98  | 0     | GSE12093 | Filtration   |
  | Sample 99  | 0     | GSE12093 | Filtration   |
  | ...        |       |          |              |
  | Sample 511 | 1     | GSE1456  | Validation   |
  | Sample 512 | 1     | GSE1456  | Validation   |
</details>


For survival analysis, annotation table should contain binary event indicator and time to event:
<details>
  <summary>Example</summary>
  
  |            | Event | Time to event | Dataset  | Dataset type |
  | ---------- | ----- | ------------- | -------- | ------------ |
  | Sample 1   | 1     | 100.1         | GSE3494  | Training     |
  | Sample 2   | 0     | 500.2         | GSE3494  | Training     |
  | ...        |       |               |          |              |
  | Sample 98  | 0     | 623.9         | GSE12093 | Filtration   |
  | Sample 99  | 0     | 717.1         | GSE12093 | Filtration   |
  | ...        |       |               |          |              |
  | Sample 511 | 1     | 40.5          | GSE1456  | Validation   |
  | Sample 512 | 1     | 66.7          | GSE1456  | Validation   |
</details>


Table with *n* / *k* grid for exhaustive feature selection:  
*n* is a number of selected features, *k* is a length of each features subset.  

If you are not sure what values for *n* *k* to use, see [Step 3: defining a *n*, *k* grid](#step-3-defining-a-n-k-grid)  

<details>
  <summary>Example</summary> 
   
  | n   | k   |  
  | --- | --- |  
  | 100 | 1   |  
  | 100 | 2   |  
  | ... | ... |  
  | 20  | 5   |  
  | 20  | 10  |  
  | 20  | 15  |  
</details>


## Step 2: creating configuration file

Configuration file is a json file containing all customizable parameters for the model (classification and survival analysis)  

<details>
  <summary>Available parameters</summary> 

  🔴!NOTE! - All paths to files / directories can be either relative to the configuration file directory or absolute paths 
  * `data_path`
      Path to csv table of the data.

  * `annotation_path`
      Path to csv table of the data annotation.

  * `n_k_path`
      Path to a *n*/*k* grid file.

  * `output_dir`
      Path to directory for output files. If it doesn't exist, it will be created.

  * `feature_pre_selector`  
      Name of feature pre-selection function from [feature pre-selectors section](#functions-and-classes).

  * `feature_pre_selector_kwargs`  
      Object/Dictionary of keyword arguments for feature pre-selector function.

  * `feature_selector`  
      Name of feature selection function from [feature selectors section](#functions-and-classes).

  * `feature_selector_kwargs`  
      Object/Dictionary of keyword arguments for feature selector function. Boolean `use_filtration` indicates whether to use *Filtration* dataset besides *Training* dataset for the selector function.

  * `preprocessor`
      Name of class for data preprocessing from [sklearn.preprocessing](https://scikit-learn.org/stable/modules/preprocessing.html).

  * `preprocessor_kwargs`
      Object/Dictionary of keyword arguments for preprocessor class initialization.  
      If you are using `sklearn` model, use `kwargs` parameters from the documentation of the model.

  * `model`  
      Name of class for classification / survival analysis from [Classifiers / Regressors section](#functions-and-classes).

  * `model_kwargs`
      Object/Dictionary of keyword arguments for model initialization.  
      If you are using `sklearn` model, use `kwargs` parameters from the documentation of the model.

  * `model_CV_ranges`
      Object/Dictionary defining model parameters which should be cross-validated. Keys are parameter names, values are lists for grid search.

  * `model_CV_folds`
      Number of folds for K-Folds cross-validation.

  * `scoring_functions`
      List with names for scoring functions (from [Accuracy scores section](#functions-and-classes)) which will be calculated for each model. If you need to pass parameters to the function (e.g. `year` in dynamic auc score), you can use object {"name": `function name`, "kwargs": `parameters object`}.

  * `main_scoring_function`
      Key from scoring_functions dict defining the "main" scoring function which will be optimized during cross-validation and will be used for model filtering.

  * `main_scoring_threshold`
      A number defining threshold for model filtering: models with score below this threshold on training/filtration sets will not be further evaluated.

  * `n_processes`
      Number of processes / threads to run on.
  
  * `random_state`
      Random seed (set to an arbitrary integer for reproducibility).

  * `verbose`
      If *true*, print running time for each pair of *n*, *k*.
</details>

## Step 3: running the pipeline

## Step 4: generating report


# Functions and classes


