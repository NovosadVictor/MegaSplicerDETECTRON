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

## Step 2: creating configuration file

## Step 3: running the pipeline

## Step 4: generating report


# Functions and classes


