# COLLAGENE: Toolbase for Building Collaborative Biomedical Data Analysis Pipelines

This repository contains the release version of the COLLAGENE library that provides tools for building collaborative and federated biomedical data analysis with a focus on data security and individual-level data privacy.

## What is COLLAGENE and how do I use it?
COLLAGENE is a collection of tools that makes it easy to build collaborative analysis pipelines. 

For extensive discussion about Motivation and Overview, please see *Docs/*.

### Installation
Please follow the instructions under *installation/* for building the executables and setting up the environment and running COLLAGENE through docker. 

## Main Components of Analysis by COLLAGENE
COLLAGENE provides number of components that can be used for building collaborative data analysis pipelines. Most functionalities are implemented in the wrapper named COLLAGENE.sh under *scripts/* directory.

### 1) Key Generation, encryption, collective decryption, ciphertext vitals, and masking of encrypted data:
KeyMaker is a cloud-based service that is used for generating and sharing secret keys among collaborating sites. It can be accessed [here](https://www.secureomics.org/KeyMaker). This step simplifies key-generation step. As KeyMaker is central for generation of the keys, sites must trust the KeyMaker service. It should be noted that KeyMaker does not take part in data processing steps and the keys can be used in local settings without any more interaction with KeyMaker after they are generated.

COLLAGENE provides the tools to encrypt and collectively decrypt data matrices (e.g., genotypes, phenotypes, covariates), modules to keep track of the ciphertext information and security level (128-bits) that are important to ensure that they can be operated on.

Usage examples of these tools are included under *examples/* directory.

### 2) Networking and File I/O:
COLLAGENE's current approach is based on exchanging data among sites using file transfers. The networking and file I/O simplifies online computation requirements where sites can perform local computations. We provide a template script that provides the basic functionalities for performing network file I/O using SCP protocol or via AWS buckets. 

### 3) Matrix Processing Library:
COLLAGENE provides a suite of plaintext and encrypted matrix processing tools that can be used to perform matrix operations such as matrix arithmetic.

As majority of the methods in bioinformatics rely on extensive matrix operations, the matrix library establishes a basic set of tools that can be used for processing small-to-medium sizes matrices.

## Usage Examples for Different Components
We provide examples of using key-setup, encryption/decryption, matrix library, network file I/O (and others) under *examples/* directory. These examples demonstrate the basic blocks for building more complex collaborative analysis pipelines.

## Use Cases:
We provide a full implementation of federated GWAS and meta-analysis for binary traits under *use_cases/* directory. The implementation includes all components that can be readily deployed. A public docker image is provided to install and test these implementations (See *installation/* and *docker/* for more details.).

## API and CLI Documentation
Documentation for the options that can be used in secure pipelines can be found under *API/* directory.

## How do I build new collaborative analysis pipelines with COLLAGENE?
As with any new tool development, Any pipeline that will be implemented in a collaborative/federated setting must be first formulated so that it can be implemented with COLLAGENE's existing modules.

Thus it is necessary for developers to convert their implementations into blocks of modules that can be implemented using COLLAGENE's tools. *Examples/* and *use_cases/* include numerous demonstrations of how these modules can be developed among a number of collaborators.

It is important to not share any individual level data (even when they are encrypted), use masking of aggregated datasets when they are decrypted as intermediate datasets.

We are currently working on providing further recommendations for secure pipeline development for different types of analyses.


