
<!-- README.md is generated from README.Rmd. Please edit that file -->
learnMotifs
===========

[learnMotifs](https://github.com/mPloenzke/learnMotifs) implements the weight constrained and regularized convolutional filter layers for learning sequence motifs as information gain matrices (IGMs) or position weight matrices (PWMs) directly. The package relies upon [Keras](https://keras.rstudio.com/) to perform much of the heavy lifting. It is being updated to include more vignettes, features, and broader generalizability so bugs may exist. If you encounter a bug, please file a minimal reproducible example on [github](https://github.com/mPloenzke/learnMotif/issues).

Installation
------------

Begin by installing the [Keras](https://keras.rstudio.com/) package.

``` r
install.packages('keras')
```

Additionally, we recommend installing the TensorFlow backend.

``` r
library(keras)
```

Install keras if it has not already been installed.

``` r
install_keras()
```

You may then install learnMotif from github with:

``` r
# install.packages("devtools")
devtools::install_github("mPloenzke/learnMotifs")
```

Getting Started
---------------

To get started, please read the intro vignette: `vignette("learn_toy_motifs", package = "learnMotif")`.

The core convolution functions create Keras layers to be used in any Keras sequential model. Several resources are available for learning how to define a sequential model, including [this](https://keras.rstudio.com/articles/tutorial_basic_classification.html) and [this](https://keras.rstudio.com/articles/sequential_model.html).

The and layers provided within [learnMotifs](https://github.com/mPloenzke/learnMotifs) may simply be substituted into the sequential model definition where the from the [Keras](https://keras.rstudio.com/) package would be. The first may be used to learn motifs *de novo* while the latter may be used in conjunction with previously-defined motifs, including those imported from databases such as [JASPAR](jaspar.genereg.net/) or [Hocomoco](http://hocomoco11.autosome.ru/). Helper functions are available to assist in importing motif annotations.

Learning sequence motifs directly is achieved through weight constraints and peri-training weight transformations. Weight regularization is employed to discourage learning redundant motifs. Various function definitions are included to these ends and called by and . may be used in place of max or average pooling and provides extended max-pooling functionality to calculate the top *N* maximums.

Several additional functions are available to aid in the analysis pipeline and improve model interpretability. One is encouraged to view their use cases in the vignettes provided. A brief listing is provided below:

1.  Callback definitions - These are used as model wrappers to output diagnostics peri-training.
2.  Helper functions - Simple functions to import database motifs, set up log directories, one-hot encode nucleotide sequences, etc.
3.  Activation functions - Model interpretability is improved by assessing filter activations. To this end functions are included to manipulate the activation values produced by each motif and optionally plot them with [ggplot2](https://ggplot2.tidyverse.org/).
