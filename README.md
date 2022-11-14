
<!-- README.md is generated from README.Rmd. Please edit that file -->

# orthoVisualizer

Rpackage orthoVisualizer allows identification/annotation/quantification
of user-specified ortholog subsequence within each of the DNA sequences
provided in a multi-FASTA file.

## Description

`orthoVisualizer` is a R package to explore the abundance of orthologous
gene subsequence across multiple DNA sequences, and provide
annotations/visualizations for intuitive analysis. Modern orthologue
detection tools specialize in comparing two sequences (species) at a
time and refer to database to look for corresponding orthologues.
However, such tools will not work with newly-discovered sequences that
are yet to be recorded to the database. `orthoVisualizer` addresses this
issue by allowing the user to specify which particular subsequence to
look for in multiple of species (sequences) at once. Such preliminary
analysis and getting the overview of abundance among multiple sequences
will give researchers the idea of which species to look into (those with
higher abundance ratio) for researching their target orthologue, which
is an essential step in orthologous gene / motif exploration.

The `orthoVisualizer` package was developed using
`R version 4.1.1 (2021-08-10)`, `Platform: x86_64-w64-mingw32` and
`Running under: mingw32`.

## Installation

To install the latest version of the package:

``` r
require("devtools")
devtools::install_github("karynkomatsu/orthoVisualizer", build_vignettes = TRUE)
library("orthoVisualizer")
```

Shiny app: Under Construction

## Overview

Provide the following commands, customized to your R package. Then
provide an overview to briefly describe the main components of the
package. Include one image illustrating the overview of the package,
that shows the inputs and outputs. Ensure the image is deposited in the
correct location, as discussed in class. Point the user to vignettes for
a tutorial of your package. E.g., <br> <br> <br>

``` r
ls("package:orthoVisualizer")
data(package = "orthoVisualizer") 
browseVignettes("orthoVisualizer")
```

`orthoVisualizer` contains 4 functions for exploratory analysis of
orthologous gene/motif occurrence in a multi-FASTA file. The
*annotateSeq* function annotates the occurrence of target subsequence
(orthologous gene/motif) in each DNA sequence. Specifically, it prints
all lines of DNA sequences in FASTA file and display an underline of
“X====Y” for each subsequence occurrence. “=” are the letter that
matches the subsequence, and “X” and “Y” are nucleotide that comes
before / after the subseqence, correspondingly. The *quantSeq* returns
tibble with quantity values related to the subsequence occurrence for
each of the DNA sequences in FASTA file. (ex. “How many times the
subsequence occurred”). *freqSeq* is a function that generates barplot
showing frequency of subsequence occurrence, and each bar represents DNA
sequence in FASTA file (ie. If there were 5 sequences in FASTA file,
there would be 5 bars where first sequence in the FASTA is leftmost bar,
second sequence in the FASTA is second-left bar…). The *freqRatioSeq*
function works similar, except each bar height represents the frequency
ratio (“Number of occurrence” / “Length of sequence over Length of
target subsequence”). Refer to package vignettes for more details. An
overview of the package is illustrated below.

![](./inst/extdata/SILVA_A_A1.png)

## Contributions

contributions from the author and contributions from other
packages/sources for each function. Remember your individual
contributions to the package are important. E.g., The author of the
package is Karyn Komatsu. All functions makes use of `biocManager` R
package to install Bioconductor packages and `BiocGenerics` to re-format
strings of DNA sequence in FASTA file for easier analysis. The
*annotateSeq* function makes use of map function from `mclust` R package
to generate information criteria values. The Integrated Complete
Likelihood (ICL) values are calculated using a function written by the
author. The `stats` R package is used for generating multinomially
distributed random number vectors. Part of the code for
*InfCriteriaCalculation* function has been taken from `<NamePackage>` R
package. (Section of the borrowed code should be clearly indicated and
referenced in the InfCriteriaCalculation R script). The
*InfCriteriaPlot* function makes use of the `graphics` R package.
*NormFactors* function uses Trimmed Mean of M-values (TMM) as
implemented in `edgeR` R package.

## References

Provide full references for all sources used, including for the packages
mentioned under ‘Contributions’, in one format. E.g., <br> <br> <br>
Akaike, H. (1973). Information theory and an extension of the maximum
likelihood principle. In *Second International Symposium on Information
Theory*, New York, USA, 267–281. Springer Verlag.
<https://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15>.

Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture
model for clustering with the integrated classification likelihood.
*IEEE Transactions on Pattern Analysis and Machine Intelligence* 22.
<https://hal.inria.fr/inria-00073163/document>

BioRender. (2020). Image created by Silva, A. Retrieved October 30,
2020, from <https://app.biorender.com/>

McCarthy, D. J., Chen Y. and Smyth, G. K. (2012). Differential
expression analysis of multifactor RNA-Seq experiments with respect to
biological variation. *Nucleic Acids Research* 40. 4288-4297.
<https://pubmed.ncbi.nlm.nih.gov/22287627/>

R Core Team (2021). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>

Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
Statistics* 6, 461–464.
<https://projecteuclid.org/euclid.aos/1176344136>.

Scrucca, L., Fop, M., Murphy, T. B. and Raftery, A. E. (2016) mclust 5:
clustering, classification and density estimation using Gaussian finite
mixture models. *The R Journal* 8(1), 289-317.
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096736/>

Wickham, H. and Bryan, J. (2019). *R Packages* (2nd edition). Newton,
Massachusetts: O’Reilly Media. <https://r-pkgs.org/>

## Acknowledgements

Provide the following text, customized to your R package. E.g., <br>
<br> <br> This package was developed as part of an assessment for
2019-2022 BCB410H: Applied Bioinformatics course at the University of
Toronto, Toronto, CANADA. `orthoVisualizer` welcomes issues, enhancement
requests, and other contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/orthoVisualizer/issues). Many
thanks to those who provided feedback to improve this package.
