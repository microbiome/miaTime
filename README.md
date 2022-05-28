# Microbiome time series manipulation with miaTime

This R package can be used to analyse time series data for microbial
communities. The package is part of [miaverse](https://microbiome.github.io/), 
and is based on the `(Tree)SummarizedExperiment` data container.

## Installation
 
The package can be directly installed from R command line.


```r
devtools::install_github("microbiome/miaTime")
```

```
## Downloading GitHub repo microbiome/miaTime@HEAD
```

```
## These packages have more recent versions available.
## It is recommended to update all of them.
## Which would you like to update?
## 
##  1: All                                           
##  2: CRAN packages only                            
##  3: None                                          
##  4: GenomeInfoDb (1.32.1     -> 1.32.2    ) [CRAN]
##  5: BiocParallel (1.30.0     -> 1.30.2    ) [CRAN]
##  6: scuttle      (1.6.0      -> 1.6.2     ) [CRAN]
##  7: V8           (4.1.0      -> 4.2.0     ) [CRAN]
##  8: circlize     (0.4.14     -> 0.4.15    ) [CRAN]
##  9: openssl      (2.0.0      -> 2.0.2     ) [CRAN]
## 10: edgeR        (3.38.0     -> 3.38.1    ) [CRAN]
## 11: limma        (3.52.0     -> 3.52.1    ) [CRAN]
## 12: RcppArmad... (0.11.0.0.0 -> 0.11.1.1.0) [CRAN]
## 13: mia          (1.3.24     -> 1.4.0     ) [CRAN]
## 
## GenomeInfoDb (1.32.1     -> 1.32.2    ) [CRAN]
## BiocParallel (1.30.0     -> 1.30.2    ) [CRAN]
## scuttle      (1.6.0      -> 1.6.2     ) [CRAN]
## V8           (4.1.0      -> 4.2.0     ) [CRAN]
## circlize     (0.4.14     -> 0.4.15    ) [CRAN]
## openssl      (2.0.0      -> 2.0.2     ) [CRAN]
## edgeR        (3.38.0     -> 3.38.1    ) [CRAN]
## limma        (3.52.0     -> 3.52.1    ) [CRAN]
## RcppArmad... (0.11.0.0.0 -> 0.11.1.1.0) [CRAN]
## mia          (1.3.24     -> 1.4.0     ) [CRAN]
```

```
## Installing 10 packages: GenomeInfoDb, BiocParallel, scuttle, V8, circlize, openssl, edgeR, limma, RcppArmadillo, mia
```

```
## Updating HTML index of packages in '.Library'
```

```
## Making 'packages.html' ... done
```

```
##      checking for file ‘/tmp/Rtmpw1Y8ZI/remotes2826e599bbab0/microbiome-miaTime-8b8485a/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/Rtmpw1Y8ZI/remotes2826e599bbab0/microbiome-miaTime-8b8485a/DESCRIPTION’
##   ─  preparing ‘miaTime’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
## ─  looking to see if a ‘data/datalist’ file should be added
##   ─  building ‘miaTime_0.1.6.tar.gz’
##      
## 
```

```r
library(miaTime)
```

### Contributions and acknowledgments

You can find us online from [Gitter](https://gitter.im/microbiome/miaverse).

Contributions are very welcome through issues and pull requests at the
[development site](https://github.com/microbiome/miaTime). We follow a git
flow kind of approach. Development version should be done against the
`main` branch and then merged to `release` for release.
(https://guides.github.com/introduction/flow/)

**Kindly cite this work**. For citation details, see R command

```r
citation("miaTime")  
```

```
## 
## Kindly cite the miaTime R package as follows:
## 
##   (C) Yagmur Simsek and Leo Lahti. miaTime R package Version 0.1.7
##   Package URL: microbiome.github.io/miaTime
## 
## A BibTeX entry for LaTeX users is
## 
##   @Misc{,
##     title = {miaTime R package},
##     author = {Yagmur Simsek and Leo Lahti},
##     url = {microbiome.github.io/miaTime},
##     note = {Version 0.1.7},
##   }
```

# Code of conduct

Please note that the project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
