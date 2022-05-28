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
## Skipping install of 'miaTime' from a github remote, the SHA1 (8b8485ad) has not changed since last install.
##   Use `force = TRUE` to force installation
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
