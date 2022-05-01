# Microbiome time series manipulation with miaTime

This R package can be used to analyse time series data for microbial
communities. The package is part of [miaverse](microbiome.github.io), 
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
## Skipping 1 packages ahead of CRAN: mia
```

```
##      checking for file ‘/tmp/RtmprfvIlK/remotesc2d4406f3bb3/microbiome-miaTime-d3ad218/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/RtmprfvIlK/remotesc2d4406f3bb3/microbiome-miaTime-d3ad218/DESCRIPTION’
##   ─  preparing ‘miaTime’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
## ─  looking to see if a ‘data/datalist’ file should be added
##   ─  building ‘miaTime_0.1.4.tar.gz’
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
  list(title = "miaTime R package", author = list(list(given = "Yagmur", family = "Simsek", role = NULL, email = "yagmur.simsek.98@gmail.com", comment = NULL), list(given = "Leo", family = "Lahti", role = NULL, email = NULL, comment = NULL)), url = "microbiome.github.io/miaTime", note = "Version 0.1.5").

# Code of conduct

Please note that the project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
