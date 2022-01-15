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
## These packages have more recent versions available.
## It is recommended to update all of them.
## Which would you like to update?
## 
## 1: All                            
## 2: CRAN packages only             
## 3: None                           
## 4: fansi (0.5.0  -> 1.0.2 ) [CRAN]
## 5: Rcpp  (1.0.7  -> 1.0.8 ) [CRAN]
## 6: rjson (0.2.20 -> 0.2.21) [CRAN]
## 
## fansi (0.5.0  -> 1.0.2 ) [CRAN]
## Rcpp  (1.0.7  -> 1.0.8 ) [CRAN]
## rjson (0.2.20 -> 0.2.21) [CRAN]
```

```
## Skipping 1 packages ahead of CRAN: vegan
```

```
## Installing 3 packages: fansi, Rcpp, rjson
```

```
## Updating HTML index of packages in '.Library'
```

```
## Making 'packages.html' ... done
```

```
##      checking for file ‘/tmp/RtmphtfJZi/remotesb409667951dc/microbiome-miaTime-84b9ffd/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/RtmphtfJZi/remotesb409667951dc/microbiome-miaTime-84b9ffd/DESCRIPTION’
##   ─  preparing ‘miaTime’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
## ─  looking to see if a ‘data/datalist’ file should be added
##   ─  building ‘miaTime_0.1.3.tar.gz’
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
  list(title = "miaTime R package", author = list(list(given = "Yagmur", family = "Simsek", role = NULL, email = "yagmur.simsek.98@gmail.com", comment = NULL), list(given = "Leo", family = "Lahti", role = NULL, email = NULL, comment = NULL)), url = "microbiome.github.io/miaTime", note = "Version 0.1.3").

# Code of conduct

Please note that the project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
