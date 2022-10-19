Statial
======================================================

Overview
--------

**Statial** is a suite of functions for identifying changes in cell state. 
The functionality provided by Statial provides robust quantification of cell 
type localisation which are invariant to changes in tissue structure. In 
addition to this Statial uncovers changes in marker expression associated with
varying levels of localisation. These features can be used to explore how the
structure and function of different cell types may be altered by the agents 
they are surrounded with.


Installation
--------
Installing the package from Bioconductor.
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("Statial")
```

Otherwise, install the development version from GitHub.

```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("SydneyBioX/Statial")
library(Statial)
```

### Submitting an issue or feature request

`Statial` is still under active development. We would greatly appreciate any and 
all feedback related to the package.

* R package related issues should be raised [here](https://github.com/SydneyBioX/Statial/issues).
* For general questions and feedback, please contact us directly via [ellis.patrick@sydney.edu.au](mailto:ellis.patrick@sydney.edu.au).


## Author

* **Farhan Ameen**
* **Sourish Iyengar**
* **Shila Ghazanfar** - [@shazanfar](https://twitter.com/shazanfar)
* **Ellis Patrick**  - [@TheEllisPatrick](https://twitter.com/TheEllisPatrick)
