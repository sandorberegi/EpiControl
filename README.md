# EpiControl

This code is related to the paper (in preparation):
Sandor Beregi, Sangeeta Bhatia, Anne Cori and Kris V. Parag, **EpiControl: automating and optimising epidemic interventions**

Sandor Beregi, Sangeeta Bhatia, Anne Cori and Kris V. Parag, **EpiControl: a data-driven tool for optimising epidemic interventions and automating scenario planning to support real-time response**, MEDRXIV/2025/340271, [Link](https://www.medrxiv.org/content/10.1101/2025.11.17.25340271v1)

The realted methodology is also used in:
Sandor Beregi, Kris V. Parag, **Optimal algorithms for controlling infectious diseases in real time using noisy infection data**. 
PLOS Computational Biology 21 (9), e1013426, [Link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013426)

This code was intially adapted from the [Epicont.jl]([src/EpiCont.jl](https://github.com/sandorberegi/Epidemic-control-with-noisy-real-time-data)) Julia package.

Dependencies are the following:

The code requires R version >4.0.0.

Packages:

| Package     | Version |
|-------------|---------|
| dplyr       | 1.1.4   |
| ggplot2    | 3.5.1   |
| knitr      | 1.45    |
| pbapply    | 1.7-2   |
| rmarkdown  | 2.26    |
| stats      | 4.3.1   |
| testthat   | 3.2.1.1 |
| utils      | 4.3.1   |
| VGAM       | 1.1-10  |


This software was created using MacOS, latest tested version: MacOS Tahoe 26.2.
If you run into issues with this software under any other OS please do not hesitate to report it to us.

## Installation

This package is developed as an **R package hosted on GitHub**.

To install it, you need the **remotes** package:

```r
install.packages("remotes")

remotes::install_github("sandorberegi/EpiControl")

```
If you want to build vignettes locally:

```r
remotes::install_github(
  "sandorberegi/EpiControl",
  build_vignettes = TRUE
)
```

On a standard desktop or laptop computer with a recent version of R (≥ 4.3) and an internet connection, installation typically takes:

1–3 minutes if all dependencies are already installed

3–5 minutes if dependencies need to be downloaded and compiled
