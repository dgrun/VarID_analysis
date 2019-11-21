# VarID analysis

from: Grün D (2019) Revealing Dynamics of Gene Expression Variability in Cell State Space. Nature Methods Nov 18. doi: 10.1038/s41592-019-0632-3

This repository contains R source code for reproducing the VarID
analysis presented in Grün D. XXX.

VarID is part of the RaceID package starting with RaceID v0.1.4

Two data sets were analysed:
1) mouse hematopoietic progenitors (data from Tusi, B. K. et
al. Population snapshots predict early haematopoietic and erythroid
hierarchies. Nature 555, 54–60 (2018))
2) mouse intestinal epithelial cells (data from Haber, A. L. et al. A
single-cell survey of the small intestinal epithelium. Nature 551,
333–339 (2017))

For each of these datasets, a separate R script and input data in
*.rds format are provided


## Hematopoietic progenitors

Input data:
```
inputData_hematopoiesis.rds
```
Source code:
```
VarID_hematopoiesis.R
```

## Intestinal epithelial cells

Input data:
```
inputData_intestine.rds
```
Source code:
```
VarID_intestine.R
```

## Reference:

Grün D (2019) Revealing Dynamics of Gene Expression Variability in Cell State Space. Nature Methods Nov 18. doi: 10.1038/s41592-019-0632-3
