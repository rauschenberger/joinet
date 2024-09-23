
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/joinet?svg=true)](https://ci.appveyor.com/project/rauschenberger/joinet)
[![R-CMD-check](https://github.com/rauschenberger/joinet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rauschenberger/joinet/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rauschenberger/joinet/graph/badge.svg)](https://app.codecov.io/gh/rauschenberger/joinet)

## Scope

Multivariate elastic net regression through stacked generalisation (extending the [R](https://cran.r-project.org) package [glmnet](https://CRAN.R-project.org/package=glmnet)).

## Installation

Install the current release from
[CRAN](https://CRAN.R-project.org/package=joinet):

``` r
install.packages("joinet")
```

or the latest development version from
[GitHub](https://github.com/rauschenberger/joinet):

``` r
#install.packages("remotes")
remotes::install_github("rauschenberger/joinet")
```

## Reference

Armin Rauschenberger
[![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801)
and Enrico Glaab
[![EG](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3977-7469)
(2021). "Predicting correlated outcomes from molecular data”. *Bioinformatics* 37(21):3889–3895.
[doi: 10.1093/bioinformatics/btab576](https://doi.org/10.1093/bioinformatics/btab576).

[![CRAN version](https://www.r-pkg.org/badges/version/joinet)](https://CRAN.R-project.org/package=joinet)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/joinet)](https://CRAN.R-project.org/package=joinet)
[![Total CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/joinet)](https://CRAN.R-project.org/package=joinet)

## Disclaimer

The R package `joinet` implements multivariate elastic net regression through stacked generalisation ([Rauschenberger et al., 2021](https://doi.org/10.1093/bioinformatics/btab576)).

Copyright &copy; 2019 Armin Rauschenberger, University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
