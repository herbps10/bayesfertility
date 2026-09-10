# bayesfertility
R package for Bayesian fertility models.

## Installation

`bayesfertility` uses [`cmdstanr`](https://mc-stan.org/cmdstanr/) to fit its
Stan models. The models are compiled when the package is installed, so you need
a working `cmdstanr` / `CmdStan` setup in place *before* installing this
package. Follow these steps, in order:

### 1. Install `cmdstanr`

`cmdstanr` is not on CRAN, so install it instead from the Stan R-universe:

```r
install.packages(
    "cmdstanr",
    repos = c("htps://stan-dev.r-universe.dev", getOption("repos"))
)
```

### 2. Install CmdStan

This downloads and installs the CmdStan command-line tool:

```r
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()
```

If `check_cmdstan_toolchain()` reports any problems, follow the [cmdstanr
toolchain
guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html#installing-cmdstan).

### 3. Install `bayesfertility`

```r
# install.packages("remotes")
remotes::install_github("herbps10/bayesfertility")
```
