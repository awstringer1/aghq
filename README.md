[![Build Status](https://travis-ci.com/awstringer1/aghq.svg?branch=master)](https://travis-ci.com/awstringer1/aghq)

# AGHQ: Adaptive Gauss Hermite Quadrature for Bayesian Inference

This package implements Bayesian inference using Adaptive Gauss-Hermite Quadrature, as specified in [our paper](https://arxiv.org/abs/2102.06801) **Stochastic Convergence Rates and Applications of Adaptive Quadrature in Bayesian Inference** with Yanbo Tang and Blair Bilodeau. See also the accompanying [vignette](https://arxiv.org/abs/2101.04468), available on arXiv.

You can install the development version from Github:

```R
install.packages('devtools')
devtools::install_github('awstringer1/aghq')
```

You can also install the stable version from [CRAN](https://CRAN.R-project.org/package=aghq):

```R
install.packages('aghq')
```

The two papers linked above give a comprehensive overview of the method, application, and theory.
