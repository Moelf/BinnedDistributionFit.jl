```@raw html
---
layout: home
hero:
  name: "BinnedDistributionFit.jl"
  tagline: "Maximum-likelihood fitting of binned distributions"
  image:
    src: logo.svg
    alt: BinnedDistributionFit_Logo
  actions:
    - theme: brand
      text: "Quick Start"
      link: /90-quickstart
    - theme: alt
      text: APIs
      link: /95-reference
---
```

```@meta
CurrentModule = BinnedDistributionFit
```

# BinnedDistributionFit

Documentation for [BinnedDistributionFit](https://github.com/Moelf/BinnedDistributionFit.jl).



# What is BinnedDistributionFit?
BinnedDistributionFit is a package which defines several loss functions to compare probability density functions (pdfs) or sums thereof to specific histogram data. BinnedDistributionFit supports any number of independently normalized pdfs. Either negative log likelihood or chi squared regressions may be used.

# Mathematical Background
`BinnedDistributionFit` supports two methods for calculating likelihood: Negative log likelihood and the chi squared test.

## Negative Log Likelihood Function
Given a histogram and a pdf, we can compare the differences between each observed and expected value over the support (domain) of the pdf using likelihood as a measurement of error. We use negative log likelihood (NLL) since this simplifies the calculation of the overall likelihood by reducing the product of each points' likelihood to the sum of the logs of their likelihoods. We negate this sum to format this properly to act as a minimization problem.

With histogram data ``y`` at each of a set of points ``x_i`` on the support ``S``, we can perform our evaluation like so 
```math
-\ln(\mathcal{L})=-\sum_{i=1}^{N_{\text{bins}}}y_i\ln\left(\sum_{j\ =1}^{N_{\text{pdfs}}}\frac{N_i\cdot\text{UserPdf}_i(x_i,\text{parameters})}{N_{\text{expected}}\cdot\int_S\text{UserPdf}_i(x,\text{parameters})dx}\right)-\left(N_{\text{observed}}\ln N_{\text{expected}}-N_{\text{expected}}\right).
```
where ``N_{\text{observed}}=\sum y_i`` and ``N_{\text{expected}}`` is the sum of all of the norms (``N_i``) of each pdf. ``N_{\text{bins}}`` and ``N_{\text{pdfs}}`` represent the number of bins and the number of pdfs respectively. 

When only one pdf is defined this equation reduces to 
```math
-\ln(\mathcal{L})=-\sum_{i=1}^{N_{\text{bins}}}y_i\ln\left(\frac{\text{UserPdf}_i(x_i,\text{parameters})}{\int_S\text{UserPdf}_i(x,\text{parameters})dx}\right)-\left(N_{\text{observed}}\ln N_{\text{expected}}-N_{\text{expected}}\right).
```

## Chi Squared Function
Compared to the NLL function, the chi squared measure of error is much simpler. However, it has the downside that it requires some measure of the standard error of each bin. This package uses the built-in `Hist1D` field `sumw2' to define the variance (``\sigma^2``) of the histogram. The formula for this method is
```math
\chi^2 = \sum_{i=1}^{N_{\text{bins}}}\frac{\left(y_i-\sum_{j=1}^{N_{\text{pdfs}}}\frac{N_i^2\cdot\text{UserPdf}_i(x_i,\text{parameters})}{N_{\text{expected}}\cdot\int_S\text{UserPdf}_i(x,\text{parameters})dx}\right)^2}{\sigma_i^2}.
```
This formulation uses the same variables as the previous equation.

