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

## Negative Log Likelihood function

The negative log likelihood function evaluates
```math
-\text{log}(\mathcal{L})=-\sum_{i=1}^{N_{bins}}y_i\text{log}\left(\frac{\text{UserPdf}(x_i,\text{parameters})}{\int_S\text{UserPdf}(x_i,\text{parameters})dx}\right)-\left(n_{\text{observed}}\text{log}N_{\text{expected}}-N_{\text{expected}}\right)
```
where ``n_{\text{observed}}=\sum y_i`` and ``N_{\text{expected}}`` is the sum of all of the norms of each pdf.
