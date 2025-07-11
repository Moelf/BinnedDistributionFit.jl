# BinnedDistributionFit

[![Build Status](https://github.com/Moelf/BinnedDistributionFit.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/Moelf/BinnedDistributionFit.jl/actions/workflows/Test.yml?query=branch%3Amain)

## [Work In Progress]
Everything is subject to change, including type names and type hierarchy.

## The Scope of this Package
We try to cover a core set of usages common in HEP, namely, fitting one or multiple distributions to binned data, using
Chi2  or Likelihood objective.

Conceptually, given user defined function ("shape"), we can treat it as if it's a PDF by numerically normalizing it at
every evaluating point. We then package that with a user defined histogram using the FHist package which can be used with an optimization or plotting package.

## Quick Start
The following code provides a framework for optimizing a single function with any number of parameters. The guess must be relatively close to the real parameters. 
```jldoctest
using BinnedDistributionFit, FHist, Optimization, ForwardDiff, ComponentArrays

function fit_pdf_to_hist(hist::Hist1D, pdf, guess::Vector)
    support = extrema(binedges(hist))

    NLL = LikelihoodSpec(ExtendPdf(pdf, support), hist)

    opt_f = OptimizationFunction(NLL, AutoForwardDiff())
    opt_p = OptimizationProblem(opt_f, ComponentArray(norms = integral(hist), p1 = guess))
    sol = solve(opt_p, Optimization.LBFGS())

    return sol.u
end

# Example

hist = Hist1D(; binedges=1:2:21, bincounts=[1:4; 1; 6:10])
pdf_input(x, ps) = ps[1]*x + ps[2]
guess = [5, 5]

fit_pdf_to_hist(hist, pdf_input, guess)

# output

ComponentVector{Float64}(norms = 51.0, p1 = [6.975319732308675, -0.8718290560996027])
```

## Math for One PDF fitting to Data: (i.e. `ExtendPdf`)

For main.py and main.jl (just the ExtendPdf)

if you change the observed bincounts from `[2.0, 4.0]` to `[4.0, 8.0]`, the NLL becomes `1.3655798792934464`

If on top of that, you also change `p0` from 2.0 to 4.0, you get NLL `-4.952186287425897`

![image](https://github.com/user-attachments/assets/9af8974a-74f9-406d-b822-dc23973fd13f)


## Math for sum of two PDFs fitting to data: (i.e. `SUM::` in RooFit)

![image](https://github.com/user-attachments/assets/5ec32c86-69df-4895-b8b1-7c179c05fb67)
