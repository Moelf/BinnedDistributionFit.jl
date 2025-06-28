# BinnedDistributionFit

[![Build Status](https://github.com/Moelf/BinnedDistributionFit.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/Moelf/BinnedDistributionFit.jl/actions/workflows/Test.yml?query=branch%3Amain)

## [Work In Progress]
Everything is subject to change, including type names and type hierarchy.

## The Scope of this Package
We try to cover a core set of usages common in HEP, namely, fitting one or multiple distributions to binned data, using
Chi2  or Likelihood objective.

Conceptually, given user defined function ("shape"), we can treat it as if it's a PDF by numerically normalizing it at
every evaluating point.

## Sum of user-defined PDFs
```julia
using BinnedDistributionFit, FHist

# user defined functions
julia> f(x,ps) = x + ps[1] + ps[2]*x

julia> g(x,ps) = ps[1] + ps[3]

# SumOfPdfs
julia> sd = ExtendPdf(f, (1,3)) + ExtendPdf(g, (1,3))

# data histogram
julia> h = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
edges: [1.0, 2.0, 3.0]
bin counts: [2.0, 4.0]
total count: 6.0

julia> NF = RooFitNLL_functor(sd, h; num_integrator=BinnedDistributionFit.QuadGKIntegrator())

# the functor takes vector of N+1 vectors as input, the first vector keeps track of the normalization for each of the N pdfs, the rest are the parameters for each user-defined pdf
julia> NF([[2.0, 3.0], [6.0, 7.0], [1.0, 5.0, 9.0]])
-0.627546320921633
```

## Math for One PDF fitting to Data: (i.e. `ExtendPdf`)

For main.py and main.jl (just the ExtendPdf)

if you change the observed bincounts from `[2.0, 4.0]` to `[4.0, 8.0]`, the NLL becomes `1.3655798792934464`

If on top of that, you also change `p0` from 2.0 to 4.0, you get NLL `-4.952186287425897`

![image](https://github.com/user-attachments/assets/9af8974a-74f9-406d-b822-dc23973fd13f)


## Math for sum of two PDFs fitting to data: (i.e. `SUM::` in RooFit)

![image](https://github.com/user-attachments/assets/5ec32c86-69df-4895-b8b1-7c179c05fb67)
