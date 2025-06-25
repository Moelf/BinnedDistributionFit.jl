module BinnedDistributionFit

import Distributions: pdf, cdf

export ExtendPdf, SumOfPdfs, RooFitNLL, RooFitNLL_functor, pdf

using FHist: Hist1D, binedges, bincenters, bincounts
using QuadGK: quadgk
using Trapz: trapz

"""
    struct ExtendPdf{F, S}
        func::F
        support::S
    end

Type representing an extended pdf where normalization does not need to equal to 1.
`func` can be any function that can be evaluated within `support`, when calculating Maximum Likelihood, `func` will be numerically normalized over `support`.

## Examples
```julia
julia> f(x, ps) = x^ps[1] + x*ps[2]

julia> d = BinnedDistributionFit.ExtendPdf(f, (1, 2))

julia> BinnedDistributionFit.scalar_eval(d, 2.0, [1,2])
6.0

julia> BinnedDistributionFit.vector_eval(d, [2.0, 3.0], [1,2])
2-element Vector{Float64}:
 6.0
 9.0
```
"""
struct ExtendPdf{F, S}
    func::F
    support::S
end

function scalar_eval(d::ExtendPdf, x, ps=(); kw...)
    d.func(x, ps; kw...)
end

function vector_eval(d::ExtendPdf, x::AbstractVector, ps=(); kw...)
    d.func.(x, Ref(ps); kw...)
end

"""
    struct SumOfPdfs{V, S}
        pdfs::V
        support::S
    end

Type representing a sum of pdfs.

## Examples

```julia
julia> d1 = ExtendPdf((x, _)->x, (1,3))

julia> d2 = ExtendPdf((x, _)->x^2, (1,3))

julia> sd = d1+d2
SumOfPdfs{Vector{ExtendPdf{F, Tuple{Int64, Int64}} where F}, Tuple{Int64, Int64}}(ExtendPdf{F, Tuple{Int64, Int64}} where F[ExtendPdf{var"#13#14", Tuple{Int64, Int64}}(var"#13#14"(), (1, 3)), ExtendPdf{var"#15#16", Tuple{Int64, Int64}}(var"#15#16"(), (1, 3))], (1, 3))

julia> BinnedDistributionFit.scalar_eval(sd, 2.0)
6.0

julia> BinnedDistributionFit.vector_eval(sd, [2.0, 3.0])
2-element Vector{Float64}:
 6.0
12.0

julia> sd2 = d1+d2+d3

julia> BinnedDistributionFit.scalar_eval(sd2, 2.0, [[], [], [1,2]])
ERROR: BoundsError: attempt to access 2-element Vector{Any} at index [3]

julia> BinnedDistributionFit.scalar_eval(sd2, 2.0, [[], [], [1,2,3]])
9.0
```
"""
struct SumOfPdfs{V, S}
    pdfs::V
    support::S
end

function Base.:+(a::ExtendPdf, b::ExtendPdf)
    if a.support != b.support
        error("+: Supports of the two ExtendPdfs are not equal")
    end
    SumOfPdfs([a, b], a.support)
end

function Base.:+(a::SumOfPdfs, b::ExtendPdf)
    if a.support != b.support
        error("+: Supports of the SumOfPdfs $a and ExtendPdf $b are not equal")
    end
    SumOfPdfs([a.pdfs; b], a.support)
end
function Base.:+(a::ExtendPdf, b::SumOfPdfs)
    if a.support != b.support
        error("+: Supports of the ExtendPdf $a and SumOfPdfs $b are not equal")
    end
    SumOfPdfs([a; b.pdfs], a.support)
end
function Base.:+(a::SumOfPdfs, b::SumOfPdfs)
    if a.support != b.support
        error("+: Supports of the two SumOfPdfs $a and $b are not equal")
    end
    SumOfPdfs([a.pdfs; b.pdfs], a.support)
end

"""
    scalar_eval(d::SumOfPdfs, x, vps::Vector{Vector}; kw...)

Evaluate the sum of pdfs at a single point. `ps` needs to be a vector of vectors of parameters because each pdf in the sum has its own set of parameters.

For example, if we have two pdfs `pdf1` and `pdf2` with parameters `params1` and `params2`, respectively, we can evaluate the sum of the pdfs at a point `x` as follows:

```julia
params1 = [1.0, 2.0]
params2 = [3.0, 4.0]
scalar_eval(sumofpdfs([pdf1, pdf2]), x, [params1, params2])
```
"""
function scalar_eval(d::SumOfPdfs, x, vps=fill(nothing, length(d.pdfs)); kw...)
    mapreduce(+, d.pdfs, vps) do d, ps
        scalar_eval(d, x, ps; kw...)
    end
end

function vector_eval(d::SumOfPdfs, x::AbstractVector, vps=fill(nothing, length(d.pdfs)); kw...)
    mapreduce(+, d.pdfs, vps) do d, ps
        vector_eval(d, x, ps; kw...)
    end
end

include("./num_integrator.jl")
include("./RooFitNLL.jl")
include("./Chi2.jl")

end
