"""
    RooFitNLL(normalized_predictions, observations; normalization)

Computes the negative log-likelihood for a binned distribution fit.
- `normalized_predictions`: A vector of normalized prediction at each bin-center.
- `observations`: A vector of observed values at each bin-center.
- `normalization`: A scalar value for normalization. Conceptually equal to the pdf generating `predictions` integrated over the domain.
This is an internal function which is automatically used in `LikelihoodSpec`.
"""
function RooFitNLL_internal(normalized_predictions, observations; normalization)
    if any(≤(0), normalized_predictions)
        return Inf
    end
    per_bin_term = sum(@. observations * log(normalized_predictions))
    extend_term = sum(observations) * log(normalization) - normalization

    return -per_bin_term - extend_term
end

#=
"""
    RooFitNLL_functor(d::ExtendPdf, data_hist::Hist1D; num_integrator=SimpleSumIntegrator())

Given an pdf and a data histogram, return a callable function that evaluates to the negative log likelihood with respect to the data histogram within the `d.support` domain.

The callable function should be called with N+1 parameter, the +1 being the first argument represening the overall # of events. For example, if user-defined function is `f(x, p) = x*p[1] + x^2*p[2]`, then the functor should be called with `[norm, p1, p2]`.

## Examples

```julia
julia> f(x,p) = x + p[1] + p[2]*x

julia> d = ExtendPdf(f, (1,3))

julia> BinnedDistributionFit.scalar_eval(d, 1, [2,3])
6

julia> NF = RooFitNLL_functor(d, h);

# the 10. is overall normalization
# [2.0, 3.0] are passed as `p` to `f(x,p)`
julia> NF([10., 2.0, 3.0])
0.060373400847997694
```
"""

"""
    RooFitNLL_functor(sumpdfs::SumOfPdfs, data_hist::Hist1D; num_integrator=SimpleSumIntegrator())

Given an `SumOfPdfs` and a data histogram, return a callable function that evaluates to the negative log likelihood with
respect to the data histogram within the `d.support` domain.

The callable function should be called with a vector of vector as input, specifically `[norms, ps1, ps2, ...]`, where `norms` is the normalization for each of the pdfs in order, and `ps1, ps2, ...` are the parameters for each pdf function, also in order.

## Examples

```julia
julia> f(x,ps) = x + ps[1] + ps[2]*x
julia> g(x,ps) = ps[1] + ps[3]

julia> d1 = ExtendPdf(f, (1,3))

julia> d2 = ExtendPdf(g, (1,3))

julia> sd = d1+d2

julia> h = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
edges: [1.0, 2.0, 3.0]
bin counts: [2.0, 4.0]
total count: 6.0

julia> NF = RooFitNLL_functor(sd, h; num_integrator=BinnedDistributionFit.QuadGKIntegrator())

julia> NF([[2.0, 3.0], [6.0, 7.0], [1.0, 5.0, 9.0]])
-0.627546320921633
```
"""
=#


abstract type LossFunction end
struct NLL <: LossFunction end
struct CSQ <: LossFunction end

"""
    LikelihoodSpec(pdf, d_hist; loss_type=NLL(), num_int=SimpleSumIntegrator())

A type constructor for a set of pdfs and binned histogram data. Method of calculating loss function and integration method may be specified in keyword arguments.
`pdf` must be an `ExtendPdf` or `SumOfPdfs`. `d_hist` 
"""
struct LikelihoodSpec{LF <: LossFunction, p <: AbstractPdf, h <: Hist1D, V <: AbstractVector, NI <: NumericalIntegrator} <: Function
    loss_type::LF
    pdf::p
    d_hist::h
    num_int::NI
    bcs::V
    function LikelihoodSpec(pdf, d_hist; loss_type = NLL(), num_int = SimpleSumIntegrator())
        bcs = bincenters(d_hist)
        p = typeof(pdf)
        h = typeof(d_hist)
        V = typeof(bcs)
        NI = typeof(num_int)
        LF = typeof(loss_type)
        return new{LF, p, h, V, NI}(loss_type, pdf, d_hist, num_int, bcs)
    end
end

function get_pdf(pdfs::ExtendPdf)
    return [pdfs]
end

function get_pdf(pdfs::SumOfPdfs)
    return pdfs.pdfs
end

function (NLL2::LikelihoodSpec{<:NLL})(nps::Vector; kw...)
    norms, vps... = nps
    overall_norm = sum(norms)
    fractions_of_pdfs = norms ./ overall_norm
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)
    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vps) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vps) - 1)"))
    end
    integrals_of_pdfs = map(pdfs, vps) do d, p0
        oneD_func(x) = scalar_eval(d, x, p0; kw...)
        _integrate(oneD_func, NLL2.d_hist, NLL2.num_int; kw...)
    end
    predictions_of_pdfs = map(pdfs, vps) do d, p0
        vector_eval(d, NLL2.bcs, p0, kw...)
    end
    normalized_predictions = sum(fractions_of_pdfs .* predictions_of_pdfs ./ integrals_of_pdfs)
    return RooFitNLL_internal(normalized_predictions, NLL2.d_hist.bincounts; normalization = overall_norm)
end

function (NLL2::LikelihoodSpec{<:NLL})(nps::ComponentVector; kw...)
    vks = valkeys(nps)
    norms = nps[vks[begin]]
    overall_norm = sum(norms)
    fractions_of_pdfs = norms ./ overall_norm
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)
    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end
    integrals_of_pdfs = map(pdfs, vks[(begin + 1):end]) do d, ps
        p0 = getproperty(nps, ps)
        oneD_func(x) = scalar_eval(d, x, p0; kw...)
        _integrate(oneD_func, NLL2.d_hist, NLL2.num_int; kw...)
    end
    predictions_of_pdfs = map(pdfs, vks[(begin + 1):end]) do d, ps
        p0 = getproperty(nps, ps)
        vector_eval(d, NLL2.bcs, p0, kw...)
    end
    normalized_predictions = sum(fractions_of_pdfs .* predictions_of_pdfs ./ integrals_of_pdfs)
    return RooFitNLL_internal(normalized_predictions, NLL2.d_hist.bincounts; normalization = overall_norm)
end

function (NLL2::LikelihoodSpec{<:NLL})(nps::Vector, _dummy; kw...)
    norms, vps... = nps
    overall_norm = sum(norms)
    fractions_of_pdfs = norms ./ overall_norm
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)
    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vps) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vps) - 1)"))
    end
    integrals_of_pdfs = map(pdfs, vps) do d, p0
        oneD_func(x) = scalar_eval(d, x, p0; kw...)
        _integrate(oneD_func, NLL2.d_hist, NLL2.num_int; kw...)
    end
    predictions_of_pdfs = map(pdfs, vps) do d, p0
        vector_eval(d, NLL2.bcs, p0, kw...)
    end
    normalized_predictions = sum(fractions_of_pdfs .* predictions_of_pdfs ./ integrals_of_pdfs)
    return RooFitNLL_internal(normalized_predictions, NLL2.d_hist.bincounts; normalization = overall_norm)
end

function (NLL2::LikelihoodSpec{<:NLL})(nps::ComponentVector, _dummy; kw...)
    vks = valkeys(nps)
    norms = nps[vks[begin]]
    overall_norm = sum(norms)
    fractions_of_pdfs = norms ./ overall_norm
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)
    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end
    integrals_of_pdfs = map(pdfs, vks[(begin + 1):end]) do d, ps
        p0 = getproperty(nps, ps)
        oneD_func(x) = scalar_eval(d, x, p0; kw...)
        _integrate(oneD_func, NLL2.d_hist, NLL2.num_int; kw...)
    end
    predictions_of_pdfs = map(pdfs, vks[(begin + 1):end]) do d, ps
        p0 = getproperty(nps, ps)
        vector_eval(d, NLL2.bcs, p0, kw...)
    end
    normalized_predictions = sum(fractions_of_pdfs .* predictions_of_pdfs ./ integrals_of_pdfs)
    return RooFitNLL_internal(normalized_predictions, NLL2.d_hist.bincounts; normalization = overall_norm)
end
