using FHist

"""
    chi2(obs, expected, sd)
Computes the reduced chi squared statistic for given vectors of observed and expected values. Returns `Inf` if any standard deviation is not a positive number.

This function is primarily for internal use.
"""
function chi2(obs, expected, sd)
    return all(i -> i > 0, sd) ? sum(@. abs2((obs - expected) / sd)) : Inf
end

"""
    (NLL2::LikelihoodSpec{<:CSQ})(nps::Vector; kw...)
Method of `LikelihoodSpec` structs. Evaluates a given `LikelihoodSpec` via the reduced chi squared statistic with a specific set of parameters `nps`. This method takes `nps` to be a `Vector{Vector}`.
"""
function (NLL2::LikelihoodSpec{<:CSQ})(nps::Vector; kw...)
    norms, vps... = nps
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)

    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end

    normed_predictions = mapreduce(+, pdfs, vps, abs2.(norms) / sum(norms)) do d, ps, fop
        p0 = getproperty(nps, ps)
        iop = _integrate(x -> scalar_eval(d, x, p0; kw...), NLL2.d_hist, NLL2.num_int; kw...)

        (fop .* vector_eval(d, NLL2.bcs, p0; kw...) ./ iop)
    end

    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, binerrors(NLL2.d_hist))
end

function (NLL2::LikelihoodSpec{<:CSQ})(nps::Vector, _dummy; kw...)
    norms, vps... = nps
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)

    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end

    normed_predictions = mapreduce(+, pdfs, vps, abs2.(norms) / sum(norms)) do d, ps, fop
        p0 = getproperty(nps, ps)
        iop = _integrate(x -> scalar_eval(d, x, p0; kw...), NLL2.d_hist, NLL2.num_int; kw...)

        (fop .* vector_eval(d, NLL2.bcs, p0; kw...) ./ iop)
    end

    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, binerrors(NLL2.d_hist))
end
#=
function (NLL2::LikelihoodSpec{<:CSQ})(nps::Vector, fixed::Vector; kw...)
    for i in fixed
        nps[i[1]][i[2]] = i[3]
    end
    norms, vps... = nps
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)

    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end

    normed_predictions = mapreduce(+, pdfs, vps, abs2.(norms) / sum(norms)) do d, ps, fop
        p0 = getproperty(nps, ps)
        iop = _integrate(x -> scalar_eval(d, x, p0; kw...), NLL2.d_hist, NLL2.num_int; kw...)

        (fop .* vector_eval(d, NLL2.bcs, p0; kw...) ./ iop)
    end

    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, binerrors(NLL2.d_hist))
end
=#

"""
    (NLL2::LikelihoodSpec{<:CSQ})(nps::Vector; kw...)
Method of `LikelihoodSpec` structs. Evaluates a given `LikelihoodSpec` via the reduced chi squared statistic with a specific set of parameters `nps`. This method takes `nps` to be a `ComponentArray`.
"""
function (NLL2::LikelihoodSpec{<:CSQ})(nps::ComponentVector; kw...)
    vks = valkeys(nps)
    norms = nps[vks[begin]]
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)

    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end

    normed_predictions = mapreduce(+, pdfs, vks[(begin + 1):end], abs2.(norms) / sum(norms)) do d, ps, fop
        p0 = getproperty(nps, ps)
        iop = _integrate(x -> scalar_eval(d, x, p0; kw...), NLL2.d_hist, NLL2.num_int; kw...)

        (fop .* vector_eval(d, NLL2.bcs, p0; kw...) ./ iop)
    end

    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, binerrors(NLL2.d_hist))
end

function (NLL2::LikelihoodSpec{<:CSQ})(nps::ComponentVector, _dummy; kw...)
    vks = valkeys(nps)
    norms = nps[vks[begin]]
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)

    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end

    normed_predictions = mapreduce(+, pdfs, vks[(begin + 1):end], abs2.(norms) / sum(norms)) do d, ps, fop
        p0 = getproperty(nps, ps)
        iop = _integrate(x -> scalar_eval(d, x, p0; kw...), NLL2.d_hist, NLL2.num_int; kw...)

        (fop .* vector_eval(d, NLL2.bcs, p0; kw...) ./ iop)
    end

    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, binerrors(NLL2.d_hist))
end
#=
function (NLL2::LikelihoodSpec{<:CSQ})(nps::ComponentVector, fixed::Vector; kw...)
    for i in fixed
        nps[i[1]] = i[2]
    end
    vks = valkeys(nps)
    norms = nps[vks[begin]]
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)

    if length(norms) != Npdfs
        throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
    end
    if length(vks) - 1 != Npdfs
        throw(ArgumentError("Expected $Npdfs set$((Npdfs > 1) ? "s" : "") of parameters, got $(length(vks) - 1)"))
    end

    normed_predictions = mapreduce(+, pdfs, vks[(begin + 1):end], abs2.(norms) / sum(norms)) do d, ps, fop
        p0 = getproperty(nps, ps)
        iop = _integrate(x -> scalar_eval(d, x, p0; kw...), NLL2.d_hist, NLL2.num_int; kw...)

        (fop .* vector_eval(d, NLL2.bcs, p0; kw...) ./ iop)
    end

    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, binerrors(NLL2.d_hist))
end
=#
