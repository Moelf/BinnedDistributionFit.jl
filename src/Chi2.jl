using FHist

"""
    chi2(obs, expected, sd)
computes stuff in stuff and stuff
"""
function chi2(obs, expected, sd)
    return sum(@. abs2((obs - expected) / sd))
end

function (NLL2::LikelihoodSpec{<:CSQ})(nps::ComponentVector; kw...)
    ber = binerrors(NLL2.d_hist)
    vks = valkeys(nps)
    norms = nps[vks[begin]]
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)
    if length(norms) != length(pdfs)
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
    overall_norms = sum(norms)
    normed_predictions = sum(abs2.(norms) / overall_norms .* predictions_of_pdfs ./ integrals_of_pdfs)
    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, ber)
end

function (NLL2::LikelihoodSpec{<:CSQ})(nps::ComponentVector, _dummy; kw...)
    ber = binerrors(NLL2.d_hist)
    vks = valkeys(nps)
    norms = nps[vks[begin]]
    pdfs = get_pdf(NLL2.pdf)
    Npdfs = length(pdfs)
    if length(norms) != length(pdfs)
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
    overall_norms = sum(norms)
    normed_predictions = sum(abs2.(norms) / overall_norms .* predictions_of_pdfs ./ integrals_of_pdfs)
    return BinnedDistributionFit.chi2(NLL2.d_hist.bincounts, normed_predictions, ber)
end