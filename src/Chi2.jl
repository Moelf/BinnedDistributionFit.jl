using FHist

function Chi2_functor(d::ExtendPdf, data_hist::Hist1D; num_integrator = SimpleSumIntegrator())
    bes, bcs = binedges(data_hist), bincenters(data_hist)
    obs = bincounts(data_hist)
    @assert extrema(bes) == extrema(d.support) "The support of the distribution must match the bin edges of the histogram."
    
    function (norm_and_ps; kw...)
        norm, ps... = norm_and_ps
        expected = vector_eval(d, bcs, ps; kw...)
        oneD_func(x) = scalar_eval(d, x, ps; kw...)
        normedExp =  expected ./ _integrate(oneD_func, data_hist, num_integrator)
        return chi2(obs, norm * normedExp, binerrors(data_hist))
    end
end

function chi2(obs, expected, sd) 
    return sum(@. abs2((obs - expected)/sd))
end

function Chi2_functor(d::SumOfPdfs, data_hist::Hist1D; num_integrator = SimpleSumIntegrator())
    bes, bcs = binedges(data_hist), bincenters(data_hist)
    obs = bincounts(data_hist)
    @assert extrema(bes) == extrema(d.support) "The support of the distribution must match the bin edges of the histogram."
    function (norms_and_vps; kw...)
        norms, vps... = norms_and_vps

        if length(norms) != length(d.pdfs)
            throw(ArgumentError("Expected $(length(d.pdfs)) normalizations, got $(length(norms))"))
        end
        integrals_of_pdfs = map(d.pdfs, vps) do d, ps
            oneD_func(x) = scalar_eval(d, x, ps; kw...)
            _integrate(oneD_func, data_hist, num_integrator)
        end
        predictions_of_pdfs =  map(d.pdfs, vps) do d, ps
            vector_eval(d, bcs, ps; kw...)
        end
        overall_norms = sum(norms)
        normed_predictions = sum(abs2.(norms) / overall_norms .* predictions_of_pdfs ./ integrals_of_pdfs)
        return chi2(obs, normed_predictions, binerrors(data_hist))
    end
end