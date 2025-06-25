using FHist

function Chi2_functor(d::ExtendPdf, data_hist::Hist1D; num_integrator=SimpleSumIntegrator())
    bes, bcs = binedges(data_hist), bincenters(data_hist)
    obs = bincounts(data_hist)
    @assert extrema(bes) == extrema(d.support) "The support of the distribution must match the bin edges of the histogram."
    
    function (norm_and_ps)
        norm, ps... = norm_and_ps
        expected = vector_eval(d, bcs, ps)
        oneD_func(x) = scalar_eval(d, x, ps)
        normedExp =  expected ./ _integrate(oneD_func, data_hist, num_integrator)
        return chi2(obs, norm * normedExp, binerrors(data_hist))
    end
end

function chi2(obs, expected, sd) 
    return sum(@. abs2((obs - expected)/sd))
end

