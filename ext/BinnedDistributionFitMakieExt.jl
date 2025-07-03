module BinnedDistributionFitMakieExt
using BinnedDistributionFit
import BinnedDistributionFit: plotthing, plotthing!

using Makie, ComponentArrays

@recipe PlotThing (LS, ps) begin
    cycle = [:color]
    Makie.mixin_generic_plot_attributes()...
end

function get_pdf(pdfs::ExtendPdf)
    return [pdfs]
end

function get_pdf(pdfs::SumOfPdfs)
    return pdfs.pdfs
end

function Makie.plot!(input::PlotThing)
    LS = input.LS[]
    ps = input.ps[]
    stairs!(input, LS.d_hist; label="Data", color=:black)
    pdfs = get_pdf(LS.pdf)
    vks = valkeys(ps)

    normalizations = ps[vks[begin]]

    predictions_of_pdfs = map(pdfs, vks[(begin + 1):end]) do ex_pdf, pas
        x = LS.bcs
        p0 = getproperty(ps, pas)
        oneD_func(x) = BinnedDistributionFit.scalar_eval(ex_pdf, x, p0)
        numerical_int = BinnedDistributionFit._integrate(oneD_func, LS.d_hist, BinnedDistributionFit.QuadGKIntegrator())
        predictions = BinnedDistributionFit.vector_eval(ex_pdf, x, p0)
        return predictions ./ numerical_int
    end

    for (N, ys) in zip(normalizations, predictions_of_pdfs)
        lines!(input, LS.bcs, N * ys)
    end

    if length(pdfs) ≠ 1
        lines!(input, LS.bcs, sum(predictions_of_pdfs))
    end
    return input
end


end

