module BinnedDistributionFitMakieExt
using BinnedDistributionFit
import BinnedDistributionFit: plotthing, plotthing!
# whatever the else goes here idk

using Makie, ComponentArrays

@recipe PlotThing (LS, ps) begin
    d = 1
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
    plot!(input, LS.d_hist)
    pdfs = get_pdf(LS.pdf)
    vks = valkeys(ps)
    predictions_of_pdfs =  map(pdfs, vks[begin+1:end]) do ex_pdf, pas
        BinnedDistributionFit.vector_eval(ex_pdf, LS.bcs, getproperty(ps, pas))
    end
    for ys in predictions_of_pdfs
         lines!(input, LS.bcs, ys)
    end

    if length(pdfs) ≠ 1
        lines!(input, LS.bcs, sum(predictions_of_pdfs); label="Sum")
    end
    return input
end


end