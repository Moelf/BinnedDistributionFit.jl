@testitem "Cross-check chi2 and NLL" begin

using Optimization, ForwardDiff, FHist

function compare_chi2_NLL(d, h1, ps)
    support = extrema(binedges(h1))
    NF = BinnedDistributionFit.Chi2_functor(ExtendPdf(d, support), h1)
    NF_wrapper(x,_) = NF(x)
    
    optf = OptimizationFunction(NF_wrapper, AutoForwardDiff())
    prob = OptimizationProblem(optf, [integral(h1); ps])
    sol = solve(prob, Optimization.LBFGS())
    
    NLL = BinnedDistributionFit.RooFitNLL_functor(ExtendPdf(d, support), h1)
    NLL_wrapper(x,_) = NLL(x)
    
    optf2 = OptimizationFunction(NLL_wrapper, AutoForwardDiff())
    prob2 = OptimizationProblem(optf2, [integral(h1); ps])
    sol2 = solve(prob2, Optimization.LBFGS())

    return sol.u, sol2.u
end

d(x, ps) = ps[1]x+ps[2]
h1 = Hist1D(; binedges=1:2:21, bincounts=[1:4; 1; 6:10], sumw2=1:10)
@test isapprox(compare_chi2_NLL(d, h1, [1,1])...) rtol=0.01

end
