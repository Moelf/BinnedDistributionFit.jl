@testitem "Cross-check chi2 and LikelihoodSpec" begin

    using Optimization, ForwardDiff, FHist, ComponentArrays

    function compare_chi2_NLL(d, h1, ps)
        support = extrema(binedges(h1))
        NF = BinnedDistributionFit.LikelihoodSpec(ExtendPdf(d, support), h1; loss_type = BinnedDistributionFit.CSQ())

        optf = OptimizationFunction(NF, AutoForwardDiff())
        prob = OptimizationProblem(optf, ComponentArray(norms = [integral(h1)], p1 = ps))
        sol = solve(prob, Optimization.LBFGS())

        NLL = BinnedDistributionFit.LikelihoodSpec(ExtendPdf(d, support), h1)
        optf2 = OptimizationFunction(NLL, AutoForwardDiff())
        prob2 = OptimizationProblem(optf2, ComponentArray(norms = [integral(h1)], p1 = ps))
        sol2 = solve(prob2, Optimization.LBFGS())

        return sol.u, sol2.u
    end

    d(x, ps) = ps[1]x + ps[2]

    h1 = Hist1D(; binedges = 1:2:21, bincounts = [1:4; 1; 6:10], sumw2 = 1:10)

    @test isapprox(compare_chi2_NLL(d, h1, [1, 1])...) rtol = 0.01

end