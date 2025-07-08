@testitem "Cross-check chi2 and LikelihoodSpec" begin

using Optimization, ForwardDiff, FHist, ComponentArrays

function compare_chi2_NLL(d, h1, ps)
    support = extrema(binedges(h1))
    NF = BinnedDistributionFit.chi2_functor(ExtendPdf(d, support), h1)
    NF_wrapper(x,_) = NF(x)
    
    optf = OptimizationFunction(NF_wrapper, AutoForwardDiff())
    prob = OptimizationProblem(optf, [integral(h1); ps])
    sol = solve(prob, Optimization.LBFGS())
    
    NLL = BinnedDistributionFit.LikelihoodSpec(ExtendPdf(d, support), h1)
    NLL_wrapper(x,_) = NLL(x)
    optf2 = OptimizationFunction(NLL_wrapper, AutoForwardDiff())
    prob2 = OptimizationProblem(optf2, ComponentArray(a=[integral(h1)],b=ps))
    sol2 = solve(prob2, Optimization.LBFGS())

    return sol.u, sol2.u
end

d(x, ps) = ps[1]x+ps[2]

h1 = Hist1D(; binedges=1:2:21, bincounts=[1:4; 1; 6:10], sumw2=1:10)

@test isapprox(compare_chi2_NLL(d, h1, [1,1])...) rtol=0.01

end

@testitem "optimization test" begin
    using FHist, ComponentArrays, Optimization, ForwardDiff, Distributions
    N = Normal(2, 1.5)
    dat1 = rand(N, 1000)
    E = Exponential(2)
    dat2 = rand(E, 10000)
    h1 = Hist1D(append!(dat1, dat2); binedges= -2:0.05:3)
    f1(x, ps) = (1/ps[1]) * exp(-x/ps[1])
    f2(x, ps) = (1/sqrt(2 *pi*abs2(ps[2]))) * exp(-abs2(x-ps[1])/(2*abs2(ps[2])))
    pdfs = ExtendPdf(f1, (-2,3)) + ExtendPdf(f2, (-2,3))
    NLL = BinnedDistributionFit.LikelihoodSpec(pdfs, h1)
    NLL_wrapper(x,_) = NLL(x)
    optf2 = OptimizationFunction(NLL_wrapper, AutoForwardDiff())
    prob2 = OptimizationProblem(optf2, ComponentArray(norms=[1.0, 1.0], p1=[2.0], p2=[2.0, 1.5]))
    sol = solve(prob2, Optimization.LBFGS())
    @test sol == 0
end