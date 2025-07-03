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
    println(ComponentArray(a=integral(h1),b=ps))
    optf2 = OptimizationFunction(NLL_wrapper, AutoForwardDiff())
    prob2 = OptimizationProblem(optf2, ComponentArray(a=integral(h1),b=ps))
    sol2 = solve(prob2, Optimization.LBFGS())

    return sol.u, sol2.u
end

d(x, ps) = ps[1]x+ps[2]

h1 = Hist1D(; binedges=1:2:21, bincounts=[1:4; 1; 6:10], sumw2=1:10)

@test isapprox(compare_chi2_NLL(d, h1, [1,1])...) rtol=0.01

end


@testitem "LikelihoodSpec ExtendPdf" begin
    using FHist, ComponentArrays
    d1 = ExtendPdf((x, _)->x, (1,3))
    h1 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF1 = BinnedDistributionFit.LikelihoodSpec(d1, h1; num_int=BinnedDistributionFit.QuadGKIntegrator())
    @test NF1(ComponentVector(norms=[2.0], p1=0)) ≈ 1.6827899396467232

    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    h2 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF2 = BinnedDistributionFit.LikelihoodSpec(d2, h2; num_int=BinnedDistributionFit.QuadGKIntegrator())
    @test NF2(ComponentVector(norms=[2.0], p1=0)) ≈ 1.84583612533

    d3 = ExtendPdf((x, _)-> exp(x) - abs2(x) + 3, (1,3))
    h3 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF3 = BinnedDistributionFit.LikelihoodSpec(d3, h3; num_int=BinnedDistributionFit.QuadGKIntegrator())
    @test NF3(ComponentVector(norms=[2.0], p1=0)) ≈ 1.90019114214
end

@testitem "LikelihoodSpec SumOfPdfs" begin
    using FHist, ComponentArrays
    d1 = ExtendPdf((x, _) -> x, (1,3))
    d2 = ExtendPdf((x, _) -> x^2, (1,3))
    d3 = ExtendPdf((x, _) -> exp(x) - abs2(3x - 1) + 3, (1,3))
    h1 = Hist1D(; binedges = 1:3, bincounts = [2.0, 4.0], sumw2 = [2.0, 4.0])
    h2 = Hist1D(; binedges = 1:3, bincounts = [10.0, 1.0], sumw2 = [1.0, 100.0])
    NF1 = BinnedDistributionFit.LikelihoodSpec(d1 + d2, h1; num_int = BinnedDistributionFit.QuadGKIntegrator())
    NF2 = BinnedDistributionFit.LikelihoodSpec(d1 + d3, h1; num_int = BinnedDistributionFit.QuadGKIntegrator())
    NF3 = BinnedDistributionFit.LikelihoodSpec(d2 + d3, h1; num_int = BinnedDistributionFit.QuadGKIntegrator())

    @test NF1(ComponentVector(norms=[2.0,0.5], p1=0, p2=0)) ≈ 0.8497340452248006
    @test NF2(ComponentVector(norms=[2.0,0.5], p1=0, p2=0)) ≈ 0.850801892114
    @test NF3(ComponentVector(norms=[2.0,0.5], p1=0, p2=0)) ≈ 1.07158581388
end