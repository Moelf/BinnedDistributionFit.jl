using TestItems, TestItemRunner
@run_package_tests verbose=true

@testsnippet hist_dep begin
    using FHist
end
@testitem "ExtendPdf" begin
    d1 = ExtendPdf((x, _)->x, (1,3))
    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    d3 = ExtendPdf((x, _) -> xor((x << 2), x) % (x^2 + 1), (1,3))
    
    @test BinnedDistributionFit.scalar_eval(d1, 3) == 3
    @test BinnedDistributionFit.vector_eval(d1, [3, 5]) == [3, 5]

    @test BinnedDistributionFit.scalar_eval(d2, 3) == 9
    @test BinnedDistributionFit.vector_eval(d2, [3, 5]) == [9, 25]

    @test BinnedDistributionFit.scalar_eval(d3, 3) == 5
    @test BinnedDistributionFit.vector_eval(d3, [3, 5]) == [5, 17]
    end

@testitem "SumOfPdfs" begin
    d1 = ExtendPdf((x, _)->x, (1,3))
    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    d3 = ExtendPdf((x, _) -> xor((x << 2), x) % (x^2 + 1), (1,3))

    @test BinnedDistributionFit.scalar_eval(d1+d2, 3) == 12
    @test BinnedDistributionFit.vector_eval(d1+d2, [3, 5]) == [12, 30]

    @test BinnedDistributionFit.scalar_eval(d1+d3, 3) == 8
    @test BinnedDistributionFit.vector_eval(d1+d3, [3, 5]) == [8, 22]

    @test BinnedDistributionFit.scalar_eval(d2+d3, 3) == 14
    @test BinnedDistributionFit.vector_eval(d2+d3, [3, 5]) == [14, 42]
    end

@testitem "RooFitNLL ExtendPdf" setup=[hist_dep] begin
    # https://www.desmos.com/calculator/2e28c86a6e finds manual solutions to RooFitNLL
    d1 = ExtendPdf((x, _)->x, (1,3))
    h1 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF1 = RooFitNLL_functor(d1, h1; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF1([2.0]) ≈ 1.6827899396467232

    d2 = ExtendPdf((x, _)->abs2(x), (1,3))
    h2 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF2 = RooFitNLL_functor(d2, h2; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF2([2.0]) ≈ 1.84583612533

    d3 = ExtendPdf((x, _)-> exp(x) - abs2(x) + 3, (1,3))
    h3 = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF3 = RooFitNLL_functor(d3, h3; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF3([2.0]) ≈ 1.90019114214
end

@testitem "RooFit SumOfPdfs" setup=[hist_dep] begin
    # Using the same Desmos as before (https://www.desmos.com/calculator/2e28c86a6e) to generate solutions
    d1 = ExtendPdf((x, _) -> x, (1,3))
    d2 = ExtendPdf((x, _) -> x^2, (1,3))
    d3 = ExtendPdf((x, _) -> exp(x) - abs2(3x - 1) + 3, (1,3))
    h1 = Hist1D(; binedges = 1:3, bincounts = [2.0, 4.0], sumw2 = [2.0, 4.0])
    h2 = Hist1D(; binedges = 1:11, bincounts = [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 5.0, 8.0, 9.0, 10.0])
    NF1 = RooFitNLL_functor(d1 + d2, h1; num_integrator = BinnedDistributionFit.QuadGKIntegrator())
    NF2 = RooFitNLL_functor(d1 + d3, h1; num_integrator = BinnedDistributionFit.QuadGKIntegrator())
    NF3 = RooFitNLL_functor(d2 + d3, h1; num_integrator = BinnedDistributionFit.QuadGKIntegrator())

    @test NF1([[2.0, 0.5], [], []]) ≈ 0.8497340452248006
    @test NF2(([2.0, 0.5], [], [])) ≈ 0.850801892114
    @test NF3(([2.0, 0.5], [], [])) ≈ 1.07158581388
end
@testitem "chi2" begin
    o1 = [1, 2, 3]
    o2 = [3, 2, 1]
    o3 = [10, 10, 10]
    e1 = [1, 2, 3]
    e2 = [30, 20, 10]
    e3 = [5, 5, 5]
    sd1 = [1, 2, 3]
    sd2 = [0.01, 0.01, 0.01]
    sd3 = [1, 10, 0.1]
    @test BinnedDistributionFit.chi2(o1, e1, sd1) == 0
    @test BinnedDistributionFit.chi2(o2, e1, sd1) ≈ 4.4444444444444
    @test BinnedDistributionFit.chi2(o3, e1, sd1) ≈ 102.44444444444
    @test BinnedDistributionFit.chi2(o1, e2, sd1) ≈ 927.44444444444 
    @test BinnedDistributionFit.chi2(o2, e2, sd1) == 819 
    @test BinnedDistributionFit.chi2(o3, e2, sd1) == 425
    @test BinnedDistributionFit.chi2(o1, e3, sd1) ≈ 18.694444444444
    @test BinnedDistributionFit.chi2(o2, e3, sd1) ≈ 8.0277777777777
    @test BinnedDistributionFit.chi2(o3, e3, sd1) ≈ 34.027777777777
    @test BinnedDistributionFit.chi2(o1, e1, sd2) == 0
    @test BinnedDistributionFit.chi2(o2, e1, sd2) == 80000
    @test BinnedDistributionFit.chi2(o3, e1, sd2) == 1.94e6
    @test BinnedDistributionFit.chi2(o1, e2, sd2) == 1.214e7
    @test BinnedDistributionFit.chi2(o2, e2, sd2) == 1.134e7
    @test BinnedDistributionFit.chi2(o3, e2, sd2) == 5e6
    @test BinnedDistributionFit.chi2(o1, e3, sd2) == 290000
    @test BinnedDistributionFit.chi2(o2, e3, sd2) == 290000
    @test BinnedDistributionFit.chi2(o3, e3, sd2) == 750000
    @test BinnedDistributionFit.chi2(o1, e1, sd3) == 0
    @test BinnedDistributionFit.chi2(o2, e1, sd3) == 404
    @test BinnedDistributionFit.chi2(o3, e1, sd3) == 4981.64
    @test BinnedDistributionFit.chi2(o1, e2, sd3) == 5744.24
    @test BinnedDistributionFit.chi2(o2, e2, sd3) == 8832.24
    @test BinnedDistributionFit.chi2(o3, e2, sd3) == 401
    @test BinnedDistributionFit.chi2(o1, e3, sd3) == 416.09
    @test BinnedDistributionFit.chi2(o2, e3, sd3) == 1604.09
    @test BinnedDistributionFit.chi2(o3, e3, sd3) == 2525.25
end

@testitem "Chi2 ExtendPdf" setup=[hist_dep] begin
    # https://www.desmos.com/calculator/mbj4cblgzr finds manual solutions to BinnedDistributionFit.chi2
    d1 = ExtendPdf((x, _) -> x-0.5, (1,3))
    h1 = Hist1D(; binedges = 1:3, bincounts = [1, 2], sumw2 = [1, 1])
    x1 = BinnedDistributionFit.Chi2_functor(d1, h1; num_integrator = BinnedDistributionFit.QuadGKIntegrator())
    @test x1([3]) == 0

    d2 = ExtendPdf((x, _) -> (x^2 - 3x + 1), (1,5))
    h2 = Hist1D(; binedges = 1:5, bincounts = [1, 2, 3, 4], sumw2 = [6, 2, 1, 5])
    x2 = BinnedDistributionFit.Chi2_functor(d2, h2; num_integrator = BinnedDistributionFit.QuadGKIntegrator())
    @test x2([2]) ≈ 9.2824829932

    d3 = ExtendPdf((x, _) -> (2^x - 3x^3 + 4 - 1/x), (1,11))
    h3 = Hist1D(; binedges = 1:11, bincounts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], sumw2 = [6, 2, 1, 5, 10, 0.1, 3, 6, 9, 3.2])
    x3 = BinnedDistributionFit.Chi2_functor(d3, h3; num_integrator = BinnedDistributionFit.QuadGKIntegrator())
    @test x3([2]) ≈ 415.959048958
end

@testitem "Chi2 SumOfPdfs" setup=[hist_dep] begin
    # https://www.desmos.com/calculator/mbj4cblgzr finds manual solutions to BinnedDistributionFit.chi2
    d1 = ExtendPdf((x, _) -> x-0.5, (1, 4))
    d2 = ExtendPdf((x, _) -> 3x^2, (1, 4))
    d3 = ExtendPdf((x, _) -> exp(2x)+sin(5x-3x^2), (1, 11))
    d4 = ExtendPdf((x, _) -> 3x^3 - 4x + 7, (1, 11))
    h1 = Hist1D(; binedges = 1:4, bincounts = [1, 2, 3], sumw2 = [1, 2, 1])
    h2 = Hist1D(; binedges = 1:11, bincounts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], sumw2 = [5, 4, 3, 2, 1, 2, 3, 4, 5, 10])
    x1 = BinnedDistributionFit.Chi2_functor(d1 + d2, h1; num_integrator = BinnedDistributionFit.QuadGKIntegrator())
    x2 = BinnedDistributionFit.Chi2_functor(d3 + d4, h2; num_integrator = BinnedDistributionFit.QuadGKIntegrator())
    @test x1([[3, 4], [], []]) ≈ 1.73774816049
    @test x2([[5, 2.5], [], []]) ≈ 105.138115485
end

include("optimization.jl")

@testitem "_" begin
    using Aqua
    Aqua.test_all(BinnedDistributionFit)
end
