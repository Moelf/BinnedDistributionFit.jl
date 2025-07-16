# [Quick Start](@id quickstart)
```julia
using FHist, Optimize, ForwardDiff
# User inputs
d(x, ps...) = x^ps[1]+ps[2]                         # <-- Enter function here
h = Hist1D(; bincounts = 1:100, binedges = 1:100)   # <-- Enter histogram here

support = extrema(binedges(h1))

NLL = LikelihoodSpec(ExtendPdf(d, support), h1)

# Guess of the values of ps
ps = [1, 2]
# Optimization implementation
optf = OptimizationFunction(NLL_wrapper, AutoForwardDiff())
prob = OptimizationProblem(optf, ComponentArray(norms=integral(h),p1 = ps))
sol = solve(prob, Optimization.LBFGS())

@show sol.u
```
