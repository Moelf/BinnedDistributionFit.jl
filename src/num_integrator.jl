abstract type NumericalIntegrator end
struct QuadGKIntegrator <: NumericalIntegrator end
struct TrapzIntegrator <: NumericalIntegrator end
struct SimpleSumIntegrator <: NumericalIntegrator end

function _integrate(func, hist::Hist1D, ::QuadGKIntegrator)
    lo, hi = extrema(binedges(hist))
    return quadgk(func, lo, hi)[1]
end
function _integrate(func, hist::Hist1D, ::TrapzIntegrator)
    lo, hi = extrema(binedges(hist))
    sample_points = range(; start = lo, stop = hi, length = 10 * nbins(hist))
    return trapz(sample_points, func.(sample_points))
end
function _integrate(func, hist::Hist1D, ::SimpleSumIntegrator)
    bcs = bincenters(hist)
    return sum(func.(bcs))
end
ß