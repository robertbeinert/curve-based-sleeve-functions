struct Curve 
    value
    tangent
    normal
    samples
    dim
end

function Curve(val, tan, nor; n::Signed=100)
    sam = [val(t / (n-1)) for t=0:n-1]
    dim = length(val(0))
    Curve(val, tan, nor, sam, dim)
end


function (γ::Curve)(t)
    γ.value(t)
end

function projection(x, γ::Curve; ϵ=1e-10, K=100)
    id = argmin([norm(x - γ.samples[n]) for n=1:length(γ.samples)])
    t = (id - 1) / (length(γ.samples) - 1)
    k = 0
    tUpdate = Inf
    tDerivative = Inf
    while k < K && abs(tUpdate) > ϵ && abs(tDerivative) > ϵ
        γt = γ.value(t)
        γpt = γ.tangent(t)
        γppt = γ.normal(t)
        tDerivative = -2dot(x - γt, γpt)
        tUpdate = tDerivative / (2norm(γpt)^2 - 2dot(x - γt, γppt))
        if (t ≤ ϵ && tDerivative ≥ 0.) || (t ≥ 1 - ϵ && tDerivative ≤ 0.)
            break    
        end
        t = min(max(t - tUpdate, 0), 1)
        k += 1
    end
    γ.value(t)
end

function dist(x, γ::Curve)
    y = projection(x, γ)
    norm(y - x)
end

struct Profile
    value
    derivative
end

struct Sleeve
    γ::Curve
    g::Profile
    h
end

function Sleeve(γ::Curve, g::Profile; h=0.0)
    Sleeve(γ, g, h)
end

function (f::Sleeve)(x)
    g.value(dist(x, f.γ)^2)
end

function ∇(f::Sleeve, x)
    if f.h == 0.0
        y = projection(x, f.γ)
        return 2f.g.derivative(norm(x - y)^2) .* (x - y)
    else
        return finiteDifferences(f, x)
    end
end

function finiteDifferences(f::Sleeve, x)
    d = f.γ.dim
    y = zeros(d)
    for k = 1:d
        e = zeros(d)
        e[k] = 1.
        y[k] = (f(x + f.h * e) - f(x - f.h * e)) / 2f.h
    end
    return y
end