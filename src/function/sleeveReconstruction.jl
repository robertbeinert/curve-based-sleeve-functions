struct ApproxProfile
    knots
    samples
end

struct ApproxSleeve 
    γ
    g2::ApproxProfile
end

function sleeveApproximation(f::Sleeve, ρ, E, σ)
    x0 = zeros(f.γ.dim)
    g2 = initialProfile(f, ρ, σ, x0)
    γ = curveApproximation(f, g2, ρ, E)
    extendProfile!(g2, f, σ, ρ, γ)
    return ApproxSleeve(γ, g2)
end

function initialProfile(f::Sleeve, ρ, σ, x0)
    g2, k, l = sampleProfile(f, ρ, σ, x0)
    offset = determineProfileOffset(f, x0, (k-1) * σ) - (k-1) * σ
    if offset > 0
        knots = [0.0; offset .+ σ * (0 : k + l - 3)]
        g2 = [0.0; g2[2:end-2]]
    else
        knots = [0.0; offset .+ σ * (1 : k + l - 3)]
        g2 = [0.0; g2[3:end-2]]
    end
    ApproxProfile(knots, g2)
end

function sampleProfile(f::Sleeve, ρ, σ, x0)
    ∇fx0 = ∇(f, x0)
    dir = ∇fx0 ./ norm(∇fx0)
    k = 1
    g2 = [f(x0 - σ * dir), f(x0)]
    while g2[1] < g2[2]
        k += 1
        prepend!(g2, [f(x0 - k * σ * dir)])
    end
    l = 0
    while g2[end-1] < g2[end] && σ * (k + l) < ρ
        l += 1
        append!(g2, [f(x0 + l * σ * dir)])
    end
    return g2, k, l
end

function determineProfileOffset(f::Sleeve, x0, t0; ϵ=1e-10, K=100)
    ∇fx0 = ∇(f, x0)
    y = ∇fx0 ./ norm(∇fx0)
    t = t0
    k = 0
    tUpdate = Inf
    tDerivative = -dot(∇(f,x0 - t * y), y)
    while k < K && abs(tUpdate) > ϵ && abs(tDerivative) > ϵ
        tDerivative = -dot(∇(f,x0 - t * y), y)
        tUpdate = f(x0 - t * y) / tDerivative
        t = t - tUpdate
        k += 1
    end
    return t
end

function curveApproximation(f::Sleeve, g2::ApproxProfile, ρ, E)
    η = min(ρ, sqrt(ρ^2 - (ρ - E)^2))
    P0 = approxProjection(zeros(f.γ.dim), f, g2)
    P1 = generateSecondPoint(P0, f, g2, η)
    Pforward = sampleCurve(P0, P1, f, g2, ρ, η)
    Pbackward = sampleCurve(P1, P0, f, g2, ρ, η)
    return [reverse(Pbackward); Pforward]
end

function generateSecondPoint(P0, f::Sleeve, g2::ApproxProfile, η; ϵ=1e-6)
    P1 = P0
    while norm(P0 - P1) < ϵ
        v = randn(f.γ.dim)
        v *= (η / 2 / norm(v))
        P1 = approxProjection(P0 + v, f, g2)
    end
    return P1
end

function sampleCurve(P0, P1, f::Sleeve, g2::ApproxProfile, ρ, η; ϵ=1e-6)
    P = Vector{Vector{Float64}}(undef, 0)
    Plast = P0
    Pnew = P1
    while norm(Plast - Pnew) > ϵ
        append!(P, [Pnew])
        s = (η^2 + 2η * ρ) / (2ρ + 2η + norm(Plast - Pnew))
        v = Pnew - Plast
        v *= s / norm(v)
        Plast = Pnew
        Pnew = approxProjection(Plast + v, f, g2)
    end
    return P
end

function extendProfile!(g2::ApproxProfile, f::Sleeve, σ, ρ, γ)
    id = argmax(norm.(γ))
    x = γ[id]
    y = (1 + ρ / norm(x)) * x
    x = approxProjection(y, f, g2)
    dir = (y - x) / norm(y- x)
    t = g2.knots[end]
    while t < 1
        t += σ
        append!(g2.knots, [t])
        append!(g2.samples, [f(x + t * dir)])
    end
end

function approxProjection(x, f::Sleeve, g2; ϵ=1e-10)
    z = f(x)
    if z <= ϵ
        return x
    else
        v = ∇(f, x)
        return x - (inverseProfile(z, g2) / norm(v)) * v 
    end
end

function inverseProfile(z, g2::ApproxProfile)
    id = min(searchsortedfirst(g2.samples, z), length(g2.samples))
    if id == 1
        return 0.
    else
        return g2.knots[id-1] + (z - g2.samples[id-1]) / 
        (g2.samples[id] - g2.samples[id-1]) * (g2.knots[id] - g2.knots[id-1])
    end
end

function (f::ApproxSleeve)(x)
    f.g2(distToPolygonalChain(x, f.γ))
end

function (g2::ApproxProfile)(t)
    id = min(searchsortedfirst(g2.knots, t), length(g2.samples))
    if id == 1
        return g2.samples[1]
    else
        return g2.samples[id-1] + (g2.samples[id] - g2.samples[id-1]) /
        (g2.knots[id] - g2.knots[id-1]) * (t - g2.knots[id-1])
    end
end

function distToPolygonalChain(x, P)
    distances = zeros(length(P) - 1)
    for k = 1 : length(P) - 1
        y = x - P[k]
        v = P[k+1] - P[k]
        h = norm(v)
        v /= h
        λ = dot(y, v)
        if λ <= 0  
            distances[k] = norm(x - P[k])
        elseif λ >= h
            distances[k] = norm(x - P[k+1])
        else
            distances[k] = norm(y - λ * v)
        end 
    end
    minimum(distances)
end