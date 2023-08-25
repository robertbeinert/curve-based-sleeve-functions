include("../../startup.jl")

a = 4π
b = 5π
c(t) = [
        (1 + 0.5cos(a * t)) * cos(b * t),
        (1 + 0.5cos(a * t)) * sin(b * t),
        sin(a * t)
]
cp(t) = [
        -0.5a * sin(a * t) * cos(b * t) - b * (1 + 0.5cos(a * t)) * sin(b * t),
        -0.5a * sin(a * t) * sin(b * t) + b * (1 + 0.5cos(a * t)) * cos(b * t),
        a * cos(a * t)
]
cpp(t) = [
        -0.5a^2 * cos(a * t) * cos(b * t) + a * b * sin(a * t) * sin(b * t) - 
        b^2 * (1 + 0.5cos(a * t)) * cos(b * t),
        -0.5a^2 * cos(a * t) * sin(b * t) - a * b * sin(a * t) * cos(b * t) - 
        b^2 * (1 + 0.5cos(a * t)) * sin(b * t),
        -a^2 * sin(a * t)
]

γ = Curve(c, cp, cpp)
#g = Profile(x -> x, x -> 1)
g = Profile(x -> tan(1.5x), x -> 1.5 * (1 + tan(1.5x)^2))
#g = Profile(x -> sin(0.5π * x), x -> 0.5π * cos(0.5π * x))
f = Sleeve(γ, g)


fRec = sleeveApproximation(f, 0.125, 1e-2, 1e-4)

display(fRec.γ)

x = [fRec.γ[k][1] for k=1:length(fRec.γ)]
y = [fRec.γ[k][2] for k=1:length(fRec.γ)]
z = [fRec.γ[k][3] for k=1:length(fRec.γ)]

# plotlyjs()
# display(plot(fRec.g2.knots.^2, fRec.g2.samples))
# png("sim/space/profile3d")
# display(plot(x, y, z, ls=:solid, shape=:auto, ms=2))
# png("sim/space/curve3d")


pgfplotsx()
display(plot(fRec.g2.knots.^2, fRec.g2.samples, label="", tex_output_standalone=true))
savefig("sim/space/profile3d.tex")
display(plot(x, y, z, ls=:solid, label="", tex_output_standalone=true))
savefig("sim/space/curve3d.tex")
