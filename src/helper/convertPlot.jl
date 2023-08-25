function convertPlotToFile(file, x, y)
    if length(x) ≠ length(y)
        error("dimension mismatch")
    end
    io = open(file, "w")
    @printf(io, "%25s %25s\n", "x", "y")
    for k = 1 : length(x)
        @printf(io, "%25.6e %25.6e\n", x[k], y[k])
    end
    close(io)
end