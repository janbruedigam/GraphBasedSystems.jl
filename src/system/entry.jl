mutable struct Entry{ET}
    value::ET

    function Entry{T}(dims...; static = true) where T
        static ? value = szeros(T,dims...) : value = zeros(T,dims...)
        new{typeof(value)}(value)
    end

    Entry(dims...; static = true) = Entry{Float64}(dims...; static = static)
end

function randomize!(entry::Entry)
    value = entry.value
    entry.value = randn(eltype(value), size(value))
end