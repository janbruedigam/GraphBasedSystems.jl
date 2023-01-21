mutable struct Entry{ET}
    value::ET
    isinverted::Bool # for solver functions


    Entry(value::ET, isinverted::Bool) where ET = new{ET}(value, isinverted)

    function Entry{T}(dims...; static = true) where T
        static ? value = szeros(T, dims...) : value = zeros(T, dims...)
        Entry(value, false)
    end
    Entry(dims...; static = true) = Entry{Float64}(dims...; static)
    
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, entry::Entry)
    println(io, "Entry with value:")
    show(io, mime, entry.value)
end


Base.zero(::Entry{ET}) where ET = Entry(zero(ET), false)
Base.zero(::Type{Entry}) = Entry(nothing, false)

Base.:+(E1::Entry, E2::Entry) = Entry(E1.value + E2.value, false)
Base.:-(E1::Entry, E2::Entry) = Entry(E1.value - E2.value, false)
Base.:*(E1::Entry, E2::Entry) = Entry(E1.value * E2.value, false)
Base.:*(E::Entry, r::Real) = Entry(E.value * r, false)
Base.:*(r::Real, E::Entry) = E*r
Base.:/(E1::Entry, E2::Entry) = Entry(E1.value / E2.value, false)
Base.:\(E1::Entry, E2::Entry) = Entry(E1.value \ E2.value, false)
Base.:/(E::Entry, r::Real) = Entry(E.value / r, false)
Base.:\(r::Real, E::Entry) = E/r
Base.inv(E::Entry) = Entry(inv(E.value), false)
LinearAlgebra.pinv(E::Entry) = Entry(pinv(E.value), false)
LinearAlgebra.conj(E::Entry) = Entry(conj(E.value), false)


function initialize!(entry::Entry, init_function = randn)
    value = entry.value
    entry.value = init_function(eltype(value), size(value))
end