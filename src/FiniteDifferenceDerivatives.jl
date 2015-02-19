module FiniteDifferenceDerivatives

using ArrayViews


include("uniform.jl")
include("nonuniform.jl")


export fdd!, fdd, fddmatrix, fddat


# returns a sparse matrix D such that D*u is a der-derivative of specified k-point scheme
function fddmatrix{T<:Number}(x::AbstractVector{T},der::Int,k::Int;args...)
    npts = length(x)
    m = zeros(T,npts,npts)
    f = eye(T,npts)
    for i = 1:npts
        m[:,i]=fdd(f[:,i],x,der,k;args...)
    end
    return sparse(m)
end


# compute the der-derivative using k-point scheme
function fdd!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T},der::Int,k::Int;uniform=true)
    if uniform
        fdd_uniform!(df,f,x,der,k)
    else
        fdd_nonuniform!(df,f,x,der,k)
    end
end


function fdd{T<:Number}(f::AbstractVector{T},x::AbstractVector{T},der::Int,k::Int;args...)
    df = Array(T,length(f))
    fdd!(df,f,x,der,k;args...)
    return df
end


end # module
