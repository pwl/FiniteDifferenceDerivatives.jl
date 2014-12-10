module FiniteDifferenceDerivatives

using ArrayViews

include("specialized.jl")

export fdd!, fdd, fddmatrix, fddat


# returns a sparse matrix D such that D*u is a der-derivative of specified order
function fddmatrix{T<:Number}(x::AbstractVector{T},der::Int,order::Int)
    npts = length(x)
    m = zeros(T,npts,npts)
    f = eye(T,npts)
    for i = 1:npts
        m[:,i]=fdd(f[:,i],x,der,order)
    end
    return sparse(m)
end


# compute the der-derivative using order-point scheme
function fdd!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T},der::Int,order::Int)
    npts = length(x)
    if order < der
        error("Order can not be smaller than der")
    elseif npts < order
        error("Number of mesh points shouldn't be smaller than order")
    end

    # specialized implementation for lower orders
    if der == 1 && order == 3
        fdd13!(df,f,x)
    elseif der == 2 && order == 3
        fdd23!(df,f,x)
    elseif der == 1
        fdd1n!(df,f,x,order)
    else
        # proceed with general implementation

        c = Array(T,order,der+1)

        for N = 1:npts
            # N1 is the leftmost point of the stencil
            N1    = min(max(1,N-div(order-1,2)),npts-order+1)
            df[N] = fddcore(f,x,der,x[N],N1,order,c)
        end
    end

end


function fdd{T<:Number}(f::AbstractVector{T},x::AbstractVector{T},der::Int,order::Int)
    df = Array(T,length(f))
    fdd!(df,f,x,der,order)
    return df
end


function fddcore{T<:Number}(f::AbstractVector{T},
                            x::AbstractVector{T},
                            der::Int,
                            x0::Real,
                            N1::Int,
                            order::Int,
                            c::Matrix{T})

    generatec!(c,x0,x,N1)

    dfN = zero(T)
    @simd for j = 1:order
        @inbounds dfN += c[j,der+1]*f[N1-1+j]
    end
    return dfN
end

function fddat{T<:Number}(f::AbstractVector{T},
                          x::AbstractVector{T},
                          der::Int,
                          x0::Real)

    if length(f) != length(x)
        error("x and f must have the same size")
    end

    N1 = 1
    order = length(x)
    c = Array(T,order,der+1)
    return fddcore(f,x,der,x0,N1,order,c)
end


# generate the coefficients c
function generatec!{T}(c::Matrix{T},x0::Real,x::AbstractVector{T},N1::Int)
    order = size(c,1)
    der   = size(c,2)-1

    if N1 < 1 || N1+order-1 > length(x)
        error("N1=$N1 and order=$order are out of bounds")
    end

    c1 = one(T)
    @inbounds c4 = x[N1] - x0
    c[:] = zero(T)
    c[1] = one(T)
    for i=1:order-1
        mn = min(i,der)
        c2 = one(T)
        c5 = c4
        @inbounds c4 = x[i+N1] - x0
        j = 0
        while j <= i-1
            @inbounds c3 = x[i+N1] - x[j+N1]
            c2 = c2*c3
            if j == i-1
                s = mn
                while s >= 1
                    @inbounds c[i+1,s+1] = c1*(s*c[i,s] - c5*c[i,s+1])/c2
                    s-=1
                end
                @inbounds c[i+1,1] = -c1*c5*c[i,1]/c2
            end
            s = mn
            while s >= 1
                @inbounds c[j+1,s+1] = (c4*c[j+1,s+1] - s*c[j+1,s])/c3
                s-=1
            end
            @inbounds c[j+1,1] = c4*c[j+1,1]/c3
            j+=1
        end
        c1 = c2
    end
end


end # module
