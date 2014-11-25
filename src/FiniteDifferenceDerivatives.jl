module FiniteDifferenceDerivatives

export fdd!, fdd, fddmatrix, fddat


# returns a sparse matrix D such that D*u is a der-derivative of specified order
function fddmatrix{T<:Number}(der::Int,order::Int,x::AbstractVector{T})
    npts = length(x)
    m = zeros(T,npts,npts)
    f = eye(T,npts)
    for i = 1:npts
        m[:,i]=fdd(der,order,x,f[:,i])
    end
    return sparse(m)
end


# compute the der-derivative using order-point scheme
function fdd!{T<:Number}(df::AbstractVector{T},der::Int,order::Int,x::AbstractVector{T},f::AbstractVector{T})
    npts = length(x)
    if order < der
        error("Order can not be smaller than der")
    elseif npts < order
        error("Number of mesh points shouldn't be smaller than order")
    end

    # specialized implementation for lower orders
    if der == 1 & order == 3
        fdd13!(df,x,f)
        return
    end

    c    = zeros(T,order,der+1)

    for N = 1:npts

        N1 = min(max(1,N-div(order-1,2)),npts-order+1)
        x0 = x[N]

        # essentially the following code is the contents of fddcoeffs
        # inlined almost without a change
        c1 = one(T)
        c4 = x[N1] - x0
        c[:] = zero(T)
        c[1] = one(T)
        for i=1:order-1
            mn = min(i,der)
            c2 = one(T)
            c5 = c4
            @inbounds c4 = x[i+N1] - x0
            j = zero(Int)
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

        dfN = zero(T)
        @simd for j = 1:order
            @inbounds dfN += c[j,der+1]*f[N1+j-1]
        end
        df[N] = dfN
    end
end


function fdd{T<:Number}(der::Int,order::Int,x::AbstractVector{T},f::AbstractVector{T})
    df = zero(f)
    fdd!(df,der,order,x,f)
    return df
end


function fdd13!{T<:Number}(df::AbstractVector{T},x::AbstractVector{T},f::AbstractVector{T})
    # first derivative using a three point stencil
    npts=length(x)
    for i = 1:npts
        if i == 1
            df2=f[i+1]-f[i]
            df3=f[i+2]-f[i]
            h2=x[i+1]-x[i]
            h3=x[i+2]-x[i]
            df[i]=df2/h2 + (df3-df2)/(h2-h3) + df3/h3
        elseif i == npts
            df1=f[i-1]-f[i]
            df2=f[i-2]-f[i]
            h1=x[i-1]-x[i]
            h2=x[i-2]-x[i]
            df[i]=df1/h1 + (df2-df1)/(h1-h2) + df2/h2
        else
            df1=f[i-1]-f[i]
            df3=f[i+1]-f[i]
            h1=x[i-1]-x[i]
            h3=x[i+1]-x[i]
            df[i]=df1/h1 + (df3-df1)/(h1-h3) + df3/h3
        end
    end
end


function fddat{T<:Number}(f::Vector{T},x::Vector{T},der::Int,x0::Real)
    order = length(x)

    if order != length(f)
        error("x and f should have the same length")
    end

    N1 = 1
    c    = zeros(T,order,der+1)
    generatec!(c,x0,x,N1)

    dfN = zero(T)
    @simd for j = 1:order
        @inbounds dfN += c[j,der+1]*f[N1+j-1]
    end
    return dfN
end


# generate the coefficients c
function generatec!{T}(c::Matrix{T},x0::Real,x::Vector{T},N1::Int)
    if N1 < 1 | N1 > length(x)
        error("N1 is out of bounds")
    end

    order = size(c,1)
    der   = size(c,2)-1

    c1 = one(T)
    c4 = x[N1] - x0
    c[1] = one(T)
    for i=1:order-1
        mn = min(i,der)
        c2 = one(T)
        c5 = c4
        @inbounds c4 = x[i+N1] - x0
        j = zero(Int)
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
