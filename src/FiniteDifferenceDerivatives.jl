module FiniteDifferenceDerivatives

export fdd!, fddcoeffs, fdd13!

function fddcoeffs{T<:Number}(c::Matrix{T},k::Int,x0::T,x::AbstractVector{T})
    n = length(x)

    if k >= n
        error("length(x) must be larger than k")
    end

    if size(c,1) < n || size(c,2) < k+1
        error("The size of c has to be at least $((n,k+1))")
    end

    c1 = one(T)
    c4 = x[1] - x0
    c[:] = zero(T)
    c[1] = one(T)
    for i=1:n-1
        mn = min(i,k)
        rng = mn:-1:1
        c2 = one(T)
        c5 = c4
        c4 = x[i+1] - x0
        for j=0:i-1
            c3 = x[i+1] - x[j+1]
            c2 = c2*c3
            if j==i-1
                for s=rng
                    c[i+1,s+1] = c1*(s*c[i,s] - c5*c[i,s+1])/c2
                end
                c[i+1,1] = -c1*c5*c[i,1]/c2
            end
            for s=rng
                c[j+1,s+1] = (c4*c[j+1,s+1] - s*c[j+1,s])/c3
            end
            c[j+1,1] = c4*c[j+1,1]/c3
        end
        c1 = c2
    end
end


function fdd!{T<:Number}(df::AbstractVector{T},der::Int,order::Int,x::AbstractVector{T},f::AbstractVector{T})
    if order < der
        error("Order can not be smaller than der")
    end

    npts = length(x)
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
    return df
end


function fdd13!{T<:Number}(df::AbstractVector{T},x::AbstractVector{T},f::AbstractVector{T})
    # first derivative of order = 2
    npts=length(x)
    for i = 1:npts
        df[i] = zero(T)
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

    return df
end

end # module
