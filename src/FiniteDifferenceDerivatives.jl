module FiniteDifferenceDerivatives

using ArrayViews

export fdd, fddcoeffs, fdd13, fdd_inline

function fdd{T<:Number}(der::Int,order::Int,x::AbstractVector{T},f::AbstractVector{T})
    if order < der
        error("Order can not be smaller than der")
    end

    npts = length(x)
    c    = zeros(T,order,der+1)
    df   = zero(f)

    for n = 1:npts
        i1 = int(min(max(1,(n-floor((order-1)/2))),npts-order+1))
        x0 = x[n]

        fddcoeffs(c,der,x0,view(x,i1:i1+order-1))

        for j = 1:order
            df[n] += c[j,der+1]*f[i1+j-1]
        end
    end
    return df
end

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


function fdd_inline{T<:Number}(der::Int,order::Int,x::AbstractVector{T},f::AbstractVector{T})
    if order < der
        error("Order can not be smaller than der")
    end

    npts = length(x)
    c    = zeros(T,order,der+1)
    df   = zero(f)

    for N = 1:npts
        i1 = int(min(max(1,(N-floor((order-1)/2))),npts-order+1))
        x0 = x[N]

        c1 = one(T)
        c4 = x[i1] - x0
        c[:] = zero(T)
        c[1] = one(T)
        for i=1:order-1
            mn = min(i,der)
            rng = mn:-1:1
            c2 = one(T)
            c5 = c4
            c4 = x[i+i1] - x0
            for j=0:i-1
                c3 = x[i+i1] - x[j+i1]
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
        for j = 1:order
            df[N] += c[j,der+1]*f[i1+j-1]
        end
    end
    return df
end


function fdd13{T<:Number}(x::AbstractVector{T},f::AbstractVector{T})
    # first derivative of order = 2
    npts=length(x)
    df = zero(f)
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

    return df
end

end # module
