# compute the der-derivative using order-point scheme
function fdd_nonuniform!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T},der::Int,order::Int)
    npts = length(x)
    if order < der
        error("Order can not be smaller than der")
    elseif npts < order
        error("Number of mesh points shouldn't be smaller than order")
    end

    # specialized implementation for lower orders
    if der == 1 && order == 3
        fdd13_nonuniform!(df,f,x)
    elseif der == 2 && order == 3
        fdd23_nonuniform!(df,f,x)
    elseif der == 1
        fdd1n_nonuniform!(df,f,x,order)
    else
        # proceed with general implementation

        c = Array(T,order,der+1)

        for N = 1:npts
            # N1 is the leftmost point of the stencil
            N1    = min(max(1,N-div(order-1,2)),npts-order+1)
            df[N] = fdd_nonuniform_core(f,x,der,x[N],N1,order,c)
        end
    end

end


function fdd_nonuniform_core{T<:Number}(f::AbstractVector{T},
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

function fddat_nonuniform{T<:Number}(f::AbstractVector{T},
                                     x::AbstractVector{T},
                                     der::Int,
                                     x0::Real)

    if length(f) != length(x)
        error("x and f must have the same size")
    end

    N1 = 1
    order = length(x)
    c = Array(T,order,der+1)
    return fdd_nonuniform_core(f,x,der,x0,N1,order,c)
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

# specialized implementations

function fdd13_nonuniform!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T})
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

function fdd23_nonuniform!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T})
    # first derivative using a three point stencil
    npts=length(x)
    for i = 1:npts
        if i == 1
            df[i] = 2((f[i+2]-f[i+1])/(x[i+2]-x[i+1])-(f[i+1]-f[i])/(x[i+1]-x[i]))/(x[i+2]-x[i])
        elseif i == npts
            df[i] = 2((f[i]-f[i-1])/(x[i]-x[i-1])-(f[i-1]-f[i-2])/(x[i-1]-x[i-2]))/(x[i]-x[i-2])
        else
            df[i] = 2((f[i+1]-f[i])/(x[i+1]-x[i])-(f[i]-f[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])
        end
    end
end

function fdd1n_nonuniform!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T},order)
    # first derivative using an n-point stencil
    npts=length(x)
    for i = 1:npts
        i1    = min(max(1,i-div(order-1,2)),npts-order+1)
        j     = i-i1+1          # i=1 => i1=1 => j=1; i=npts => i1=npts-order+1 => j = order
        df[i] = fdd1n_point_nonuniform(f,x,i1-1,j,order)
    end
end

function fdd1n_point_nonuniform{T}(f::AbstractVector{T},x::AbstractVector{T},shift,j,n)
    df=zero(T)
    for i=1:n
        coef=one(T)
        if i!=j
            for k=1:n
                if k!=i
                    coef /= x[i+shift]-x[k+shift]
                end
                if k!=i && k!=j
                    coef *= x[j+shift]-x[k+shift]
                end
            end
            df += coef*(f[i+shift]-f[j+shift])
        end
    end
    return df
end
