function fdd13!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T})
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

function fdd23!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T})
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


function fdd1n!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T},order)
    # first derivative using a three point stencil
    npts=length(x)
    for i = 1:npts
        i1    = min(max(1,i-div(order-1,2)),npts-order+1)
        j     = i-i1+1          # i=1 => i1=1 => j=1; i=npts => i1=npts-order+1 => j = order
        df[i] = fdd1n_point(f,x,i1-1,j,order)
    end
end

function fdd1n_point{T}(f::AbstractVector{T},x::AbstractVector{T},shift,j,n)
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
