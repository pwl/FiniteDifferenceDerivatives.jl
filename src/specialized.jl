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


# these functions don't work well on nonuniform meshes
function fdd15!(df,f,x)
    for i = 1:length(df)
        f_ξ = fdd15_uniform(f,i)
        x_ξ = fdd15_uniform(x,i)
        df[i] = f_ξ/x_ξ
    end
end

function fdd25!(df,f,x)
    for i = 1:length(df)
        x_ξ  = fdd15_uniform(x,i)
        x_ξξ = fdd25_uniform(x,i)
        f_ξ  = fdd15_uniform(f,i)
        f_ξξ = fdd25_uniform(f,i)
        df[i] = f_ξξ/x_ξ^2-x_ξξ*f_ξ/x_ξ^3
    end
end

function fdd15_uniform(f,i)
    nx = length(f)
    if i == 1
        df = (25*f[1]-48*f[2]+36*f[3]-16*f[4]+3*f[5])/12
    elseif i == 2
        df = (-3*f[1]-10*f[2]+18*f[3]-6*f[4]+f[5])/12
    elseif i == nx-1
        df = (f[nx-4]-6*f[nx-3]+18*f[nx-2]-10*f[nx-1]-3*f[nx])/12
    elseif i == nx
        df = (3*f[nx-4]-16*f[nx-3]+36*f[nx-2]-48*f[nx-1]+25*f[nx])/12
    else
        df = (f[i+2]-8f[i+1]+8f[i-1]-f[i-2])/12
    end
    return df
end

function fdd25_uniform(f,i)
    nx = length(f)
    if i == 1
        df = (35*f[i]-104*f[1+i]+114*f[2+i]-56*f[3+i]+11*f[4+i])/12
    elseif i == 2
        df = (11*f[-1+i]-20*f[i]+6*f[1+i]+4*f[2+i]-f[3+i])/12
    elseif i == nx-1
        df = -(f[-3+i]-4*f[-2+i]-6*f[-1+i]+20*f[i]-11*f[1+i])/12
    elseif i == nx
        df = (11*f[-4+i]-56*f[-3+i]+114*f[-2+i]-104*f[-1+i]+35*f[i])/12
    else
        df = -(f[-2+i]-16*f[-1+i]+30*f[i]-16*f[1+i]+f[2+i])/12
    end
    return df
end
