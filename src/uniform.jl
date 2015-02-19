# these functions work only on uniform meshes

function fdd_uniform!{T<:Number}(df::AbstractVector{T},f::AbstractVector{T},x::AbstractVector{T},der::Int,k::Int)

    h = x[2]-x[1]

    # specialized implementation for lower ks
    for i = 1:length(f)
        if der == 1 && k == 5
            df[i] = fdd25_uniform(f,i)
        elseif der == 2 && k == 5
            df[i] = fdd25_uniform(f,i)
        elseif der == 1 && k == 3
            df[i] = fdd13_uniform(f,i)
        elseif der == 2 && k == 3
            df[i] = fdd23_uniform(f,i)
        else
            error("No generic method implemented for uniform mesh")
        end
        df[i] /= h^der
    end
end

function fdd13_uniform(f,i)
    nx = length(f)
    if i == 1
        df = (f[2]-f[1])
    elseif i == nx
        df = (f[nx]-f[nx-1])
    else
        df = (f[i+1]-f[i-1])/2
    end
    return df
end

function fdd23_uniform(f,i)
    nx = length(f)
    if i == 1
        df = (f[1]-2*f[2]+f[3])
    elseif i == nx
        df = (f[nx-2]-2*f[nx-1]+f[nx])
    else
        df = (f[i-1]-2*f[i]+f[i+1])
    end
    return df
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
