module FiniteDifferenceDerivatives

using ArrayViews

export fdd

function fdd{T<:Number}(der::Int,order::Int,x::Vector{T},f::Vector{T})
    if order < der
        error("Order can not be smaller than der")
    end

    npts  = length(x)
    c = zeros(order,der+1)
    df = zero(f)

    for n = 1:npts
        i1 = min(max(n,1),npts-order+1)
        z = x[n]

        c[:]=zero(T)            # c(:,:) = 0
        c[1]=one(T)             # c(0,0) = 1
        c1 = one(T)             # c1 = 0
        c4 = x[i1]-z            # c4 = x(0)-z
        for i = 1:order-1       # i=1,n
            mn = min(i,der)     # mn = min(i,m)
            c2 = one(T)         # c2=1
            c5 = c4             # c5 = c4
            c4 = x[i1+i]-z      # c4 = x(i)-z
            for j = 0:i-1
                c3 = x[i1+i]-x[i1+j] # x(i)-x(j)
                c2 = c2*c3
                if j == i-1
                    for k = mn:-1:1
                        c[i+1,k+1] = c1*(k*c[i,k]-c5*c[i,k+1])/c2
                    end
                    c[i+1,1] = -c1*c5/c2*c[i,1]
                end
                for k = mn:-1:1
                    c[j+1,k+1] = (c4*c[j+1,k+1]-k*c[j+1,k])/c3
                end
                c[j+1,1]*=c4/c3*c[j+1,1]
            end
            c1 = c2
        end

        for j = 1:order
            df[n] += c[j,der+1]*f[i1+j-1]
        end
    end
    return df
end

end # module
