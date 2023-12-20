# Second-order backward differences
dta = 1 / (2 * dt)
bm = -4 * dta
bmm = 1 * dta
bp = 3 * dta  # [STEP 6]

# Temporal derivatives
x0t = bp .* x0 .+ bm .* x0m .+ bmm .* x0mm


u = x0[1:ns]
duds = ds0 * x0[1:ns]
dudss = dss0 * x0[1:ns]
dudt = x0t[1:ns]
g = x0[ns+1:2*ns]
dgds = ds0* x0[ns+1:2*ns]
dgdss = dss0 * x0[ns+1:2*ns]
dgdt = x0t[ns+1:2*ns]



# Number of variables
nv = 2
# Number of symbolic derivatives
nd = 4

FAA = zeros(nv, ns)
DFAA = zeros( nv, nv*nd,ns)

# Getting analytical expressions [STEP8]
for i in 1:ns
    xs0 = [t,s0[i],u[i], duds[i], dudss[i], dudt[i], g[i], dgds[i], dgdss[i], dgdt[i],nu0,alpha10]  # Symbolic derivatives evaluation..
    # Bulk
    if i in 2:(ns-1)
        FAA[:, i] = FAAb_expr1(xs0)
        DFAA[:,:, i]= reshape(dFAAb_expr1(xs0),nv,nv*nd,1)
    end
    # Left
    if i == 1
        FAA[:, i] = FAAl_expr1(xs0) 
        c=reshape(dFAAl_expr1(xs0),nv,nv*nd,1)
        DFAA[:,:, i].= reshape(dFAAl_expr1(xs0),nv,nv*nd,1)
    end
    # Right
    if i == ns
        FAA[:, i] = FAAr_expr1(xs0)
        DFAA[:, :, i]=reshape(dFAAr_expr1(xs0),nv,nv*nd,1)
    end
end


# Getting numerical Jacobian (step 9)
ablock = [spzeros(ns, ns) for _ in 1:nv, _ in 1:nv]
xablock = Vector{SparseVector{Float64, Int64}}(undef, nv)
xablock = Array{SparseVector{Float64,Int64},1}(undef, nv)

for j in 1:nv
    C1 = FAA[j, :]
    xablock[j] = -C1
    
    for k in 1:nv
        km = (k - 1) * nd + 1
        kp = k * nd
        C = DFAA[j, km:kp, :]
       
        c0 = C[1, :]
        c1 = C[2, :]
        c2 = C[3, :]
        c3 = C[4, :]
        B = spdiagm(0 => c0) + spdiagm(0 => c1) * ds0 + spdiagm(0 => c2) * dss0 + bp * spdiagm(0 => c3)
        ablock[j, k] = B
    end
end
# Mouting BLOCK 


b=vcat(xablock...)
b=Array(b)
# Convert lists of sparse matrices into sparse matrices
a = reduce(vcat, [reduce(hcat, ablock[i, :]) for i in 1:nv])