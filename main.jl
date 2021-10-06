# finite element analysis for 1D bar problem
# tuthor: @WuJC
# problem: EA*d²u/dx² = 0.,   x∈(0,1)
#          u(0) = 0.
#          EAdu/dx(1) = 1.
# exact solution: u(x) = x

# length of bar
Lb = 1.
# material coefficients
EA = 1.

# num of nodes
np = 11

# num of cells
nel = np - 1

# nodes
x = zeros(np)
for i in 1:nel
    x[i+1] = i*Lb/nel
end

# init stiffness matirx and force vector
k = zeros(np,np)
f = zeros(np)

# main loop
for i = 1:nel
    x1 = x[i]
    x2 = x[i+1]
    L = x2 - x1

    # vector N = [N₁,N₂] = [1-ξ,1+ξ]/2
    # ε = ∑Bᵢdᵢ
    # Bᵢ = dNᵢ/dξ * dξ/dx
    # In each cell, vector B = [B₁,B₂]
    # [dN₁/dξ,dN₂/dξ] = [-1,1]/2
    # dξ/dx = 2/L
    B = [-1,1]/L
    # local stiffness in each cell:ke = ∫BEAB'dx
    # or can be rewritten as ke = [1 -1;-1 1]/L
    ke = B*EA*B'*L
    # assemble into global stiffness matrix
    k[i:i+1,i:i+1] = k[i:i+1,i:i+1] + ke
end

# apply natural boundary
# EAdu/dx(1) = 1.
f[np] += 1.

# apply essential boundary
k[1,:] .= 0.
k[:,1] .= 0.
k[1,1] = 1.
f[1] = 0.

# solve coefficients
d = k\f
