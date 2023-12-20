using Pkg
using Symbolics
using LinearAlgebra
using StaticArrays
using ModelingToolkit
using SparseArrays 
using Serialization
using Plots 
using Dates 
#Pkg.instantiate()
gr() # grafics
# Number of derivatives 
nd = 4
# NUmber of variables
nv = 2
le=1 #to compile the symbolic expressions
#computing symbolic?
if (le==1)
include("blockA.jl")
else
    FAAb_expr1 = deserialize(open("jacobians/FAAb_expr1.jls"))
    FAAr_expr1 = deserialize(open("jacobians/FAAr_expr1.jls"))
    FAAl_expr1 = deserialize(open("jacobians/FAAl_expr1.jls"))
    dFAAb_expr1 = deserialize(open("jacobians/dFAAb_expr1.jls"))
    dFAAr_expr1 = deserialize(open("jacobians/dFAAr_expr1.jls"))
    dFAAl_expr1 = deserialize(open("jacobians/dFAAl_expr1.jls"))
end   

#numerical part
alpha10=0
nu0=0.1



# Spatial discretization
ns = 20000
H = 1.0
dt=10^10
# Collocation matrices (spatial discretization using 2d order finites differences)
include("collocation/finites2th.jl")
s0, ds0, dss0 = finites2th(ns, H)
#initial time
t=0.
# Initial conditions
u0 = 1 .- 2 * s0
g0 = s0
duds0 = ds0 * u0
dudss0 = dss0 * u0
dudtau0 = zeros(ns)
dgds0 = ds0 * g0
dgdss0 = dss0 * g0
dgdtau0 = zeros(ns)

#intial guess solution
x0=[u0;g0] #vector [2N,1];
x0m=x0
x0mm=x0



# Newton method (step 10)
 errorNewton = 1e9
 iter = 0


 while errorNewton  >= 1e-4
 global iter 
 global x0
 global errorNewton
  
  #matrix a and vector b: a*dx=b
  include("matrixAB.jl")
 
 
 # inverting the Jacobian
 dxa = a \ b
 # computing error
 errorNewton = norm(dxa)
 println("error: $errorNewton")
 iter=iter+1 
 x0=x0+dxa
 end



   # Making variables readable
   u = x0[1:ns]'
   g = x0[ns+1:2*ns]'

   #plotting results
   plot(g',u',color=:red,linestyle=:dash,xlabel="x",ylabel="u",title="Velocity")