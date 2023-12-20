#simbolic part++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Symbolic variables
@variables s tau u duds dudss dudtau g dgds dgdss dgdtau nu alpha1

# 1. Mapping s to x
dudx = duds / dgds
dudxx = (dudss / dgds - duds * dgdss / dgds^2) / dgds
dgdx = dgds / dgds
dgdxx = (dgdss / dgds - dgds * dgdss / dgds^2) / dgds
dudt = dudtau - dgdtau / dgds * duds
dgdt = dgdtau - dgdtau / dgds * duds
M=1/(alpha1*dudx^2+0.4);
Ms=Symbolics.derivative(M, duds )*dudss+Symbolics.derivative(M, dgds )*dgdss

n=(dgds^2)^1.5
# Step 2: Mounting vectors and Jacobian tensors
xs = [u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau]
# Step 1: Equation for bulk
FAAb = [dudt + u * dudx - nu * dudxx,
dgds*dgdss-n*Ms]
# Boundary conditions
FAAl = [u - 1, g]    #left
FAAr = [u + 1, g - 1] #right
#analical Jacobians
dFAAb=Symbolics.jacobian(FAAb, xs)
dFAAl=Symbolics.jacobian(FAAl, xs)
dFAAr=Symbolics.jacobian(FAAr, xs)

#all the simbolic variables
xsa=[s,tau,u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau,nu,alpha1]
#Bulding functions from symbolics
FAAl_expr1 = eval(build_function(FAAl, xsa,expression=Val{false})[1])
FAAb_expr1 = eval(build_function(FAAb, xsa,expression=Val{false})[1])
FAAr_expr1 = eval(build_function(FAAr, xsa,expression=Val{false})[1])
dFAAl_expr1 = eval(build_function(dFAAl, xsa,expression=Val{false})[1])
dFAAb_expr1 = eval(build_function(dFAAb, xsa,expression=Val{false})[1])
dFAAr_expr1 = eval(build_function(dFAAr, xsa,expression=Val{false})[1])


# Guardar la función en un archivo
open("jacobians/FAAb_expr1.jls", "w") do f
    serialize(f, FAAb_expr1)
end
open("jacobians/FAAl_expr1.jls", "w") do f
    serialize(f, FAAl_expr1)
end
open("jacobians/FAAr_expr1.jls", "w") do f
    serialize(f, FAAr_expr1)
end

open("jacobians/dFAAb_expr1.jls", "w") do f
    serialize(f, dFAAb_expr1)
end
open("jacobians/dFAAl_expr1.jls", "w") do f
    serialize(f, dFAAl_expr1)
end
open("jacobians/dFAAr_expr1.jls", "w") do f
    serialize(f, dFAAr_expr1)
end


#New faster way (bulding+evaluating using Sybolic)
#FAAl_expr2 = eval(Symbolics.build_function(FAAl,xsa,
#            parallel=Symbolics.MultithreadedForm())[1])
#FAAb_expr2 = eval(Symbolics.build_function(FAAb,xsa,
#             parallel=Symbolics.MultithreadedForm())[1])
#FAAr_expr2 = eval(Symbolics.build_function(FAAr,xsa,
#             parallel=Symbolics.MultithreadedForm())[1])

#dFAAl_expr2 = eval(Symbolics.build_function(dFAAl,xsa,
#             parallel=Symbolics.MultithreadedForm())[1])
#dFAAb_expr2 = eval(Symbolics.build_function(dFAAb,xsa,
#              parallel=Symbolics.MultithreadedForm())[1])
#dFAAr_expr2 = eval(Symbolics.build_function(dFAAr,xsa,
#              parallel=Symbolics.MultithreadedForm())[1])