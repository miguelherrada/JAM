import sympy as sp
import pickle
# Define symbolic variables
w, dwds, dwdss, dwdy, dwdyy, dwdys, dwdt0 = sp.symbols('w dwds dwdss dwdy dwdyy dwdys dwdt0')
u, duds, dudss, dudy, dudyy, dudys, dudt0 = sp.symbols('u duds dudss dudy dudyy dudys dudt0')
p, dpds, dpdss, dpdy, dpdyy, dpdys, dpdt0 = sp.symbols('p dpds dpdss dpdy dpdyy dpdys dpdt0')
f, dfds, dfdss, dfdy, dfdyy, dfdys, dfdt0 = sp.symbols('f dfds dfdss dfdy dfdyy dfdys dfdt0')
# Define space (s,y) and time (tau)
s, y, t0 = sp.symbols('s y t0')
# Define parameters
Cmu, Bo = sp.symbols('Cmu  Bo')

# Mapping to (r,z)
# radial position
r = y * f
# axial position
z = s
# mapping factors
fs = dfds
fss = dfdss
mi = 1. / f
mis = -fs / f**2
miss = -fss / f**2 + 2 * fs**2 / f**3
r1 = r * mis
r2 = r1**2
r3 = r * miss
# Radial derivatives
dudr = dudy * mi
dpdr = dpdy * mi
dwdr = dwdy * mi
dudrr = dudyy * mi**2
dwdrr = dwdyy * mi**2
# Axial derivatives      
dwdz = dwds + r1 * dwdy
dudz = duds + r1 * dudy
dpdz = dpds + r1 * dpdy
dwdzz = dwdss + 2 * r1 * dwdys + r2 * dwdyy + r3 * dwdy
dudzz = dudss + 2 * r1 * dudys + r2 * dudyy + r3 * dudy
# temporal mapping 
ft = dfdt0
tt = -ft * r / f**2
dudt = dudt0 + tt * dudy     
dwdt = dwdt0 + tt * dwdy   

nv = 4 # number of equations
# tangential and normal vector+curvature
n = 1 + fs**2
# tangent vector
t0z = 1. /n**0.5
t0r = fs /n**0.5     
# normal vector
n0z = t0r
n0r = -t0z     
# stress tensor
taurr = 2 * dudr
tauzr = dwdr + dudz
tauzz = 2 * dwdz
# total stresses
# normal
tauni = n0r * taurr * n0r + n0z * tauzz * n0z + 2 * n0z * tauzr * n0r 
# tangential
tauti = t0r * taurr * n0r + t0z * tauzz * n0z + t0z * tauzr * n0r + n0z * tauzr * t0r 
# Interface curvature
k = (f * fss - 1 - fs**2) / (f * (1 + fs**2)**1.5)

# Create the symbolic vector using list concatenation
xs = [w, dwds, dwdss, dwdy, dwdyy, dwdys, dwdt0] 
xs = xs + [u, duds, dudss, dudy, dudyy, dudys, dudt0]
xs = xs + [p, dpds, dpdss, dpdy, dpdyy, dpdys, dpdt0]
xs = xs + [f, dfds, dfdss, dfdy, dfdyy, dfdys, dfdt0]
#Equations.......
# Define   bulk
FAb = [
sp.Eq(-(dwdt+u*dwdr+w*dwdz+dpdz)+(dwdzz+dwdrr+dwdr/r)*Cmu,0),
sp.Eq(-(dudt+u*dudr+w*dudz+dpdr)+(dudzz+dudrr+dudr/r-u/r**2)*Cmu,0),
sp.Eq(dwdz+dudr+u/r,0),
sp.Eq(dfdyy,0)
]
# Define   left
FAl = [
sp.Eq(w,0),
sp.Eq(u,0),
sp.Eq(dpdz,0),
sp.Eq(f-1,0)
]
# Define   right
FAr = [
sp.Eq(w,0),
sp.Eq(u,0),
sp.Eq(dpdz,0),
sp.Eq(f-1,0)
]
# Define   top (interface)
FAt = [
sp.Eq(p+Bo*z+k-tauni*Cmu,0), #normal balance
sp.Eq(Cmu*tauti,0), #tangential balance
sp.Eq(dpdr,0), #pressure
sp.Eq(ft+w*fs-u,0) #Kinematic equation
]

# Define   bottom (axis)
FAd = [
sp.Eq(dwdr,0), 
sp.Eq(u,0),
sp.Eq(dpdr,0), 
sp.Eq(dfdy,0) 
]
#Analytical  Jacobians
# Calculate the Jacobian of the equation with respect to the vector of variables.
DFAl = sp.Matrix([[sp.diff(eq.lhs, var) for var in xs] for eq in FAl])
DFAr = sp.Matrix([[sp.diff(eq.lhs, var) for var in xs] for eq in FAr])
DFAb = sp.Matrix([[sp.diff(eq.lhs, var) for var in xs] for eq in FAb])
DFAd = sp.Matrix([[sp.diff(eq.lhs, var) for var in xs] for eq in FAd])
DFAt = sp.Matrix([[sp.diff(eq.lhs, var) for var in xs] for eq in FAt])


#saving using pickle

with open('jacobians/FAr.pkl', 'wb') as file:
    pickle.dump(FAr, file)    
with open('jacobians/FAl.pkl', 'wb') as file:
    pickle.dump(FAl, file) 
with open('jacobians/FAb.pkl', 'wb') as file:
    pickle.dump(FAb, file) 
with open('jacobians/FAt.pkl', 'wb') as file:
    pickle.dump(FAt, file)  
with open('jacobians/FAd.pkl', 'wb') as file:
    pickle.dump(FAd, file)  

with open('jacobians/DFAr.pkl', 'wb') as file:
    pickle.dump(DFAr, file)
with open('jacobians/DFAl.pkl', 'wb') as file:
    pickle.dump(DFAl, file) 
with open('jacobians/DFAb.pkl', 'wb') as file:
    pickle.dump(DFAb, file) 
with open('jacobians/DFAt.pkl', 'wb') as file:
    pickle.dump(DFAt, file)  
with open('jacobians/DFAd.pkl', 'wb') as file:
    pickle.dump(DFAd, file)     