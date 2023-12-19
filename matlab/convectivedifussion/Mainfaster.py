#importing from phyton...
import numpy as np
import pickle
import sympy as sp
from sympy import symbols
from scipy.sparse import lil_matrix
from scipy.sparse import coo_matrix #to fill a sparse matrix
from scipy.sparse import bmat #to create sparse block matrices
from scipy.sparse.linalg import spsolve #to solve linear systems.
from scipy.sparse import spdiags, csc_matrix
from scipy.sparse import diags
import scipy.sparse as sps
from scipy.sparse import hstack, vstack

#importing from my files...
from matricescolocacion import collocacion #all the collocation matrices
# number of derivatives
nd=int('4')
NDA=nd
#number of variables
nv=int('2')
NVA=nv
# spatial discretization
ns =int(40000)
ntA=ns
H = 1.0
#
s0, ds0, dss0 =collocacion.finitas4(ns, H)
# Convert the row vector to a column vector
#s0 = s0.reshape(-1, 1)ns


# Define symbolic variables
u, duds, dudss, dudtau ,g, dgds, dgdss, dgdtau, nu, alpha, t, s = sp.symbols('u duds duss dudtau g dgds dgdss dgdtau nu alpha t s')
#we charge the analytical Jacobians
with open('jacobianos/FAb.pkl', 'rb') as archivo:
 FAb = pickle.load(archivo)
with open('jacobianos/FAl.pkl', 'rb') as archivo:
 FAl = pickle.load(archivo)
with open('jacobianos/FAr.pkl', 'rb') as archivo:
 FAr = pickle.load(archivo) 
# Obtain the left-hand expressions of the equations
FAb = [eq.lhs for eq in FAb]
FAl = [eq.lhs for eq in FAl]
FAr = [eq.lhs for eq in FAr]

with open('jacobianos/DFAb.pkl', 'rb') as archivo:
 DFAb = pickle.load(archivo)
with open('jacobianos/DFAl.pkl', 'rb') as archivo:
 DFAl = pickle.load(archivo)
with open('jacobianos/DFAr.pkl', 'rb') as archivo:
 DFAr = pickle.load(archivo)
 
#print([FAb])
#Making the function faster using lambdify
FAb_func = [sp.lambdify((u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau, nu, alpha, t, s), expr, 'numpy') for expr in FAb]
FAl_func = [sp.lambdify((u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau, nu, alpha, t, s), expr, 'numpy') for expr in FAl]
FAr_func = [sp.lambdify((u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau, nu, alpha, t, s), expr, 'numpy') for expr in FAr]

DFAb_func = [sp.lambdify((u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau, nu, alpha, t, s), expr, 'numpy') for expr in DFAb]
DFAr_func = [sp.lambdify((u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau, nu, alpha, t, s), expr, 'numpy') for expr in DFAr]
DFAl_func = [sp.lambdify((u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau, nu, alpha, t, s), expr, 'numpy') for expr in DFAl]




#We define the numerical variables
u0=1-2*s0
duds0=ds0@u0
dudss0=dss0@u0
dudtau0=0*u0
g0=s0
dgds0=ds0@g0
dgdss0=dss0@g0
dgdtau0=0*g0

#time step
dt=10**10
dt1=10**10
# Time derivative 
dt2=dt1+dt; 
 #second order backwar differencesd
bm= -(dt2/dt)/(dt2-dt); 
bmm= (dt/dt2)/(dt2-dt); 
bp=-((dt2/dt)**2-1)/((dt2/dt)*(dt-dt2));  #[STEP 6]

#define numerical parameters
nu0=0.1
alpha0=0.01

import time

# Start timer
timeninit = time.time()


errorf=10**10

while errorf>0.01:
 result = [func(u0, duds0, dudss0, dudtau0, g0, dgds0, dgdss0, dgdtau0, nu0, alpha0, 0, s0) for func in FAb_func]
 result0 = [func(u0[0], duds0[0], dudss0[0], dudtau0[0], g0[0], dgds0[0], dgdss0[0], dgdtau0[0], nu0, alpha0, 0, s0[0]) for func in FAl_func]
 result1 = [func(u0[ns-1], duds0[ns-1], dudss0[ns-1], dudtau0[ns-1], g0[ns-1], dgds0[ns-1], dgdss0[ns-1], dgdtau0[ns-1], nu0, alpha0, 0, s0[ns-1]) for func in FAr_func]
 #reshape(2,8)
 FAA = np.array(result)
 FAA[:,0]=np.array(result0)
 FAA[:,ns-1]=np.array(result1)
 result=[func(u0, duds0, dudss0, dudtau0, g0, dgds0, dgdss0, dgdtau0, nu0, alpha0, 0, s0) for func in DFAb_func]
 resultf = [np.zeros(ns) if isinstance(x, int) and x == 0 else x for x in result]
 result0 = [func(u0[0], duds0[0], dudss0[0], dudtau0[0], g0[0], dgds0[0], dgdss0[0], dgdtau0[0], nu0, alpha0, 0, s0[0]) for func in DFAl_func]
 result1 = [func(u0[ns-1], duds0[ns-1], dudss0[ns-1], dudtau0[ns-1], g0[ns-1], dgds0[ns-1], dgdss0[ns-1], dgdtau0[ns-1], nu0, alpha0, 0, s0[ns-1]) for func in DFAr_func]
 #cambiamos 0 por arrays
 DFAA=np.array(resultf)
 DFAA[:,0]=np.array(result0)
 DFAA[:,ns-1]=np.array(result0)

 DFAA=DFAA.reshape(nv,nv*nd,ns)
 
  # Getting numerical Jacobian (step 9)
 ablock = [[None] * NVA for _ in range(NVA)]
 xablock = [None] * NVA

 for j in range(NVA):
    C1 = FAA[j, :]
    xablock[j] = -C1
    
    for k in range(NVA):
        
        km = k * NDA
        kp = (k + 1) * NDA
        C = DFAA[j, range(km,kp), :]
       
        c0=C[0,:].reshape(1,-1)
        c1=C[1,:].reshape(1,-1)
        c2=C[2,:].reshape(1,-1)
        c3=C[3,:].reshape(1,-1)
        B = spdiags(c0, 0, ntA, ntA)+spdiags(c1, 0, ntA, ntA).dot(ds0)+spdiags(c2, 0, ntA, ntA).dot(dss0)+bp*spdiags(c3, 0, ntA, ntA)
        ablock[j][k]=B
 b= np.concatenate(xablock)
 # Convert sparse array lists to sparse arrays
 a = vstack([hstack([m for m in row if m is not None]) for row in ablock if any(m is not None for m in row)])
#invert matrix
 dx = spsolve(a, b)

 u0=u0+dx[0:ns]
 g0=g0+dx[ns:2*ns]
 duds0=ds0@u0
 dudss0=dss0@u0
 dudtau0=0*u0
 dgds0=ds0@g0
 dgdss0=dss0@g0
 dgdtau0=0*g0
# print(duds0)
 errorf=np.linalg.norm(dx)
 print(errorf)

 

# Stop timer
timeend = time.time()
# Calcular el tiempo transcurrido en segundos
tiempo_transcurrido = timeend - timeninit
print(f"Tiempo transcurrido: {tiempo_transcurrido} segundos")
