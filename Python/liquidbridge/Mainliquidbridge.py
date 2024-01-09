#importing from phyton...
import numpy as np
import pickle
import sympy as sp
from sympy import symbols
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse import coo_matrix #to fill a sparse matrix
from scipy.sparse import bmat #to create sparse block matrices
from scipy.sparse.linalg import spsolve #to solve linear systems.
from scipy.sparse import spdiags, csc_matrix
from scipy.sparse import diags
from scipy.sparse import kron, eye
import scipy.sparse as sps
from scipy.sparse import hstack, vstack
from scipy.optimize import fsolve
import time
import warnings
#importing from my files...
from colocationmatrices import colocation #all the collocation matrices
import auxilarfunctions
# Ignora las advertencias de tiempo de ejecución
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Reading symbolics Jacobians

with open('jacobians/FAb.pkl', 'rb') as archivo:
 FAb = pickle.load(archivo)
with open('jacobians/FAl.pkl', 'rb') as archivo:
 FAl = pickle.load(archivo)
with open('jacobians/FAr.pkl', 'rb') as archivo:
 FAr = pickle.load(archivo) 
with open('jacobians/FAt.pkl', 'rb') as archivo:
 FAt = pickle.load(archivo) 
with open('jacobians/FAd.pkl', 'rb') as archivo:
 FAd = pickle.load(archivo) 
# Obtain the left-hand expressions of the equations
FAb = [eq.lhs for eq in FAb]
FAl = [eq.lhs for eq in FAl]
FAr = [eq.lhs for eq in FAr]
FAt = [eq.lhs for eq in FAt]
FAd = [eq.lhs for eq in FAd]

with open('jacobians/DFAb.pkl', 'rb') as archivo:
 DFAb = pickle.load(archivo)
with open('jacobians/DFAl.pkl', 'rb') as archivo:
 DFAl = pickle.load(archivo)
with open('jacobians/DFAr.pkl', 'rb') as archivo:
 DFAr = pickle.load(archivo)
with open('jacobians/DFAt.pkl', 'rb') as archivo:
 DFAt = pickle.load(archivo)
with open('jacobians/DFAd.pkl', 'rb') as archivo:
 DFAd = pickle.load(archivo)
#Making the function faster using lambdify
 

 # Define symbolic variables
w, dwds, dwdss, dwdy, dwdyy, dwdys, dwdt0 = sp.symbols('w dwds dwdss dwdy dwdyy dwdys dwdt0')
u, duds, dudss, dudy, dudyy, dudys, dudt0 = sp.symbols('u duds dudss dudy dudyy dudys dudt0')
p, dpds, dpdss, dpdy, dpdyy, dpdys, dpdt0 = sp.symbols('p dpds dpdss dpdy dpdyy dpdys dpdt0')
f, dfds, dfdss, dfdy, dfdyy, dfdys, dfdt0 = sp.symbols('f dfds dfdss dfdy dfdyy dfdys dfdt0')
# Define  symbolicspace (s,y) and symbolic time (t0)
s, y, t0 = sp.symbols('s y t0')
# Define space (s,y) and time (tau)
Cmu, Bo = sp.symbols('Cmu  Bo') 
xs = [w, dwds, dwdss, dwdy, dwdyy, dwdys, dwdt0] 
xs = xs + [u, duds, dudss, dudy, dudyy, dudys, dudt0]
xs = xs + [p, dpds, dpdss, dpdy, dpdyy, dpdys, dpdt0]
xs = xs + [f, dfds, dfdss, dfdy, dfdyy, dfdys, dfdt0]
xs = xs+[s,y]
xs = xs+[Cmu,Bo]


FAb_func = [sp.lambdify(xs, expr, 'numpy') for expr in FAb]
FAl_func = [sp.lambdify(xs, expr, 'numpy') for expr in FAl]
FAr_func = [sp.lambdify(xs, expr, 'numpy') for expr in FAr]
FAd_func = [sp.lambdify(xs, expr, 'numpy') for expr in FAd]
FAt_func = [sp.lambdify(xs, expr, 'numpy') for expr in FAt]

DFAb_func = [sp.lambdify(xs, expr, 'numpy') for expr in DFAb]
DFAr_func = [sp.lambdify(xs, expr, 'numpy') for expr in DFAr]
DFAl_func = [sp.lambdify(xs, expr, 'numpy') for expr in DFAl]
DFAt_func = [sp.lambdify(xs, expr, 'numpy') for expr in DFAt]
DFAd_func = [sp.lambdify(xs, expr, 'numpy') for expr in DFAd]


#End of the symbolic part
#Parameters (figure 4 Herrada&Montanero 2016)
Cmu=0.0349 #Capilar number
Bo=0.478 #Bond number
Lambda=1.72 #Slender
V=0.925    #Dimensionless volumem


# number of derivatives
nd=int('7')
NDA=nd
#number of variables
nv=int('4')
NVA=nv
# spatial discretization
ns =int(61) #in the s direction (\xi in the paper)
ny =int(9)  #in the y direction (\eta in the paper)
#total number of spatial points
Np=ns*ny
#y-collocation matrices (spectral)
y, dy, dyy =colocation.Chevigood(ny-1, 1.0,0.)
#s-collocation matrices (4th finite differences)
s, ds, dss =colocation.finites4th(ns, 2.*Lambda)
#s, ds, dss =collocacion.Chevigood(ns-1, 2.*Lambda,0)
#s-collocation matrices (4th finite differences)

#getting the full yxs collocation matrices
ss, yy, dds, ddss, ddy,ddyy,ddsy =auxilarfunctions.matrix_A(s,ds,dss,ns,y,dy,dyy,ny)

#getting the pointers to the Boundaries
l_l, l_r, l_d, l_t, l_b= auxilarfunctions.computingpointers(ny,ns)





#Parameters (figure 4 Herrada&Montanero 2016)
Cmu=0.0349
Bo=0.478
Lambda=1.72
V=0.925

#initial conditions
w = np.zeros((ny, ns)) #axial velocity
u = np.zeros((ny, ns)) #radial velocity
p = np.ones((ny, ns))  #radial velocity

# Getting the initial shape (r = f_i(s)) by solving the Laplace-Young equation
f = np.ones(ns)

pref = 1
x0 = np.concatenate([f, [pref]])
#solving initial shape
# Define the tolerance and maximum iterations
xtol = 1e-8
maxfev = 1000

# Set the options directly
options = {'xtol': xtol, 'maxfev': maxfev}
    
#getting the function to solve the shape   
Fshape=auxilarfunctions.Fshape
# Using fsolve to solve the equation
x = fsolve(Fshape, x0, args=(ns, s, ds, dss, Bo, V, Lambda), xtol=xtol, maxfev=maxfev)
f = x[:ns]
pref = x[ns]



# repmat(f, [ny, 1])
f = np.tile(f, (ny, 1))
# p = pref * ones(ny, ns)
p = np.ones((ny, ns)) * pref
#p=ss*yy
#w=ss
#u=yy
#f=1+ss*yy
yy=np.reshape(yy, Np)    
ss=np.reshape(ss, Np)  # vector [2N, 1]
# We construct a unique vector of unknowns x0 (guess solution) [STEP 7]
x0 = np.concatenate([
    np.reshape(w, Np),  # vector [2N, 1]
    np.reshape(u, Np),
    np.reshape(p, Np),
    np.reshape(f, Np)
])
#previous time solutions
x0m=x0
x0mm=x0
# Time steps
dt = 0.1
dt1 = 0.1
#introducing perturbation in Bo to  break the liquid bridge
Bo=Bo+0.085*Bo
#End of the simulation
t_end=130
# Start timer
timeninit = time.time()
NT=int(t_end/dt)+1

#array for time
timef=np.zeros((NT))
f0=np.zeros((NT))
for ll in range(NT):
 # Time derivative 
 dt2=dt1+dt; 
 #second order backwar differencesd
 bm= -(dt2/dt)/(dt2-dt); 
 bmm= (dt/dt2)/(dt2-dt); 
 bp=-((dt2/dt)**2-1)/((dt2/dt)*(dt-dt2));  #[STEP 6]

 #computing time to check speed

 


 errorf=10**10

 while errorf>0.0001:
   # Pointers variables
   lw = list(range(Np))
   lu = list(range(Np, 2 * Np))
   lp = list(range(2 * Np, 3 * Np))
   lf = list(range(3 * Np, 4 * Np))

   #time derivative
   x0t = bp*x0 + bm*x0m + bmm*x0mm
   #computing all the derivaties
   w=x0[lw]
   dwds=dds@x0[lw]
   dwdss=ddss@x0[lw]
   dwdy=ddy@x0[lw]
   dwdyy=ddyy@x0[lw]
   dwdsy=ddsy@x0[lw]
   dwdt0=x0t[lw]
   u=x0[lu]
   duds=dds@x0[lu]
   dudss=ddss@x0[lu]
   dudy=ddy@x0[lu]
   dudyy=ddyy@x0[lu]
   dudsy=ddsy@x0[lu]
   dudt0=x0t[lu]
   p=x0[lp]
   dpds=dds@x0[lp]
   dpdss=ddss@x0[lp]
   dpdy=ddy@x0[lp]
   dpdyy=ddyy@x0[lp]
   dpdsy=ddsy@x0[lp]
   dpdt0=x0t[lp]
   f=x0[lf]
   dfds=dds@x0[lf]
   dfdss=ddss@x0[lf]
   dfdy=ddy@x0[lf]
   dfdyy=ddyy@x0[lf]
   dfdsy=ddsy@x0[lf]
   dfdt0=x0t[lf]
   a=np.reshape(dwds,[ny,ns])
   b=np.reshape(dudy,[ny,ns])

  
   # Create the symbolic vector using list concatenation
   xs = [w, dwds, dwdss, dwdy, dwdyy, dwdsy, dwdt0] 
   xs = xs + [u, duds, dudss, dudy, dudyy, dudsy, dudt0]
   xs = xs + [p, dpds, dpdss, dpdy, dpdyy, dpdsy, dpdt0]
   xs = xs + [f, dfds, dfdss, dfdy, dfdyy, dfdsy, dfdt0]
   xs = xs+[ss,yy]
   #xs = xs+[Cmu,Bo]
 
   #Create the list at the down boundary
   xs_d = [param[l_d] for param in xs]
   #Create the list at the  top boundary
   xs_t = [param[l_t] for param in xs]
   #Create the list at the left boundary
   xs_l = [param[l_l] for param in xs]
   #Create the list at the right boundary
   xs_r = [param[l_r] for param in xs]
   #Create the list at the right boundary
   xs_b = [param[l_b] for param in xs]

 
   #Computing analytical functions
   #bulk
   result_b = [func(*xs_b,Cmu,Bo) for func in FAb_func]
   #down
   result_d = [func(*xs_d,Cmu,Bo) for func in FAd_func]
   #top
   result_t = [func(*xs_t,Cmu,Bo) for func in FAt_func] 
   #left
   result_l = [func(*xs_l,Cmu,Bo) for func in FAl_func]
   #right
   result_r = [func(*xs_r,Cmu,Bo) for func in FAr_func] 
   FAA = np.zeros((nv,Np))
   FAA[:,l_b]= np.array(result_b)
   FAA[:,l_d]=np.array(result_d)
   FAA[:,l_t]=np.array(result_t)
   FAA[:,l_l]=np.array(result_l)
   FAA[:,l_r]=np.array(result_r)


   #Computing analytical Jacobians
   #bulk
   result = [func(*xs,Cmu,Bo) for func in DFAb_func]
   result_b = [np.full(Np, x) if isinstance(x, (int, float)) else x for x in result]
   #down
   result = [func(*xs_d,Cmu,Bo) for func in DFAd_func]
   result_d = [np.full(ns, x) if isinstance(x, (int, float)) else x for x in result]
   #top
   result = [func(*xs_t,Cmu,Bo) for func in DFAt_func]
   result_t = [np.full(ns, x) if isinstance(x, (int, float)) else x for x in result]
   #left
   result = [func(*xs_l,Cmu,Bo) for func in DFAl_func]
   result_l = [np.full(ny, x) if isinstance(x, (int, float)) else x for x in result]
   #right
   result = [func(*xs_r,Cmu,Bo) for func in DFAr_func] #reshape(2,8)
   result_r = [np.full(ny, x) if isinstance(x, (int, float)) else x for x in result]
   #saving in array
   DFAA=np.zeros((nv*nv*nd,Np))
   
   DFAA[:,:]=np.array(result_b)
   #BC
   DFAA[:,l_d]=np.array(result_d)
   DFAA[:,l_t]=np.array(result_t)
   DFAA[:,l_l]=np.array(result_l)
   DFAA[:,l_r]=np.array(result_r)
   #reshaping
   DFAA=DFAA.reshape(nv,nv*nd,Np)
   

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
         c4=C[4,:].reshape(1,-1)
         c5=C[5,:].reshape(1,-1)
         c6=C[6,:].reshape(1,-1)
         B = spdiags(c0, 0, Np, Np)+spdiags(c1, 0, Np, Np)@dds+spdiags(c2, 0, Np, Np)@ddss+spdiags(c3, 0, Np, Np)@ddy+spdiags(c4, 0, Np, Np)@ddyy+spdiags(c5, 0, Np, Np)@ddsy+spdiags(c6, 0, Np, Np)*bp
         ablock[j][k]=B
       


   b= np.concatenate(xablock)
   # Convert sparse array lists to sparse arrays
   a = vstack([hstack([m for m in row if m is not None]) for row in ablock if any(m is not None for m in row)])
   #invert matrix
   dx = spsolve(a, b)
   #updating solution
   x0=x0+dx
   # print error
   errorf=np.linalg.norm(dx)
   print(errorf)
 #end of the Newton method
 #making variable readble  
 w=np.reshape(x0[lw], [ny,ns])
 u=np.reshape(x0[lu], [ny,ns])
 p=np.reshape(x0[lp], [ny,ns])
 f=np.reshape(x0[lf], [ny,ns])

 timef[ll]=ll*dt; #time
 #computing the shape at z=Lambda/2
 aa=f[ny-1, int((ns-1)/4)]
 f0[ll]=aa
 print(timef[ll])
 print(f0[ll])
 x0mm=x0m
 x0m=x0


timeend = time.time()
# Calculate the elapsed time in seconds
elapsed_time = timeend - timeninit
print(f"Elapsed_ time: {elapsed_time} seconds")


#plotting figure 4 Herrada&Montanero (2016)
plt.plot(timef, f0)
plt.xlabel('t')  # Replace with your actual x-axis label
plt.ylabel('$R(\\frac{\\Lambda}{2}, t)$')
 #plt.title('Your Plot Title')  # Replace with your actual plot title
plt.show() 
 # Stop timer
