import sympy as sp
import pickle
# Definir las variables simbólicas
u, duds, dudss, dudtau ,g, dgds, dgdss, dgdtau, nu, alpha, t, s = sp.symbols('u duds duss dudtau g dgds dgdss dgdtau nu alpha t s')
#Mapeado
dudx=duds/dgds
dudxx=(dudss/dgds-duds*dgdss/dgds**2)/dgds
dgdx=dgds/dgds
dgdxx=(dgdss/dgds-dgds*dgdss/dgds**2)/dgds
dudt=dudtau-dgdtau/dgds*duds
#Concentrador
M=1/(alpha*dudx**2+0.4)
Ms=sp.diff(M, duds)*dudss+sp.diff(M, dgds)*dgdss
n=(dgds**2)**1.5

#vertor
# Crear un vector de variables
variables = [u, duds, dudss, dudtau ,g, dgds, dgdss, dgdtau]

# Definir la ecuación 1 bulk
FAb = [
sp.Eq(dudt + u * dudx-nu * dudxx,0),
sp.Eq(dgds*dgdss-n*Ms,0)
#sp.Eq(duds,0),
#sp.Eq(dgdss,0)
]





#Guardamos las cosas
with open('jacobians/FAb.pkl', 'wb') as archivo:
    pickle.dump(FAb, archivo)


# Calcular el jacobiano de la ecuación respecto al vector de variables
DFAb = sp.Matrix([[sp.diff(eq.lhs, var) for var in variables] for eq in FAb])
DFAb=DFAb*sp.sign(g**2+1)

#Guardamos las cosas


with open('jacobians/DFAb.pkl', 'wb') as archivo:
    pickle.dump(DFAb, archivo)
  

# Definir la ecuación left
FAl = [
sp.Eq(u-1,0),
sp.Eq(g,0)
]
# Calcular el jacobiano de la ecuación respecto al vector de variables
DFAl = sp.Matrix([[sp.diff(eq.lhs, var) for var in variables] for eq in FAl])
DFAl=DFAl*sp.sign(g**2+1)
#Guardamos las cosas
with open('jacobians\DFAl.pkl', 'wb') as archivo:
    pickle.dump(DFAl, archivo)
with open('jacobians\FAl.pkl', 'wb') as archivo:
    pickle.dump(FAl, archivo)    

print(FAl)

# Definir la ecuación r
FAr = [
sp.Eq(u+1,0),
sp.Eq(g-1,0)
]
# Calcular el jacobiano de la ecuación respecto al vector de variables
DFAr = sp.Matrix([[sp.diff(eq.lhs, var) for var in variables] for eq in FAr])
DFAr=DFAr*sp.sign(g**2+1)
#Guardamos las cosas
with open('jacobians/DFAr.pkl', 'wb') as archivo:
    pickle.dump(DFAr, archivo)
with open('jacobians/FAr.pkl', 'wb') as archivo:
    pickle.dump(FAr, archivo)    


