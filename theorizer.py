from equationSystem import *
from polynomial import *

a,b,c,d,e,f,g = symbols('a,b,c,d,e,f,g')
vars = [a,b,c,d,e,f,g]
eqnSystem = EquationSystem.GenerateRandom(vars=vars, numEqns=4, maxVarsPerEqn=4)

print("System: \n" + eqnSystem.toString())

eqn = eqnSystem.replaceRandomEqnByIndex(eqnIndex=2)
print("Replacement for 3rd eqn: " + eqn.toString())

print("\n")

d1,d2,m1,m2,w,p,Fc,Fg = symbols('d1,d2,m1,m2,w,p,Fc,Fg')
vars = [d1,d2,m1,m2,w,p,Fc,Fg]
eq1 = Polynomial(d1*m1 - d2*m2)
eq2 = Polynomial(d1**2*Fg + 2*d1*d2*Fg + d2**2*Fg - m1*m2)
eq3 = Polynomial(Fc - m2*d2*w**2)
eq4 = Polynomial(Fg - Fc)
eq5 = Polynomial(w*p - 1)
eqns = [eq1, eq2, eq3, eq4, eq5]
#eqns = [eq1, eq4]
eqnSystem = EquationSystem(equations=eqns)

print("System: \n" + eqnSystem.toString())

eqn = eqnSystem.replaceRamdomEqn()

print("Revised System: \n" + eqnSystem.toString())

e,m,c = symbols('e,m,c')
eq = Polynomial(e/m*c**2)  #fractional polys are allowed!
print(eq._poly)