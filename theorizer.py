from equationSystem import *
from equation import *
from m2_functions import *
from unitOfMeasure import *
from derivative import *

    

if __name__ == "__main__":
    u = BaseUnit('u')
    U = UofM(u)
    UP = UofM(Pow(U._units, 2))
    d1, d2, m1, m2, W, p, Fc, Fg = variables('d1,d2,m1,m2,W,p,Fc,Fg')
    #d1,d2,m1,m2,w,p,Fc,Fg = symbols('d1,d2,m1,m2,w,p,Fc,Fg')   #see if it is possible to do this using variables()
    vars = [Fc,Fg,W,d1,d2,m1,m2,p]

    eqnSystem = EquationSystem.GenerateRandom(vars=vars, measuredVars=vars, numEqns=4, maxVarsPerEqn=4)
    res = eqnSystem._equations[0].isDimensionallyConsistent()
    #measured_vars = [d1,d2,m1,m2,p]
    eq2 = Equation(d1**2*Fg + 2*d1*d2*Fg + d2**2*Fg - m1*m2)
    res = eq2.isDimensionallyConsistent()
    eq2prime = Equation(d1 ** 2 * Fg + 2 * d1 * d2 * Fg + d2 ** 2 * Fg )
    res = eq2prime.isDimensionallyConsistent()
    vbles = eq2._variables
    terms = eq2.getTerms()
    m,s = base_units('m,s') #see below
    units = UofM(m/s)  #should not be done in terms of symbols - perhaps a new class called units, which would be veneer over Symbols

    v = Variable('v')
    dxdt = Derivative('dxdt')
    eq = Equation(v - dxdt)
    res = eq.isDimensionallyConsistent()


    F,m,d = variables('F,m,d')
    dxdt, d2xdt2 = derivatives('dxdt,d2xdt2')
    eq3 = Equation(F - m*d2xdt2)
    res = eq3.isDimensionallyConsistent()

    #eq3 = Polynomial(Fc - m2*d2*w**2)
    #eq4 = Polynomial(Fg - Fc)
    #eq5 = Polynomial(w*p - 1)
    #eqns = [eq2, eq3, eq4, eq5]
    #eqnSystem = EquationSystem(vars=vars, measuredVars=measured_vars, equations=eqns)

    #print("System: \n" + eqnSystem.toString() + "\n")

    #isConsistent(eqnSystem)
    #print(checkConsistencyReplacedSystems(eqnSystem, 3, 10))
    #print(checkConsistencyRandomSystems(vars, 10))
    #project(eqnSystem)

    vars = [Fc,Fg,w,d1,d2,m1,m2,p]
    # Take the measured vars to be a random subset of the variables
    #measured_vars = random.sample(vars, len(vars)//2)
    
    measured_vars = vars.copy()
    random.shuffle(measured_vars)

    EquationSystem.ProjectRandomSystems(vars, measured_vars, 10)














